# -*- coding: utf-8 -*-
"""
/***************************************************************************
 IRIS - RoutingCinemat
                                 A QGIS plugin
 Indici di RIschio Sismico (per le dighe in calcestruzzo)
                              -------------------
        begin                : 2018-01-04
        git sha              : $Format:%H$
        copyright            : (C) 2018 by L. Mancusi /RSE S.p.A
        email                : leonardo.mancusi@rse-web.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

The script, assigned a DamID calculates the propagation of the flood wave
using kinematic model and propagation of steep waves theories

Reads the coefficients and exponents of the geometrical and hydraulic characteristics
of the various river reaches from the file:
    - MatriceAexp.csv

Calculate the flood hydrogram at the first upstream section

Starting from the first section, it calculates the propagation for each time interval
of the instantaneous flow rate in the reach and the point where it reaches the front.
Calculates the output hydrogram from the reach and uses it as input for the next reach.

Makes the graph of the maximum flow rate and water depth along the progressive distance

Save results to file:
    - Q_H_max.csv

"""
import os, sys
import numpy
import math

try:
    from pyspatialite import dbapi2 as db
    import sqlite3
except:
    import sqlite3

try:
    from osgeo import ogr
    from osgeo.osr import osr
except:
    import ogr
    import osr

import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib.font_manager import FontProperties
import csv

def DamBreakHyrograph(VolumeInvaso,AltezzaDiga,LarghezzaSezioneRottura,NomeDiga):
    """
    Evaluation of the dambreak hydrograph:
    -------------------------------------
    VolumeInvaso            : Reservoir volume in millions of cubic meters
    AltezzaDiga             : height of the dam in meters
    LarghezzaSezioneRottura : Width of the breach in meters
    NomeDiga                : name of the dam
    """
    # form Mmc to mc
    VolumeInvaso=VolumeInvaso*10**6

    # reservoir volume curve : V=k*h^alfa
    # ====================================
    alfa=1.5
    kinvaso=VolumeInvaso/(AltezzaDiga**alfa)

    Vcurr=VolumeInvaso
    Vfin=VolumeInvaso*0.001

    # flow calculation for instant collapse
	# ====================================
	# gravity acceleration
    G=9.81

    Yc=4.0/9.0*AltezzaDiga

    Qmax=8.0/27.0*LarghezzaSezioneRottura*AltezzaDiga*(G*AltezzaDiga)**0.5

    Vmax=Qmax/(Yc*LarghezzaSezioneRottura)

    # Total theoretical energy height in the hypothesis of permanent motion
    Htot=Yc+Vmax**2/2./G
    # volume
    Vol_Htot=kinvaso*Htot**alfa

    dt=60.0
    Idrogramma=[]
    Time=[]
    T=0.0
    Time.append(T)
    Idrogramma.append(Qmax)

    while Vcurr>=Vol_Htot:
        dv=-Qmax*dt
        # update volume
        Vcurr+=dv
        if Vcurr<0.0:
            while Vcurr<0.0:
                dt=dt/2.0
                Vcurr+=-dv
                dv=-Qmax*dt
                Vcurr+=dv
        Idrogramma.append(Qmax)
        T+=dt
        Time.append(T)

    alfa_1=1.0/alfa
    # water depth
    hcurr=math.pow((Vcurr/kinvaso),alfa_1)

    hmin=0.5

    while hcurr>hmin:
        # critical depth for rectangular section
        Yc=2.0/3.0*hcurr
        # critical speed for rectangular section
        Uc=math.sqrt(G*Yc)
        # critical flow rate for rectangular section
        Qc=Uc*Yc*LarghezzaSezioneRottura
        # volume released in time unit
        dv=-Qc*dt
        # update the volume at the end of dt
        Vcurr+=dv
        if Vcurr<0.0:
            while Vcurr<0.0:
                dt=dt/2.0
                Vcurr+=-dv
                dv=-Qc*dt
                Vcurr+=dv
        Idrogramma.append(Qc)
        T+=dt
        Time.append(T)
        # update the water depth at the end of dt
        hcurr=math.pow((Vcurr/kinvaso),alfa_1)

    # calculate hydrogram volume
    Vol=0.0
    nn=len(Idrogramma)
    for i in range(1,nn):
        Vol+=(Time[i]-Time[i-1])*(Idrogramma[i]+Idrogramma[i-1])/2.0

    grafico=0
    if grafico>0:
        try:
            plt.subplot(1, 1, 1)

            plt.plot(Time,Idrogramma)

            plt.grid(True)
            txt='Dam: %s  - Dam-Break discharge' % NomeDiga
            plt.title(txt)
            plt.xlabel('Time (sec)')
            plt.ylabel('discharge (cm/s)')
            plt.legend()
            plt.show()
        except:
            pass

    return Idrogramma, Time


def ValoreInX(x,n,datix,datiy):

    nx=len(datix)
    ny=len(datiy)
    if nx==ny:
        if x>datix[0]:
            if x<datix[-1]:
                for i in range(1,nx):
                    if datix[i]>=x:
                        y=datiy[i-1]+(datiy[i]-datiy[i-1])*(x-datix[i-1])/(datix[i]-datix[i-1])
                        break
            else:
                y=datiy[-1]
        else:
            y=datiy[0]
    else:
        sys.exit()

    return y

def RunCinemat(mydb_path_user,DamID,PathFiles,grafico2):

    """
    Run Cinemat model for routing flood downstream
    """

    pythonver_log=sys.version
    ppp=pythonver_log.split(' ')
    py_ver=ppp[0][0:1]


    NotErr=bool('True')
    errMsg='OK'


    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    FileOutput=PathFiles+os.sep+'Q_H_max.csv'


    # =========================================
    # apertura del database sqlite
    # =========================================

    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
   # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    NomeTabella='DAMS'

    sql='SELECT ResVolMcm,Height_m,BreachWidth,Name FROM %s WHERE DamID=%d' % (NomeTabella,DamID)
    cur.execute(sql)
    DatiDiga=cur.fetchone()

    if DatiDiga!=None:

        VolumeInvaso=float(DatiDiga[0])

        # height of the dam
        AltezzaDiga = float(DatiDiga[1])

        # Width of the breach
        LarghezzaSezioneRottura = float(DatiDiga[2])

        # Name of the dam
        NomeDiga=DatiDiga[3]

    else:
        errMsg = 'In the table= %s there is no data for the dam num =%s \nPerform data saving first !' % (NomeTabellaTest,DamID)
        NotErr= bool()
        return NotErr, errMsg

    # reading the data for the calculation
    NomeTabellaMatrice='MatriceAexp'

    sql='SELECT '
    sql+=' PixDist'
    sql+=', progr_lungo_fiume'
    sql+=', distanza_fiume'
    sql+=', distanza_linea_retta'
    sql+=', pend'
    sql+=', ka'
    sql+=', ma'
    sql+=', kq'
    sql+=', mq'
    sql+=', kcel'
    sql+=', mcel'
    sql+=' FROM %s' % NomeTabellaMatrice
    sql+=' WHERE DamID=%d' % (DamID)
    sql+=' ORDER BY PixDist;'
    cur.execute(sql)
    MatriceDati=cur.fetchall()

    ListaTratti=[]
    Progressive_fiume=[]
    Distanze_fiume=[]
    # coefficient matrix: ka;ma;kq;mq;kcel;mcel for every reach
    MatriceCoeff=[]
    for row in MatriceDati:
        ListaTratti.append(int(row[0]))
        Progressive_fiume.append(float(row[1]))
        Distanze_fiume.append(float(row[2]))
        MatriceCoeff.append(row[3:])


    # Calculation of the flood hydrograph
    # ................................

    Idrogramma, T = DamBreakHyrograph(VolumeInvaso,AltezzaDiga,LarghezzaSezioneRottura,NomeDiga)


    # the first input corresponds to the hydrograph at the initial section
    Qin=numpy.array(Idrogramma,dtype =numpy.float)
    maxQ=Qin.max()

    nn=len(MatriceCoeff)
    x=0.0

    numtime=len(Qin)

    # lists of times, distances and the envelope of the maximum flow along the river
    TQX=[]
    XQX=[]
    QX=[]
    # Water depth
    HX=[]
    # water width
    BX=[]
    # water velocity
    VX=[]

    TrattoX=[]
    PendX=[]
    TQX.append(T[0])
    XQX.append(0.0)
    Progr_fiume=[]
    Progr_fiume.append(0.0)

    QX.append(Qin[0])
    TrattoX.append(0)
    PendX.append(MatriceCoeff[0][1])

    # coeff. and exponent of the formula monomer of the flow Q=kq*h^mq
    kq=float(MatriceCoeff[0][4])
    mq=float(MatriceCoeff[0][5])
    invmq=1.0/mq
    h_0=math.pow(Qin[0]/kq,invmq)
    HX.append(h_0)

    # water width = ka*ma*h^(ma-1)
    # ----------------------------------------------------
    ka=float(MatriceCoeff[0][2])
    ma=float(MatriceCoeff[0][3])

    # water width calculation
    mb=ma-1.0
    b_0=ka*ma*math.pow(h_0,mb)
    BX.append(b_0)

    # velocity
    # ------------------
    A_0=ka*math.pow(h_0,ma)
    if A_0>0:
        V_0=Qin[0]/A_0
    else:
        V_0=0.0

    VX.append(V_0)


    # settings for the chart
    #============================
    fontP = FontProperties()
    fontP.set_size('small')

    Progr_fiume=[]
    Progr_fiume.append(0.0)


    # cycle along the river reaches
    # -------------------------------
    for i in range(nn):

        # abscissa and instants of the initial point of the previous reach
        xini=XQX[-1]
        xini_fiume=Progr_fiume[-1]

        tini=TQX[-1]
        # abscissa from the beginning of reach
        x=0.0
        xfronte=0.0

        # kx and mx coefficients of geometric and hydraulic reach parameters
        tratto=ListaTratti[i]

        # distance in a straight line: used for calculating the propagation
        distanza=float(MatriceCoeff[i][0])

        # reading distance on the river
        Progr_curr=Progressive_fiume[i]
        Parz_curr=Distanze_fiume[i]
        # relationship between the distance along the river axis and the distance in a straight line
        Amplificazione=Parz_curr/distanza

        # slope
        pend=float(MatriceCoeff[i][1])
        # coeff. and exponent of the area monomular formula
        ka=float(MatriceCoeff[i][2])
        ma=float(MatriceCoeff[i][3])
        # coeff. and exponent of the monomer flow rate formula Q=kq*h^mq
        kq=float(MatriceCoeff[i][4])
        mq=float(MatriceCoeff[i][5])
        # coeff. and exponent of the monomy formula of celerity
        kcel=float(MatriceCoeff[i][6])
        mcel=float(MatriceCoeff[i][7])

        # parameters for calculating the water width B=dA/dh
        # ----------------------------------------------------------------
        # B=ka*ma*h^(m-1) = kb*h^mb
        kb=ka*ma
        mb=ma-1.0

        # calculation of the  Dt
        # ......................
        t=T[0]
        tau=T[0]
        # calculate the corresponding water depth
        invmq=1.0/mq

        for k in range(1,numtime):

            dtau=float(T[k]-T[k-1])

            Qtau=Qin[k-1]
            Qtau1=Qin[k]

            h_tau=math.pow(Qtau/kq,invmq)
            h_tau1=math.pow(Qtau1/kq,invmq)
            #  dQ/dA for tau cost = celerity
            Qprimo=kcel*math.pow(h_tau,mcel)
            # dQ/dA^2 for tau cost = derivative of celerity
            mcel2=mcel-1.0
            Qseconda=kcel*mcel*math.pow(h_tau,mcel2)
            # computation A(tau)
            Atau=ka*math.pow(h_tau,ma)
            Atau1=ka*math.pow(h_tau1,ma)
            # computation dA/dtau
            dA=(Atau1-Atau)/dtau
            # Vfronte calculation that on dry riverbed (ie Qo = 0) is Vfronte = Qtau / Atau
            Vfronte=Qtau/Atau
            # then dt calculation
            numerat=Qseconda*dA*x/Qprimo-Qprimo
            denom=-Qprimo+Vfronte
            dt=numerat/denom*dtau
            # dx calculation
            dx=dt*Vfronte

            # check the distance traveled
            if (x+dx)>distanza:
                # calculates the times to reach the end of the section
                dx1=distanza-x
                dt1=dx1/Vfronte
                dtau1=denom/numerat*dt1
                t1=t+dt1
                tau1=tau+dtau1
                # flow at instant tau1
                Q1=ValoreInX(tau1,numtime,T,Qin)
                h_1=math.pow(Q1/kq,invmq)
                # calculating the water width
                b_1=kb*math.pow(h_1,mb)

                # calculating celerity
                c1=kcel*math.pow(h_1,mcel)

                # calculates the hydrograph leaving the reach
                Tout=[]
                Qout=[]
                Tout.append(t1)
                Qout.append(Q1)

                # inserts a new point in the TQX curve
                TQX.append(t1)
                x1=xini+distanza
                XQX.append(x1)

                # progressive distance upgrade along the river axis
                x1_fiume=xini_fiume+distanza*Amplificazione
                Progr_fiume.append(x1_fiume)

                # save the results in the respective lists
                QX.append(Q1)
                HX.append(h_1)
                BX.append(b_1)
                A_1=ka*math.pow(h_1,ma)
                if A_1>0:
                    V_1=Q1/A_1
                else:
                    V_1=0.0
                VX.append(V_1)

                TrattoX.append(tratto)
                PendX.append(pend)

                # look for the next point in tau1 in the hydrograph
                for kk in range(1,numtime):
                    if T[kk]>tau1:
                        k1=kk
                        break
                for kk in range(k1,numtime):
                    dttau=T[kk]-tau1
                    tau1=tau1+dttau
                    Q2=ValoreInX(tau1,numtime,T,Qin)
                    h_2=math.pow(Q2/kq,invmq)
                    # calculating celerity
                    c2=kcel*math.pow(h_2,mcel)
                    # time difference of the two points of the hydrograph at the end of the section
                    dtt=(1.0/c2-1.0/c1)*dx1+dttau

                    t1=t1+dtt
                    Q1=Q2*1.0
                    c1=c2*1.0
                    Tout.append(t1)
                    Qout.append(Q1)
                break

            else:
                # case in which the end of reach was not reached
                x=x+dx
                tau=tau+dtau
                t=t+dt
                # inserts a new point in the TQX curve
                TQX.append(t)
                x1=xini+x
                XQX.append(x1)
                # progressive distance upgrade along the river axis
                x1_fiume=xini_fiume+x*Amplificazione
                Progr_fiume.append(x1_fiume)

                Q1=ValoreInX(tau,numtime,T,Qin)
                h_1=math.pow(Q1/kq,invmq)
                # water width calculation
                b_1=kb*math.pow(h_1,mb)

                # save results
                QX.append(Q1)
                HX.append(h_1)
                BX.append(b_1)
                A_1=ka*math.pow(h_1,ma)
                if A_1>0:
                    V_1=Q1/A_1
                else:
                    V_1=0.0
                VX.append(V_1)

                TrattoX.append(tratto)
                PendX.append(pend)

                # -----------------------------------
                # check if the hydrograph has ended
                # -----------------------------------
                if k> (numtime-2):

                    # calculates the times to reach the end of reach
                    dx1=distanza-x
                    dt1=dx1/Vfronte
                    dtau1=denom/numerat*dt1
                    t1=t+dt1
                    tau1=tau+dtau1
                    # flow at instant tau1
                    Q1=ValoreInX(tau1,numtime,T,Qin)
                    h_1=math.pow(Q1/kq,invmq)
                    # water width calculation
                    b_1=kb*math.pow(h_1,mb)

                    # calculate the celerity
                    c1=kcel*math.pow(h_1,mcel)

                    # calculates the hydrograph leaving the reach
                    Tout=[]
                    Qout=[]
                    Tout.append(t1)
                    Qout.append(Q1)

                    # inserts a new point in the TQX curve
                    TQX.append(t1)
                    x1=xini+distanza
                    XQX.append(x1)

                    # progressive distance upgrade along the river axis
                    x1_fiume=xini_fiume+distanza*Amplificazione
                    Progr_fiume.append(x1_fiume)

                    # save the results in the respective lists
                    QX.append(Q1)
                    HX.append(h_1)
                    BX.append(b_1)
                    A_1=ka*math.pow(h_1,ma)
                    if A_1>0:
                        V_1=Q1/A_1
                    else:
                        V_1=0.0
                    VX.append(V_1)

                    TrattoX.append(tratto)
                    PendX.append(pend)

                    # extends the length of the hydrograph
                    Tout.append(t1+1000.0)
                    Qout.append(Q1)



        grafico1=0
        if grafico1>0:
            # praph
            fig1 = figure()

            ax1 = fig1.add_subplot(111)
            ax1.plot(T,Qin,'-')
            ax1.plot(Tout,Qout,'-')
            nomi=[]
            nomi.append('Inflow')
            nomi.append('Outflow')

            PosizoneLegenda = 'best'
            leg = ax1.legend((nomi),
                       PosizoneLegenda, shadow=True,prop = fontP)
            ax1.grid(True)
            ax1.set_ylabel('Q (mc/s)')
            ax1.set_xlabel('Time (hrs)')
            Titolo='CINE method applied to reach %d distance pixel n: %d  x=%.3f' % (i+1,tratto,x)

            suptitle(Titolo, fontsize=16,fontstyle='italic')

            show()

        #assigns the outgoing idrograph to input to the next reach
        Qin=numpy.array(Qout,dtype =numpy.float)
        T=numpy.array(Tout,dtype =numpy.float)
        numtime=len(T)

    if grafico2>0:

        # graph QX
        fig1 = figure()

        ax1 = fig1.add_subplot(111)
        ax1.plot(Progr_fiume,QX,'o-')
        nomi=[]
        nomi.append('QX')

        ax2 = ax1.twinx()
        ax2.plot(Progr_fiume,HX,'-r')

        nomi2=[]
        nomi2.append('HX')

        PosizoneLegenda = 'upper left'
        leg = ax1.legend((nomi), loc=PosizoneLegenda,
         shadow=True,prop = fontP)
        PosizoneLegenda2 = 'upper right'
        leg2 = ax2.legend((nomi2), loc=PosizoneLegenda2,
         shadow=True,prop = fontP)


        ax1.grid(True)
        ax1.set_ylabel('Qmax (mc/s)')
        ax2.set_ylabel('Hmax (m)')
        ax1.set_xlabel('x (m)')
        Titolo='Peak flow rate reduction with distance'

        suptitle(Titolo, fontsize=16,fontstyle='italic')

        show()

    # save the results
    # --------------------------

    NomeTabella='Q_H_max'

    # cancel any previous data
    # ---------------------------------
    sql='DELETE FROM %s WHERE DamID=%d' % (NomeTabella,DamID)
    cur.execute(sql)
    conn.commit()

    fout=open(FileOutput,'w')
    txt='PixDist;Progr_fiume;Progr_retta;pend;Qmax;Hmax;Bmax;Vmax;Time\n'
    fout.write(txt)
    n2=len(QX)
    for i in range(n2):
        # salva i risultati
        txt='%d;%.2f;%.2f;%s;%.2f;%.2f;%.2f;%.2f;%.2f\n' %(TrattoX[i],Progr_fiume[i],XQX[i],PendX[i],QX[i],HX[i],BX[i],VX[i],TQX[i])
        fout.write(txt)

        sql='INSERT INTO %s (DamID' % (NomeTabella)
        sql_value=') VALUES (%d' % DamID
        sql_value+=',%d,%.2f,%.2f,%s,%.2f,%.2f,%.2f,%.2f,%.2f' %(TrattoX[i],Progr_fiume[i],XQX[i],PendX[i],QX[i],HX[i],BX[i],VX[i],TQX[i])
        sql+=', PixDist'
        sql+=', Progr_fiume'
        sql+=', Progr_retta'
        sql+=', pend'
        sql+=', Qmax'
        sql+=', Hmax'
        sql+=', Bmax'
        sql+=', Vmax'
        sql+=', Time'
        sql+='%s' % sql_value
        sql+=');'
        cur.execute(sql)

    fout.close()

    conn.commit()

    # Close communication with the database
    cur.close()
    conn.close()


    return NotErr, errMsg

if __name__ == '__main__':

    mydb_path_user='..'+ os.sep+'db'+os.sep+'USER_GeoDB.sqlite'

    # San Giuliano
    DamID=449
    PathFiles='..'+ os.sep+ str(DamID)

    grafico2=1
    NotErr, errMsg= RunCinemat(mydb_path_user,DamID,PathFiles,grafico2)

    print(NotErr,errMsg)

