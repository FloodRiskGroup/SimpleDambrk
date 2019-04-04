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

Lo script, assegnato un ID_Diga calcola la propagazione dell'onda di piena
con il metodo di Cinematico + Onda a Fronte ripido

Legge i coeffienti ed esponenti delle caratteristiche geometriche ed idrauliche
dei vari tratti dal file:
    - MatriceAexp.csv

Calcola l'idrogramma di piena alla prima sezione di monte

A partire dal primo tratto calcola per ogni intervallo temporale la propagazione
della portata instantanea nel tratto ed il punto in cui raggiunge il fronte.
Calcola l'idrogramma in uscita dal tratto e lo utilizza come input per il tratto
successivo.

Fa il grafico dell'inviluppo delle portate ed altezze massime lungo la progressiva

Salva i risultati nel file:
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
    # importa i sistemi di riferimento
    from osgeo.osr import osr
except:
    import ogr
    # importa i sistemi di riferimento
    import osr

import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib.font_manager import FontProperties
import csv

def DamBreakHyrograph(VolumeInvaso,AltezzaDiga,LarghezzaSezioneRottura,NomeDiga):
    """
    Calcolo idrogramma di dambreak
    -------------------------------------
    VolumeInvaso            : in milioni mc
    AltezzaDiga             : in metri
    LarghezzaSezioneRottura : metri
    NomeDiga                : nome della diga
    """
    # passo a mc
    VolumeInvaso=VolumeInvaso*10**6

    # curva invaso : V=k*h^alfa
    # =========================
    alfa=1.5
    kinvaso=VolumeInvaso/(AltezzaDiga**alfa)

    Vcurr=VolumeInvaso
    Vfin=VolumeInvaso*0.001

    # calcolo portata per crollo istantaneo
	# ====================================
	# accelerazione di gravita'
    G=9.81
    # altezza alla sezione di rottura (Marchi Rubatta - MECCANICA DEI FLUIDI - UTET , Torino 1981 - pag 718-719)
    # corrisponde alla massima portata defluibile nell'ipotesi di eliminazione istantanea di una paratoria
    # che sorreggeva una colonna d'acqua ferma
    Yc=4.0/9.0*AltezzaDiga
    # portata massima da un lago inizialmente in quiete
    # e' inferiore a quella di una corrente gia' in moto permanente
    Qmax=8.0/27.0*LarghezzaSezioneRottura*AltezzaDiga*(G*AltezzaDiga)**0.5

    Vmax=Qmax/(Yc*LarghezzaSezioneRottura)
##    txt='Lbreccia=%.1f,  Yc=%.2f  , Qmax= %.1f , Vmax=%.2f' % (LarghezzaSezioneRottura,Yc,Qmax,Vmax)
##    print (txt)

    # Carico totale teorico nell'ipotesi di moto permanente
    Htot=Yc+Vmax**2/2./G
    # volume teorico alla quota del carico totale teorico
    Vol_Htot=kinvaso*Htot**alfa

    dt=60.0
    Idrogramma=[]
    Time=[]
    T=0.0
    Time.append(T)
    Idrogramma.append(Qmax)
    # assumo che la portata rimanga costante fino a quando il volume istantaneo
    # dell'invaso (Vcurr) diminuendo per effetto dello svuotamento non sia
    # inferiore a quello teorico che darebbe la stessa Qmax qualora il moto
    # fosse permanente
    while Vcurr>=Vol_Htot:
        dv=-Qmax*dt
        # aggiorna il volume d'invaso
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
    # altezza d'invaso che corrisponte al volume istantaneo Vcurr
    hcurr=math.pow((Vcurr/kinvaso),alfa_1)

    hmin=0.5
    # hcurr: livello d'invaso corrente e pari all'energia della corrente
    # in moto permanente
    # da questo istante calcolo la portata nell'ipotesi di altezza critica
    # nella sezione di rottura
    while hcurr>hmin:
        # profondita' critica per sezione rettangolare
        Yc=2.0/3.0*hcurr
        # velocita' critica per sezione rettangolare
        Uc=math.sqrt(G*Yc)
        # portata critica per sezione rettangolare
        Qc=Uc*Yc*LarghezzaSezioneRottura
        # volume rilasciato nell'unita' di tempo
        dv=-Qc*dt
        # aggiorna il volume d'invaso a fine dt
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
        # aggiorna il livello d'invaso a a fine dt
        hcurr=math.pow((Vcurr/kinvaso),alfa_1)

    # calcola volume idrogramma
    Vol=0.0
    nn=len(Idrogramma)
    for i in range(1,nn):
        # volume secondo il grafico (approssimato in eccesso!)
        Vol+=(Time[i]-Time[i-1])*(Idrogramma[i]+Idrogramma[i-1])/2.0
        # volume congruente con l'approssimazione di calcolo (grafico portata a gradini)
##        Vol+=(Time[i]-Time[i-1])*Idrogramma[i]
##    txt='Volume Idrogramma = %d' % (Vol)
##    print (txt)

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
    # data una funzione con ascisse datix ed ordinate datiy
    # si calcola, per interpolazione, il valore di y nell'ascissa x
    # se x e' minore di datix[0] alloca y vale datiy[0]
    # de x e' maggiore di datix[-1] allora y vale datiy[-1]
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

##def run_script(iface,ListaDati):
def RunCinemat(mydb_path_user,ID_Diga,PathFiles,grafico2):

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
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
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

    # leggo i dati della diga dal database
    sql='SELECT Volume_mlnm3,Altezza_m,Breccia_m,Nome FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
    cur.execute(sql)
    DatiDiga=cur.fetchone()

    if DatiDiga!=None:

        VolumeInvaso=float(DatiDiga[0])

        # AltezzaDiga
        AltezzaDiga = float(DatiDiga[1])

        # LarghezzaSezioneRottura
        LarghezzaSezioneRottura = float(DatiDiga[2])

        # Nome della diga
        NomeDiga=DatiDiga[3]

    else:
        errMsg = 'Nella tabella= %s non ci sono dati per la diga num =%s \nEffettuare prima il salvataggio dei dati !' % (NomeTabellaTest,ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    # lettura i dati per il calcolo dei parametri delle grandezze geometriche
    # ed idrauliche
    # lettura dal database
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
    sql+=' WHERE ID_Diga=%d' % (ID_Diga)
    sql+=' ORDER BY PixDist;'
    cur.execute(sql)
    MatriceDati=cur.fetchall()

    ListaTratti=[]
    Progressive_fiume=[]
    Distanze_fiume=[]
    # matrice dei coefficienti: ka;ma;kq;mq;kcel;mcel per ogni tratto
    MatriceCoeff=[]
    for row in MatriceDati:
        ListaTratti.append(int(row[0]))
        Progressive_fiume.append(float(row[1]))
        Distanze_fiume.append(float(row[2]))
        MatriceCoeff.append(row[3:])


    # Calcolo dell'idrogramma di piena
    # ................................

    Idrogramma, T = DamBreakHyrograph(VolumeInvaso,AltezzaDiga,LarghezzaSezioneRottura,NomeDiga)


    # il primo input corrisponde all'idrogramma alla sezione iniziale
    Qin=numpy.array(Idrogramma,dtype =numpy.float)
    maxQ=Qin.max()

    nn=len(MatriceCoeff)
    x=0.0

    numtime=len(Qin)

    # liste dei tempi, distanze e dell'inviluppo della portata massima lungo il fiume
    TQX=[]
    # distanza lungo la linea retta
    XQX=[]
    # portata
    QX=[]
    # altezza d'acqua
    HX=[]
    # larghezza pelo libero
    BX=[]
    # velocita
    VX=[]

    TrattoX=[]
    PendX=[]
    TQX.append(T[0])
    # vettore distanze in linea retta
    XQX.append(0.0)
    # vettore distanze lungo l'asse del fiume
    Progr_fiume=[]
    Progr_fiume.append(0.0)

    QX.append(Qin[0])
    TrattoX.append(0)
    PendX.append(MatriceCoeff[0][1])

    # coeff. ed esponente della formula monomia della portata Q=kq*h^mq
    kq=float(MatriceCoeff[0][4])
    mq=float(MatriceCoeff[0][5])
    invmq=1.0/mq
    h_0=math.pow(Qin[0]/kq,invmq)
    HX.append(h_0)

    # Larghezza pelo libero = ka*ma*h^(ma-1)
    # ----------------------------------------------------
    # coeff. ed esponente della formula monomia dell'area
    ka=float(MatriceCoeff[0][2])
    ma=float(MatriceCoeff[0][3])

    # calcolo larghezza pelo libero
    mb=ma-1.0
    b_0=ka*ma*math.pow(h_0,mb)
    BX.append(b_0)

    # velocita
    # ------------------
    A_0=ka*math.pow(h_0,ma)
    if A_0>0:
        V_0=Qin[0]/A_0
    else:
        V_0=0.0

    VX.append(V_0)


    # impostazioni per il grafico
    #============================
    fontP = FontProperties()
    fontP.set_size('small')

    Progr_fiume=[]
    Progr_fiume.append(0.0)


    # ciclo lungo i tratti di fiume
    # -------------------------------
    for i in range(nn):

        # ascissa ed istanti del punto iniziale del tratto precedente
        xini=XQX[-1]
        xini_fiume=Progr_fiume[-1]

        tini=TQX[-1]
        # ascissa a partire dall'inizio del tratto
        x=0.0
        xfronte=0.0

        # coefficienti kx ed esponenti mx di parametri geometrici ed idraulici dei tratti
        tratto=ListaTratti[i]

        # distanza in linea retta: utilizzata per il calcolo della propagazione
        distanza=float(MatriceCoeff[i][0])

        # lettura distanza sul fiume
        Progr_curr=Progressive_fiume[i]
        Parz_curr=Distanze_fiume[i]
        # rapporto fra la distanza lungo l'asse del fiume e la distanza in linea retta
        Amplificazione=Parz_curr/distanza

        # pendenza
        pend=float(MatriceCoeff[i][1])
        # coeff. ed esponente della formula monomia dell'area
        ka=float(MatriceCoeff[i][2])
        ma=float(MatriceCoeff[i][3])
        # coeff. ed esponente della formula monomia della portata Q=kq*h^mq
        kq=float(MatriceCoeff[i][4])
        mq=float(MatriceCoeff[i][5])
        # coeff. ed esponente della formula monomia della celerita'
        kcel=float(MatriceCoeff[i][6])
        mcel=float(MatriceCoeff[i][7])

        # parametri per il calcolo della larghezza del pelo libero B=dA/dh
        # ----------------------------------------------------------------
        # B=ka*ma*h^(m-1) = kb*h^mb
        kb=ka*ma
        mb=ma-1.0

        # calcolo del Dt
        # ..............
        t=T[0]
        tau=T[0]
        # calcolo l'altezza corrispondente
        invmq=1.0/mq

        for k in range(1,numtime):

            dtau=float(T[k]-T[k-1])

            Qtau=Qin[k-1]
            Qtau1=Qin[k]

            h_tau=math.pow(Qtau/kq,invmq)
            h_tau1=math.pow(Qtau1/kq,invmq)
            # calcolo dQ/dA a tau costante = celerita'
            Qprimo=kcel*math.pow(h_tau,mcel)
            # calcolo di dQ/dA^2 a tau costante = derivata della celerita'
            mcel2=mcel-1.0
            Qseconda=kcel*mcel*math.pow(h_tau,mcel2)
            # calcolo A(tau)
            Atau=ka*math.pow(h_tau,ma)
            Atau1=ka*math.pow(h_tau1,ma)
            # calcolo dA/dtau
            dA=(Atau1-Atau)/dtau
            # calcolo Vfronte che su alveo asciutto (cioe' Qo=0) vale Vfronte=Qtau/Atau
            Vfronte=Qtau/Atau
            # quindi calcolo dt
            numerat=Qseconda*dA*x/Qprimo-Qprimo
            denom=-Qprimo+Vfronte
            dt=numerat/denom*dtau
            # calcolo dx
            dx=dt*Vfronte

            # controllo la distanza percorsa
            if (x+dx)>distanza:
                # calcola i tempi per giungere a fine tratto
                dx1=distanza-x
                dt1=dx1/Vfronte
                dtau1=denom/numerat*dt1
                t1=t+dt1
                tau1=tau+dtau1
                # portata all'istante tau1
                Q1=ValoreInX(tau1,numtime,T,Qin)
                h_1=math.pow(Q1/kq,invmq)
                # calcolo larghezza pelo libero
                b_1=kb*math.pow(h_1,mb)

                # calcola la celerita'
                c1=kcel*math.pow(h_1,mcel)

                # calcola l'idrogramma in uscita dal tratto
                Tout=[]
                Qout=[]
                Tout.append(t1)
                Qout.append(Q1)

                # inserisce un nuovo punto nella curva TQX
                TQX.append(t1)
                x1=xini+distanza
                XQX.append(x1)

                # aggiornamento distanza progressiva lungo l'asse del fiume
                x1_fiume=xini_fiume+distanza*Amplificazione
                Progr_fiume.append(x1_fiume)

                # salva i risultati nelle rispettive liste
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

                # cerca nell'idrogramma il punto successivo a tau1
                for kk in range(1,numtime):
                    if T[kk]>tau1:
                        k1=kk
                        break
                for kk in range(k1,numtime):
                    dttau=T[kk]-tau1
                    tau1=tau1+dttau
                    Q2=ValoreInX(tau1,numtime,T,Qin)
                    h_2=math.pow(Q2/kq,invmq)
                    # calcola la celerita'
                    c2=kcel*math.pow(h_2,mcel)
                    # differenza temporale dei due punti dell'idrogramma alla fine del tratto
                    dtt=(1.0/c2-1.0/c1)*dx1+dttau

                    t1=t1+dtt
                    Q1=Q2*1.0
                    c1=c2*1.0
                    Tout.append(t1)
                    Qout.append(Q1)
                break

            else:
                # caso in cui non e' stato raggiunta la fine del tratto
                x=x+dx
                # caso in cui non e' stato raggiunta la fine del tratto
                tau=tau+dtau
                t=t+dt
                # inserisce un nuovo punto nella curva TQX
                TQX.append(t)
                x1=xini+x
                XQX.append(x1)
                # aggiornamento distanza progressiva lungo l'asse del fiume
                x1_fiume=xini_fiume+x*Amplificazione
                Progr_fiume.append(x1_fiume)

                Q1=ValoreInX(tau,numtime,T,Qin)
                h_1=math.pow(Q1/kq,invmq)
                # calcolo larghezza pelo libero
                b_1=kb*math.pow(h_1,mb)

                # salva i risultati
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
                # controllo se terminato l'idrogramma
                # -----------------------------------
                if k> (numtime-2):

                    # calcola i tempi per giungere a fine tratto
                    dx1=distanza-x
                    dt1=dx1/Vfronte
                    dtau1=denom/numerat*dt1
                    t1=t+dt1
                    tau1=tau+dtau1
                    # portata all'istante tau1
                    Q1=ValoreInX(tau1,numtime,T,Qin)
                    h_1=math.pow(Q1/kq,invmq)
                    # calcolo larghezza pelo libero
                    b_1=kb*math.pow(h_1,mb)

                    # calcola la celerita'
                    c1=kcel*math.pow(h_1,mcel)

                    # calcola l'idrogramma in uscita dal tratto
                    Tout=[]
                    Qout=[]
                    Tout.append(t1)
                    Qout.append(Q1)

                    # inserisce un nuovo punto nella curva TQX
                    TQX.append(t1)
                    x1=xini+distanza
                    XQX.append(x1)

                    # aggiornamento distanza progressiva lungo l'asse del fiume
                    x1_fiume=xini_fiume+distanza*Amplificazione
                    Progr_fiume.append(x1_fiume)

                    # salva i risultati nelle rispettive liste
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

                    # estende la lunghezza dell'idrogramma
                    Tout.append(t1+1000.0)
                    Qout.append(Q1)



        grafico1=0
        if grafico1>0:
            # grafico
            fig1 = figure()
            # crea il riquadro principale
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
            Titolo='Metodo CINE applicato al tratto %d distanza pixel n: %d  x=%.3f' % (i+1,tratto,x)

            suptitle(Titolo, fontsize=16,fontstyle='italic')

            show()

        # assegna l'idrogramma uscente in input al tratto successivo
        Qin=numpy.array(Qout,dtype =numpy.float)
        T=numpy.array(Tout,dtype =numpy.float)
        numtime=len(T)

    if grafico2>0:

        # grafico QX
        fig1 = figure()
        # crea il riquadro principale
        ax1 = fig1.add_subplot(111)
##        ax1.plot(XQX,QX,'o-')
        ax1.plot(Progr_fiume,QX,'o-')
        nomi=[]
        nomi.append('QX')

        ax2 = ax1.twinx()
##        ax2.plot(XQX,HX,'-r')
        ax2.plot(Progr_fiume,HX,'-r')

        nomi2=[]
        nomi2.append('HX')

##        PosizoneLegenda = 'best'
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
        Titolo='Laminazione della portata con la distanza'

        suptitle(Titolo, fontsize=16,fontstyle='italic')
        #title(Titolo,fontstyle='italic')

        show()

    # salva i risultati
    # --------------------------

    NomeTabella='Q_H_max'

    # cancella eventuali dati pregressi
    # ---------------------------------
    sql='DELETE FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
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

        sql='INSERT INTO %s (ID_Diga' % (NomeTabella)
        sql_value=') VALUES (%d' % ID_Diga
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
    ID_Diga=449
    PathFiles='..'+ os.sep+ str(ID_Diga)

    grafico2=1
    NotErr, errMsg= RunCinemat(mydb_path_user,ID_Diga,PathFiles,grafico2)

    print(NotErr,errMsg)

