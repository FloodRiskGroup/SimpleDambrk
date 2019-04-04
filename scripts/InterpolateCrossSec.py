# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SimpleDamBrk - CreaSezInterpolate
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

Lo script legge:
    - lo shape file dei punti di intersezione delle sezioni con l'asse del fiume
      ed i punti sull'asse ove richiesto la creazione delle sezioni interpolate
    - lo shape file delle sezioni da cui interpolare
    - lo shape file del contorno del dominio
    - il grid del DTM_clip
crea:

    - lo shape file delle sezioni totali (originali ed interpolate) clippate lungo
      il contorno del dominio

    - uno shape file poligonale composto da trapezi fra una sezione e la successiva

    - un grid con la classificazione delle celle secondo l'id dei poligoni che
      corrisponde al numero di celle lungo l'asse del fiume distanti dal punto
      iniziale

    - un grid Tratti.tif con le classi dei tratti di fiume a partire dalla diga
      il numero del tratto rappresenta in numero di celle di
      distanza contate lungo l'asse del fiume

    - uno shape file poligonale composto da trapezi fra una sezione e la successiva
      dividendo il trapezio in sinistra (lato=0) e destra (lato=1)

    - un grid DestraSinistra.tif con l'indicazione se la cella si trova a
      sinistra (=0) o a destra (=1) del fiume

    - un grid StreamDH con in valori delle altezze del terreno rispetto al fiume
"""
import os, sys
import numpy
try:
    from pyspatialite import dbapi2 as db
    import sqlite3
except:
    import sqlite3

try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    from gdalconst import *

try:
    from osgeo import ogr
    # importa i sistemi di riferimento
    from osgeo.osr import osr
except:
    import ogr
    # importa i sistemi di riferimento
    import osr
import numpy as np
import math
import matplotlib.pyplot as plt
import time
import scipy.interpolate as il

def isLeft(a,b,c):
    # controlla se il punto c(x,y) e' a sinistra del segmento a-b
    position=((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))

    return 1*numpy.sign(position)

##def SezioneInterpolata(pt_curr,SezMonteGeom,SezValleGeom):
##
##
##    # Meta' lunghezza sezione di monte
##    Length=SezMonteGeom.Length()/2.0
##
##    # calcola la geometria interpolata
##    # --------------------------------
##    p1 = numpy.array( [SezMonteGeom.GetX(0), SezMonteGeom.GetY(0)] )
##    p2 = numpy.array( [SezMonteGeom.GetX(1), SezMonteGeom.GetY(1)] )
##
##    p3 = numpy.array( [SezValleGeom.GetX(0), SezValleGeom.GetY(0)] )
##    p4 = numpy.array( [SezValleGeom.GetX(1), SezValleGeom.GetY(1)] )
##
##
##    # intersezione delle due sezioni di monte e valle
##    # -------------------------------------------------
##    # questo e' il punto a cui deve convergere la sezione interpolata
##    # e puo' trovarsi in destra o sinistra a seconda della curvatura
##    # dell'asse del fiume
##    Intersect_0=seg_intersect( p1,p2, p3,p4)
##
##    p1 = numpy.array( [pt_curr.GetX(0), pt_curr.GetY(0)] )
##    p2 = numpy.array( [Intersect_0[0], Intersect_0[1]] )
##
##    p3 = numpy.array( [SezMonteGeom.GetX(0), SezMonteGeom.GetY(0)] )
##    p4 = numpy.array( [SezValleGeom.GetX(0), SezValleGeom.GetY(0)] )
##
##    # primo punto della sezione interpolata  : quello in sinistra
##    Intersect_1=seg_intersect( p1,p2, p3,p4)
##
##    p3 = numpy.array( [SezMonteGeom.GetX(1), SezMonteGeom.GetY(1)] )
##    p4 = numpy.array( [SezValleGeom.GetX(1), SezValleGeom.GetY(1)] )
##
##    # secondo punto della sezione interpolata : quello in destra
##    Intersect_2=seg_intersect( p1,p2, p3,p4)
##
##    # estende la lunghezza della sezione
##    # -----------------------------------
##    GeomNew=SetCrossSecPointDiretion(p1,Length,Intersect_1,Intersect_2)
##
##
##    return GeomNew

def CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select):

    # controllo il verso della linea

    if inters_destra:
        # controllo distanza a destra
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(line_select.GetX(1),line_select.GetY(1))
        dist_dx=pt_interes_sez.Distance(point)
        point.Destroy()
        if dist_dx > distBufferClip:
            # invertire
            line= flip_line(line_select)
        else:
            line=line_select
    else:
        # controllo distanza a sinistra
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(line_select.GetX(0),line_select.GetY(0))
        dist_sx=pt_interes_sez.Distance(point)
        point.Destroy()
        if dist_sx > distBufferClip:
            # invertire
            line= flip_line(line_select)
        else:
            line=line_select

    return line


def graficoPoligono(poly1):
    """
    grafico debug un poligono
    """
    nn=poly1.GetGeometryCount()
    name=poly1.GetGeometryName()
    area_debug=poly1.GetArea()
    npt=poly1.GetPointCount()
    if name=="POLYGON":
        for k in range(nn):
            ring1=poly1.GetGeometryRef(k)
            n1= ring1.GetPointCount()
            if n1>0:
                x1=[]
                y1=[]
                for i in range(n1):
                    pt=ring1.GetPoint(i)
                    x1.append(pt[0])
                    y1.append(pt[1])
                plt.plot(x1,y1)
    else:
        for k in range(nn):
            ring1=poly1.GetGeometryRef(k)
            ngeo=ring1.GetGeometryCount()
            n1= ring1.GetPointCount()
            ring2=ring1.GetGeometryRef(0)
            n2=ring2.GetPointCount()
            if n2>0:
                x1=[]
                y1=[]
                for i in range(n2):
                    pt=ring2.GetPoint(i)
                    x1.append(pt[0])
                    y1.append(pt[1])
                plt.plot(x1,y1)
    plt.show()

def graficoPoligoni(poly1,poly2):
    """
    grafico debug due poligoni
    """
    nn=poly1.GetGeometryCount()

    for k in range(nn):
        ring1=poly1.GetGeometryRef(k)
        n1= ring1.GetPointCount()
        if n1>0:
            x1=[]
            y1=[]
            for i in range(n1):
                pt=ring1.GetPoint(i)
                x1.append(pt[0])
                y1.append(pt[1])
            plt.plot(x1,y1,'-k')

    ring2=poly2.GetGeometryRef(0)
    n2= ring2.GetPointCount()

    x2=[]
    y2=[]
    for i in range(n2):
        pt=ring2.GetPoint(i)
        x2.append(pt[0])
        y2.append(pt[1])

    plt.plot(x2,y2,'-r')
    plt.show()

def perp( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return
def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

def ProdottoVettSezioni(SezMonte,SezValle):
    """
    Calcola il prodotto vettoriale fra
    due vettori a e b nel piano x,y

    a (xA,yA)
    b (xB,yB)
    """
    xA=SezMonte.GetPoint(1)[0]-SezMonte.GetPoint(0)[0]
    yA=SezMonte.GetPoint(1)[1]-SezMonte.GetPoint(0)[1]

    xB=SezValle.GetPoint(1)[0]-SezValle.GetPoint(0)[0]
    yB=SezValle.GetPoint(1)[1]-SezValle.GetPoint(0)[1]

    axb=(xA*yB-yA*xB)

    return numpy.sign(axb)

def SetCrossSecPointDiretion(p0,Length,psx,pdx):
    """
    Crea la traccia di una sezione avente una direzione definita da due punti
    e centrata su un punto
    p0 [x,y]    : punto centrale
    Length      : semilunghezza
    psx [x,y]   : punto sinistro della direzione
    pdx [x,y]   : punto destro della direzione
    """
    # trovo l'angolo
    deltaX=pdx[0]-psx[0]
    deltaY=pdx[1]-psx[1]

    if deltaX!=0:
        arco=numpy.arctan(deltaY/deltaX)
        dx=numpy.sign(deltaX)*abs(Length*numpy.cos(arco))
        dy=numpy.sign(deltaY)*abs(Length*numpy.sin(arco))

        x1=p0[0]-dx
        y1=p0[1]-dy
        x2=p0[0]+dx
        y2=p0[1]+dy
    else:
        dy=numpy.sign(deltaY)*Length
        x1=p0[0]
        y1=p0[1]-dy
        x2=p0[0]
        y2=p0[1]+dy

    GeomNew=ogr.Geometry(ogr.wkbLineString)

    GeomNew.AddPoint(x1,y1)
    GeomNew.AddPoint(x2,y2)

    return GeomNew

##def SezioneInterpolata(num,pt_curr,SezMonteGeom,SezValleGeom):
##def SezioneInterpolata(num,pt_curr,SezMonteGeom,SezValleGeom):
def SezioneInterpolata(pt_curr,SezMonteGeom,SezValleGeom):
    """
    Trova per interpolazione una sezione intermedia fra
    - SezMonteGeom: geom del segmento della sezione a monte
    - SezValleGeom: geom del segmento della sezione a monte
    - pt_curr     : punto intermendio fra la due sezioni e da cui passa
                    la sezione interpolata
                    il secondo punto della sezione interpolata e'
                    quello in cui si incontrano le rette passanti
                    per le due sezioni di monte e di valle
    """
##    Toll=0.1
##    num_controllo=4

    # Meta' lunghezza sezione di monte
    Length=SezMonteGeom.Length()/2.0*1.5

    # calcola la geometria interpolata
    # --------------------------------
    p1 = np.array( [SezMonteGeom.GetX(0), SezMonteGeom.GetY(0)] )
    p2 = np.array( [SezMonteGeom.GetX(1), SezMonteGeom.GetY(1)] )

##    if num==num_controllo:
##        plt.plot([p1[0],p2[0]],[p1[1],p2[1]],'-k')

    p3 = np.array( [SezValleGeom.GetX(0), SezValleGeom.GetY(0)] )
    p4 = np.array( [SezValleGeom.GetX(1), SezValleGeom.GetY(1)] )

##    if num==num_controllo:
##        pass
##        plt.plot([p3[0],p4[0]],[p3[1],p4[1]],'-g')


    # intersezione delle due sezioni di monte e valle
    # -------------------------------------------------
    # questo e' il punto a cui deve convergere la sezione interpolata
    # e puo' trovarsi in destra o sinistra a seconda della curvatura
    # dell'asse del fiume
    Intersect_0=seg_intersect( p1,p2, p3,p4)

    # controllo se a destra o sinistra del fiume
    orientamento=ProdottoVettSezioni(SezMonteGeom,SezValleGeom)

    # controllo sulla linea di unione delle sezioni direttici alla parte interna della curva
    LineInterna=ogr.Geometry(ogr.wkbLineString)

    if  orientamento>0:
        # si presume una intersezione in sponda sinistra
        inters_sx=bool('True')
        inters_destra=bool()
        LineInterna.AddPoint(SezMonteGeom.GetX(0), SezMonteGeom.GetY(0))
        LineInterna.AddPoint(SezValleGeom.GetX(0), SezValleGeom.GetY(0))

    else:
        inters_sx=bool()
        inters_destra=bool('True')
        LineInterna.AddPoint(SezMonteGeom.GetX(1), SezMonteGeom.GetY(1))
        LineInterna.AddPoint(SezValleGeom.GetX(1), SezValleGeom.GetY(1))


     # estende la lunghezza della sezione
    # -----------------------------------
    if  inters_sx:
        GeomNew=SetCrossSecPointDiretion(pt_curr.GetPoint(0),Length,Intersect_0,pt_curr.GetPoint(0))
    else:
        GeomNew=SetCrossSecPointDiretion(pt_curr.GetPoint(0),Length,pt_curr.GetPoint(0),Intersect_0)



    # controllo che la sezione non raggiunga il punto di intersezione delle principali
    # --------------------------------------------------------------------------------
    bufferClip=50.0
    distBufferClip=bufferClip*1.1

    # crea un buffer intorno al punto
    pt_interes_sez=ogr.Geometry(ogr.wkbPoint)
    pt_interes_sez.AddPoint(Intersect_0[0],Intersect_0[1])
    pt_poly = pt_interes_sez.Buffer(bufferClip)


    if GeomNew.Crosses(pt_poly):

        pt_curr_poly= pt_curr.Buffer(bufferClip)

        # la linea viene accorciata
        lines=GeomNew.Difference(pt_poly)

        # conta il numero di geometrie
        nn=lines.GetGeometryCount()
        if nn>1:
            for i in range(nn):
                g = lines.GetGeometryRef(i)
                # sceglie la parte di linea che interseca l'asse del fiume
                if g.Crosses(pt_curr_poly):
                    line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                    # controllare il verso
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                    break
        else:
            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

        wkt_new=line.ExportToWkt()
        GeomNew=ogr.CreateGeometryFromWkt(wkt_new)

    # controlla l'intersezione con la linea uniona delle sezioni principali
    # ----------------------------------------------------------------------
    line_poly=LineInterna.Buffer(1.0)

    if GeomNew.Crosses(LineInterna):

        pt_curr_poly= pt_curr.Buffer(bufferClip)

        # trova il punto di intersezione
        pt_interes_sez=GeomNew.Intersection(LineInterna)


        # la linea viene accorciata
        lines=GeomNew.Difference(line_poly)

        # conta il numero di geometrie
        nn=lines.GetGeometryCount()
        if nn>1:
            for i in range(nn):
                g = lines.GetGeometryRef(i)
                # sceglie la parte di linea che interseca l'asse del fiume
                if g.Crosses(pt_curr_poly):
                    line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                    # controllare il verso
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                    break
        else:
            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

        wkt_new=line.ExportToWkt()
        GeomNew=ogr.CreateGeometryFromWkt(wkt_new)


##    if num==num_controllo:
##        plt.plot([GeomNew.GetPoint(0)[0],GeomNew.GetPoint(1)[0]],[GeomNew.GetPoint(0)[1],GeomNew.GetPoint(1)[1]],'-r')
##        plt.show()

##    point_intersect = ogr.Geometry(ogr.wkbPoint)
##    point_intersect.AddPoint(Intersect_0[0],Intersect_0[1])
##
##    dist_monte=SezMonteGeom.Distance(point_intersect)
##    dist_valle=SezValleGeom.Distance(point_intersect)
##
##    chk=bool()
##    if  dist_monte<Toll or dist_valle<Toll:
##        chk=bool('True')
##        if dist_monte<dist_valle:
##            inters_monte=bool('True')
##        else:
##            inters_monte=bool()
##        # controllo se a destra o sinistra del fiume
##        dist_sx=math.sqrt((SezMonteGeom.GetX(0)-Intersect_0[0])**2+(SezMonteGeom.GetY(0)-Intersect_0[1])**2)
##        dist_dx=math.sqrt((SezMonteGeom.GetX(1)-Intersect_0[0])**2+(SezMonteGeom.GetY(1)-Intersect_0[1])**2)
##        if  dist_sx< dist_dx:
##            inters_sx=bool('True')
##        else:
##            inters_sx=bool()
##
##
##    if num==num_controllo:
##        plt.scatter(pt_curr.GetX(0), pt_curr.GetY(0),marker='o', c="b")
##        plt.scatter(Intersect_0[0], Intersect_0[1],marker='o', c="r")
##
##
##    p1 = np.array( [pt_curr.GetX(0), pt_curr.GetY(0)] )
##    p2 = np.array( [Intersect_0[0], Intersect_0[1]] )
##
##    if  chk and inters_sx:
##        if inters_monte:
##            x_medio=(SezMonteGeom.GetX(0)+SezMonteGeom.GetX(1))/2.0
##            y_medio=(SezMonteGeom.GetY(0)+SezMonteGeom.GetY(1))/2.0
##            xx=(x_medio+Intersect_0[0])/2.0
##            yy=(y_medio+Intersect_0[1])/2.0
##            p3 = np.array( [xx, yy])
##            p4 = np.array( [SezValleGeom.GetX(0), SezValleGeom.GetY(0)] )
##        else:
##            x_medio=(SezValleGeom.GetX(0)+SezValleGeom.GetX(1))/2.0
##            y_medio=(SezValleGeom.GetY(0)+SezValleGeom.GetY(1))/2.0
##            xx=(x_medio+Intersect_0[0])/2.0
##            yy=(y_medio+Intersect_0[1])/2.0
##            p3 = np.array( [SezMonteGeom.GetX(0), SezMonteGeom.GetY(0)])
##            p4 = np.array( [xx, yy])
##    else:
##        p3 = np.array( [SezMonteGeom.GetX(0), SezMonteGeom.GetY(0)])
##        p4 = np.array( [SezValleGeom.GetX(0), SezValleGeom.GetY(0)] )
##
##    # primo punto della sezione interpolata  : quello in sinistra
##    Intersect_1=seg_intersect( p1,p2, p3,p4)
##
##    if num==num_controllo:
##        plt.scatter(Intersect_1[0], Intersect_1[1],marker='o', c="g")
##
##    # controllo per il punto a destra
##    if  chk and not inters_sx:
##        if inters_monte:
##            x_medio=(SezMonteGeom.GetX(0)+SezMonteGeom.GetX(1))/2.0
##            y_medio=(SezMonteGeom.GetY(0)+SezMonteGeom.GetY(1))/2.0
##            xx=(x_medio+Intersect_0[0])/2.0
##            yy=(y_medio+Intersect_0[1])/2.0
##            p3 = np.array( [xx, yy])
##            p4 = np.array( [SezValleGeom.GetX(1), SezValleGeom.GetY(1)])
##        else:
##            x_medio=(SezValleGeom.GetX(0)+SezValleGeom.GetX(1))/2.0
##            y_medio=(SezValleGeom.GetY(0)+SezValleGeom.GetY(1))/2.0
##            xx=(x_medio+Intersect_0[0])/2.0
##            yy=(y_medio+Intersect_0[1])/2.0
##            p3 = np.array( [SezMonteGeom.GetX(1), SezMonteGeom.GetY(1)])
##            p4 = np.array( [xx, yy])
##    else:
##        p3 = np.array( [SezMonteGeom.GetX(1), SezMonteGeom.GetY(1)] )
##        p4 = np.array( [SezValleGeom.GetX(1), SezValleGeom.GetY(1)] )
##
##    # secondo punto della sezione interpolata : quello in destra
##    Intersect_2=seg_intersect( p1,p2, p3,p4)
##
##    if num==num_controllo:
##        plt.scatter(Intersect_2[0], Intersect_2[1], marker='^', c="g")
##
##    # estende la lunghezza della sezione
##    # -----------------------------------
##    # trovo l'angolo
##    deltaX=Intersect_2[0]-Intersect_1[0]
##    deltaY=Intersect_2[1]-Intersect_1[1]
##
##    if deltaX!=0:
##        arc=np.arctan(deltaY/deltaX)
##        dx=np.sign(deltaX)*abs(Length*np.cos(arc))
##        dy=np.sign(deltaY)*abs(Length*np.sin(arc))
##
##        x1=p1[0]-dx
##        y1=p1[1]-dy
##        x2=p1[0]+dx
##        y2=p1[1]+dy
##    else:
##        dy=np.sign(deltaY)*Length
##        x1=p1[0]
##        y1=p1[1]-dy
##        x2=p1[0]
##        y2=p1[1]+dy
##
##
##    if num==num_controllo:
##
##        plt.plot([Intersect_1[0],Intersect_2[0]],[Intersect_1[1],Intersect_2[1]],'-r')
##        plt.show()
##
##    GeomNew=ogr.Geometry(ogr.wkbLineString)
##
##    GeomNew.AddPoint(x1,y1)
##    GeomNew.AddPoint(x2,y2)

    return GeomNew

def ChkDimSez(line,pt):
    # controlla le distanze dei punti della sezione dal punto centrale
    pass

def PuntoMedio(pt1,pt2):

    # calcola il punto medio fra due
    pt_mean=[(pt1[0]+pt2[0])/2.0,(pt1[1]+pt2[1])/2.0]

    return pt_mean

def PuntoIntermedio(pt1,pt2,percent):

    # calcola il punto intermedio ad una cerca percentuale fra due
    pt_intermedio=[pt1[0]+percent*(pt2[0]-pt1[0]),pt1[1]+percent*(pt2[1]-pt1[1])]

    return pt_intermedio

def flip_line(line):
    # inverte i punti della linea
    GeomNew=ogr.Geometry(ogr.wkbLineString)
    for i in range(line.GetPointCount()-1,-1,-1):
        # GetPoint returns a tuple not a Geometry
        pt = line.GetPoint(i)
        GeomNew.AddPoint(pt[0], pt[1])

    return GeomNew


def SetStreamLinePoints(MultiStreamLine,End_cross_sec_line):
    """
    Crea la lista dei punti del fiume dalla prima all'ultima sezione
    Input:
        -  StreamLine           : linea del fiume
        -  End_cross_sec_line   : linea della sezione finale
    """
    # crea la lista dei punti della streamLine
    ListaPuntiStreamLine=[]
    if MultiStreamLine.GetGeometryName() == 'LINESTRING':
        StreamLine=MultiStreamLine
    else:

        # controlla il numero di geometrie
        # prende il poligono piu' grande nei Multipoligons
        numlines=MultiStreamLine.GetGeometryCount()

        if numlines>1:
            Amax=0.0
            ii_max=0
            for ii in range(numlines):
                g=MultiStreamLine.GetGeometryRef(ii)
                Aii=g.Length()
                if  Aii>Amax:
                    Amax=Aii
                    ii_max=ii
            StreamLine=MultiStreamLine.GetGeometryRef(ii_max)
            Length=StreamLine.Length()

        else:
            # prende la prima
            StreamLine=MultiStreamLine.GetGeometryRef(0)
            Length=StreamLine.Length()

    # inserisco il primo punto
    ListaPuntiStreamLine.append(StreamLine.GetPoint(0))

    for ii in range(0, StreamLine.GetPointCount()-1):
        segm_curr=ogr.Geometry(ogr.wkbLineString)
        segm_curr.AddPoint(StreamLine.GetPoint(ii)[0],StreamLine.GetPoint(ii)[1])
        segm_curr.AddPoint(StreamLine.GetPoint(ii+1)[0],StreamLine.GetPoint(ii+1)[1])

        if segm_curr.Intersect(End_cross_sec_line):
            pt_intersect=segm_curr.Intersection(End_cross_sec_line)
            pt_end=pt_intersect.GetPoint(0)
            ListaPuntiStreamLine.append(pt_end)
            break

        else:
            ListaPuntiStreamLine.append(StreamLine.GetPoint(ii+1))

    return ListaPuntiStreamLine

def SetIntermediatePoints(mydb_path_user,ID_Diga,PathFiles,fileDEM,DistanzaSezInterp=2000.0,DeltaSezPrincipale=5):
    """
    Crea la sequenza dei punti intermedi fra le sezioni principali
    """
    NotErr=bool('True')
    errMsg='OK'

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima l'editing delle sezioni principali !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    if not os.path.exists(fileDEM):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il clip del DEM a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg


    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
   # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    # controlla esistenza linea tratto a valle
    NomeTabellaLinee='LineaValleDiga'

    # codice del sistema di riferimento della tabella
    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
    cur.execute(sql)

    record=cur.fetchone()
    if record!=None:
        DominioEPSG=record[0]
    else:
        DominioEPSG=32632

    dest_srs = ogr.osr.SpatialReference()
    dest_srs.ImportFromEPSG(DominioEPSG)


    sql='SELECT ST_AsText(geom) FROM %s WHERE ID_Diga=%d' % (NomeTabellaLinee,ID_Diga)
    cur.execute(sql)
    ChkDiga=cur.fetchone()

    if ChkDiga==None:
        errMsg = "Nella tabella= %s non ci sono dati per la diga num =%s \nEffettuare prima il calcolo della linea a valle !" % (NomeTabellaLinee,ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt=ChkDiga[0]
        linea_geom=ogr.CreateGeometryFromWkt(wkt)
        TotalLength=linea_geom.Length()


    # tabella delle sezioni principali
    NomeTabella='MainCrossSec'

    sql='SELECT PKUID,ST_AsText(geom) FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
    sql+=' ORDER BY id'
    sql+=';'
    cur.execute(sql)
    ListaCrossSec=cur.fetchall()

    n=len(ListaCrossSec)

    if n <= 0:

        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni principali !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    else:

        # codice del sistema di riferimento della tabella
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        cur.execute(sql)
        record=cur.fetchone()
        if record!=None:
            OriginEPSG=record[0]
        else:
            OriginEPSG=32632


    # --------------
    # leggo il DTM
    # --------------

    gdal.AllRegister()

    indataset = gdal.Open(fileDEM, GA_ReadOnly )
    if indataset is None:
        errMsg = "Could not open " + fileDEM
        NotErr= bool()
        return NotErr, errMsg

    geotransform = indataset.GetGeoTransform()

    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    cols=indataset.RasterXSize
    rows=indataset.RasterYSize
    bands=indataset.RasterCount
    iBand = 1
    inband = indataset.GetRasterBand(iBand)
    inNoData= inband.GetNoDataValue()

    Elevation = inband.ReadAsArray(0, 0, cols, rows).astype(numpy.float32)

    inband = None
    indataset=None

    # definisco una distanza minima fra sezioni secondarie
    DistanzaSezInterp_min=3*pixelWidth

    # Crea la lista dei punti della linea del fiume
    # ---------------------------------------------
    ListaPuntiFiume=[]
    for i in range(0, linea_geom.GetPointCount()):
        ListaPuntiFiume.append(linea_geom.GetPoint(i))

    num_pt_fiume=len(ListaPuntiFiume)

    # dizionario PKUID_sec_id: memorizza i valori id da assegnare alle sezioni principali
    PKUID_sec_id={}
    listaIdSecPrincipali=[]

    # dizionario id_Points : dizionario dei punti
    id_Points={}
    # dizionario delle distanze progressive
    progr_Points={}

    # Scorre le varie sezioni
    # -----------------------
    point_stream_start=0
    num_point_curr=0
    progr_monte=0.0
    pt_start=ListaPuntiFiume[0]
    id_Points[num_point_curr]=pt_start
    progr_Points[num_point_curr]=progr_monte

    # cerco la prima sezione principale
    cross_sec_wkt=ListaCrossSec[0][1]

    # assegno id=0 al punto della prima sezione
    PKUID_sec_id[ListaCrossSec[0][0]]=num_point_curr
    listaIdSecPrincipali.append(num_point_curr)

    cross_sec_line=ogr.CreateGeometryFromWkt(cross_sec_wkt)

    # memorizza la linea della sezione interpolata corrispondente
    UpstreaSez_wkt=cross_sec_line.ExportToWkt()
    # memorizza i punti della sezione
    Point_a=cross_sec_line.GetPoint(0)
    Point_b=cross_sec_line.GetPoint(1)

    line_tratto=ogr.Geometry(ogr.wkbLineString)
    line_tratto.AddPoint(pt_start[0],pt_start[1])


    ListaSecOk=[]
    NumSecPrincipali=len(ListaCrossSec)

    for isec in range(1,NumSecPrincipali):

        next_cross_sec_wkt=ListaCrossSec[isec][1]
        next_cross_sec_line=ogr.CreateGeometryFromWkt(next_cross_sec_wkt)


        for i_fiume in range(point_stream_start,num_pt_fiume-1):
            segm_curr=ogr.Geometry(ogr.wkbLineString)
            segm_curr.AddPoint(ListaPuntiFiume[i_fiume][0],ListaPuntiFiume[i_fiume][1])
            segm_curr.AddPoint(ListaPuntiFiume[i_fiume+1][0],ListaPuntiFiume[i_fiume+1][1])

            if segm_curr.Intersect(next_cross_sec_line):

                ListaSecOk.append(isec)

                pt_intersect=segm_curr.Intersection(next_cross_sec_line)
                pt_start=pt_intersect.GetPoint(0)
                line_tratto.AddPoint(pt_start[0],pt_start[1])
                tratto_length=line_tratto.Length()

                DistanzaSezInterpCurr=tratto_length/float(DeltaSezPrincipale)

                if DistanzaSezInterp<DistanzaSezInterp_min:
                    DeltaSezPrincipaleCur=int(tratto_length/DistanzaSezInterp_min)
                    if DeltaSezPrincipaleCur<1:
                        DeltaSezPrincipaleCur=1
                else:
                    DeltaSezPrincipaleCur=DeltaSezPrincipale*1


                NumSezSecondarieTratto=DeltaSezPrincipaleCur+1

                # sequenza della distanze progressive secondarie
                xvals=numpy.linspace(0, tratto_length, NumSezSecondarieTratto)

                # generazione dei punti intermedi
                listax=[]
                listay=[]
                for kk in range(0, line_tratto.GetPointCount()):
                    x_cur=line_tratto.GetPoint(kk)[0]
                    y_cur=line_tratto.GetPoint(kk)[1]
                    listax.append(x_cur)
                    listay.append(y_cur)

                # get the cumulative distance along the line
                x0=numpy.array(listax)
                y0=numpy.array(listay)
                dist = numpy.sqrt((x0[:-1] - x0[1:])**2 + (y0[:-1] - y0[1:])**2)
                dist_along = numpy.concatenate(([0], dist.cumsum()))
                # creo il vettore delle progressive
                progr_along=xvals + progr_monte

                # trovo i punti alle distanze predefinite mediante interpolazione
                ArrayXSec=numpy.interp(xvals, dist_along, x0)
                ArrayYSec=numpy.interp(xvals, dist_along, y0)


                # creazione dei punti intermedi
                # -----------------------------
                npt=len(ArrayXSec)

                pt_valle=[ArrayXSec[-1],ArrayYSec[-1]]

                Tooll_Distanza=tratto_length/float(npt)/4.0

                for i in range(1,npt-1):

                    num_point_curr+=1

                    pt_curr=ogr.Geometry(ogr.wkbPoint)
                    pt_curr.AddPoint(ArrayXSec[i],ArrayYSec[i])

                    # controllo la posizione punto i-esimo rispetto al segmento di monte
                    c=pt_curr.GetPoint(0)
                    PosizPtCur=isLeft(Point_a,Point_b,c)

                    # controllo la posizione punto a valle rispetto al segmento di valle
                    PosizPtValle=isLeft(Point_a,Point_b,pt_valle)


                    # salva i dati del segmento di monte per il controllo del punto successiva
                    NewGeom = SezioneInterpolata(pt_curr,cross_sec_line,next_cross_sec_line)
                    Point_a=NewGeom.GetPoint(0)
                    Point_b=NewGeom.GetPoint(1)

                    # controlla la posizione reciprova attuale dei due punti
                    Congruenza=PosizPtCur*PosizPtValle

                    if Congruenza>0:

                        # controllo anche che vi sia una distanza minima dalla sezione di monte
                        # e dalla sezione direttrice di valle
                        UpstreaSez=ogr.CreateGeometryFromWkt(UpstreaSez_wkt)
                        Dist1=pt_curr.Distance(UpstreaSez)
                        Dist2=pt_curr.Distance(next_cross_sec_line)

                        if Dist1>Tooll_Distanza and Dist2>Tooll_Distanza:

                            # salvo le coordinate del punto
                            id_Points[num_point_curr]=pt_curr.GetPoint(0)
                            progr_Points[num_point_curr]=progr_along[i]

                            # salva la geometria della sezione a monte del punto i+1 -esimo
                            UpstreaSez_wkt=NewGeom.ExportToWkt()

                        else:
                            txt='--  Scartato punto n: %d delle sezioni intermedie a valle sez. direttrice : %d-esima causa distanza inferiore alla minima' % (i,isec-1)
                            errMsg+='\n%s'% txt

                            pt_curr.Destroy()

                    else:
                        txt='--  Scartato punto n: %d  delle sezioni intermedie a valle sez. direttice :%d-esima causa meandro che verso monte' % (i,isec-1)
                        errMsg+='\n%s'% txt
                        pt_curr.Destroy()

                # -----------------------------------------------------
                # Creazione del punto della sezione principale di valle
                # -----------------------------------------------------

                # incremento il numero di id del punto
                num_point_curr+=1

                # creazione del punto di valle
                PKUID_sec_id[ListaCrossSec[isec][0]]=num_point_curr
                listaIdSecPrincipali.append(num_point_curr)

                # salvo le coordinate del punto
                id_Points[num_point_curr]=pt_valle
                progr_Points[num_point_curr]=progr_along[-1]

                # aggiorno la progressiva di monte
                progr_monte+=tratto_length

                # salva la geometria della sezione a monte del punto i+1 -esimo
                UpstreaSez_wkt=ListaCrossSec[isec][1]

                cross_sec_line=ogr.CreateGeometryFromWkt(UpstreaSez_wkt)

                # memorizza i punti della sezione
                Point_a=cross_sec_line.GetPoint(0)
                Point_b=cross_sec_line.GetPoint(1)

                # inizializza il nuovo tratto
                line_tratto.Destroy()
                line_tratto=ogr.Geometry(ogr.wkbLineString)
                line_tratto.AddPoint(pt_valle[0],pt_valle[1])

                # aggiorna e esce dal ciclo
                point_stream_start=i_fiume

                break

                # azzero la memoria
                line_tratto.Destroy()

            else:
                line_tratto.AddPoint(ListaPuntiFiume[i_fiume+1][0],ListaPuntiFiume[i_fiume+1][1])

    # controlla il caso dell'ultima sezione
    # -------------------------------------
    NumTratti= NumSecPrincipali-1

    if NumTratti not in ListaSecOk:

        # caso in cui l'ultimo tratto (line_tratto) non e' stato valutato
        tratto_length=line_tratto.Length()

        DistanzaSezInterpCurr=tratto_length/float(DeltaSezPrincipale)

        if DistanzaSezInterp<DistanzaSezInterp_min:
            DeltaSezPrincipaleCur=int(tratto_length/DistanzaSezInterp_min)
            if DeltaSezPrincipaleCur<1:
                DeltaSezPrincipaleCur=1
        else:
            DeltaSezPrincipaleCur=DeltaSezPrincipale*1


        NumSezSecondarieTratto=DeltaSezPrincipaleCur+1

        # sequenza della distanze progressive secondarie
        xvals=numpy.linspace(0, tratto_length, NumSezSecondarieTratto)

        # generazione dei punti intermedi
        listax=[]
        listay=[]
        for kk in range(0, line_tratto.GetPointCount()):
            x_cur=line_tratto.GetPoint(kk)[0]
            y_cur=line_tratto.GetPoint(kk)[1]
            listax.append(x_cur)
            listay.append(y_cur)

        # get the cumulative distance along the line
        x0=numpy.array(listax)
        y0=numpy.array(listay)
        dist = numpy.sqrt((x0[:-1] - x0[1:])**2 + (y0[:-1] - y0[1:])**2)
        dist_along = numpy.concatenate(([0], dist.cumsum()))
        # creo il vettore delle progressive
        progr_along=xvals + progr_monte

        # trovo i punti alle distanze predefinite mediante interpolazione
        ArrayXSec=numpy.interp(xvals, dist_along, x0)
        ArrayYSec=numpy.interp(xvals, dist_along, y0)


        # creazione dei punti intermedi
        # -----------------------------
        npt=len(ArrayXSec)

        pt_valle=[ArrayXSec[-1],ArrayYSec[-1]]

        Tooll_Distanza=tratto_length/float(npt)/4.0

        for i in range(1,npt-1):

            num_point_curr+=1

            pt_curr=ogr.Geometry(ogr.wkbPoint)
            pt_curr.AddPoint(ArrayXSec[i],ArrayYSec[i])

            # controllo la posizione punto i-esimo rispetto al segmento di monte
            c=pt_curr.GetPoint(0)
            PosizPtCur=isLeft(Point_a,Point_b,c)

            # controllo la posizione punto a valle rispetto al segmento di valle
            PosizPtValle=isLeft(Point_a,Point_b,pt_valle)


            # salva i dati del segmento di monte per il controllo del punto successiva
            NewGeom = SezioneInterpolata(pt_curr,cross_sec_line,next_cross_sec_line)
            Point_a=NewGeom.GetPoint(0)
            Point_b=NewGeom.GetPoint(1)

            # controlla la posizione reciprova attuale dei due punti
            Congruenza=PosizPtCur*PosizPtValle

            if Congruenza>0:

                # controllo anche che vi sia una distanza minima dalla sezione di monte
                # e dalla sezione direttrice di valle
                UpstreaSez=ogr.CreateGeometryFromWkt(UpstreaSez_wkt)
                Dist1=pt_curr.Distance(UpstreaSez)
                Dist2=pt_curr.Distance(next_cross_sec_line)

                if Dist1>Tooll_Distanza and Dist2>Tooll_Distanza:

                    # salvo le coordinate del punto
                    id_Points[num_point_curr]=pt_curr.GetPoint(0)
                    progr_Points[num_point_curr]=progr_along[i]

                    # salva la geometria della sezione a monte del punto i+1 -esimo
                    UpstreaSez_wkt=NewGeom.ExportToWkt()

                else:
                    txt='--  Scartato punto n: %d delle sezioni intermedie a valle sez. direttrice : %d-esima causa distanza inferiore alla minima' % (i,isec-1)
                    errMsg+='\n%s'% txt

                    pt_curr.Destroy()

            else:
                txt='--  Scartato punto n: %d  delle sezioni intermedie a valle sez. direttice :%d-esima causa meandro che verso monte' % (i,isec-1)
                errMsg+='\n%s'% txt
                pt_curr.Destroy()

        # -----------------------------------------------------
        # Creazione del punto della sezione principale di valle
        # -----------------------------------------------------

        # incremento il numero di id del punto
        num_point_curr+=1

        # creazione del punto di valle
        PKUID_sec_id[ListaCrossSec[isec][0]]=num_point_curr
        listaIdSecPrincipali.append(num_point_curr)

        # salvo le coordinate del punto
        id_Points[num_point_curr]=pt_valle
        progr_Points[num_point_curr]=progr_along[-1]

        # azzero la memoria
        line_tratto.Destroy()

    # trova le quote di fondo dei punti
    # ---------------------------------

    # trovo le distanze progressive delle sezioni secondarie
    # -------------------------------------------------------
    XSec=[]
    YSec=[]
    ListaId=[]
    for id_pt in id_Points:
        ListaId.append(id_pt)
        XSec.append(id_Points[id_pt][0])
        YSec.append(id_Points[id_pt][1])

    ArrayXSec=numpy.array(XSec)
    ArrayYSec=numpy.array(YSec)

    # distanze progressive delle sezioni secondarie
    # -------------------------------------------------------
    dist = numpy.sqrt((ArrayXSec[:-1] - ArrayXSec[1:])**2 + (ArrayYSec[:-1] - ArrayYSec[1:])**2)
    dist_along_2 = numpy.concatenate(([0], dist.cumsum()))
    nn_along_2=len(dist_along_2)


    ArrayIdSecPrinc=numpy.array(listaIdSecPrincipali)


    QuoteAlveo=[]
    listax=[]
    listay=[]

    # effettuo la ricerca della quota minima
    # ...........................................................
    # controllando anche un certo numero di celle intorno al punto
    deltapixel=1
    steps=deltapixel*2+1
    finale=len(ArrayXSec)-1

    # dizionario sez_monte
    dic_SezMonte={}
    # dizionario sez_valle
    dic_SezValle={}
    i_monte=0
    for i in range(len(ArrayXSec)):
        listax.append(ArrayXSec[i])
        listay.append(ArrayYSec[i])
        xOffset= int((ArrayXSec[i]-originX)/pixelWidth)
        yOffset= int((ArrayYSec[i]-originY)/pixelHeight)
        if xOffset>cols-1:
            xOffset=cols-1
        elif xOffset<0:
            xOffset=0
        if yOffset>rows-1:
            yOffset=rows-1
        elif yOffset<0:
            yOffset=0

        zpunto=Elevation[yOffset, xOffset]

        # aggiorna dizionari
        if i in  listaIdSecPrincipali:
            dic_SezMonte[i]=-1
            dic_SezValle[i]=-1
        else:
            id_valle=numpy.searchsorted(ArrayIdSecPrinc, i)
            dic_SezMonte[i]=ArrayIdSecPrinc[id_valle-1]
            dic_SezValle[i]=ArrayIdSecPrinc[id_valle]

        # cercare anche le quote di un certo numero di punti intorno !
        ListaZeta=[]
        for ii in range(steps):
            for jj in range(steps):
                yyy=yOffset-deltapixel+ii
                xxx=xOffset-deltapixel+jj
                try:
                    zzz= Elevation[yyy, xxx]
                    if zzz!=inNoData:
                        ListaZeta.append(zzz)
                except:
                    pass
        nn=len(ListaZeta)
        if nn>0:
            ZetaArray=numpy.array(ListaZeta)
            # prendo il valor minimo
            z=ZetaArray.min()
            # calcolo il valor medio con il punto della sezione
            z=(zpunto+z)/2.0
        else:
            z=0.0
            txt='err sez=%d' % i
            errMsg+='\n%s'% txt
        QuoteAlveo.append(z)

    # costruisco un alveo pendente
    npt=len(ArrayXSec)

    PendMin=float(0.0001)

    if QuoteAlveo[-1]<0:
        QuoteAlveo[-1]=0.0

    for i in range(npt-1,0,-1):
        dy=float(QuoteAlveo[i-1]-QuoteAlveo[i])
        dx=float(dist_along_2[i]-dist_along_2[i-1])
        pend=dy/dx
        if pend<PendMin:
            QuoteAlveo[i-1]=QuoteAlveo[i]+PendMin*dx


    # ============================
    # salva i dati nel geodatabase
    # ============================

    TargetTabella='PathPoints'

    # codice del sistema di riferimento della tabella
    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (TargetTabella.lower())
    cur.execute(sql)
    record=cur.fetchone()
    if record!=None:
        OriginEPSG=record[0]
    else:
        OriginEPSG=32632

    # controllo la esistenza di dati pregressi
    sql='SELECT PKUID,id FROM %s WHERE ID_Diga=%d' % (TargetTabella,ID_Diga)
    cur.execute(sql)
    ListaTratti=cur.fetchall()
    if len(ListaTratti)>0:
        # cancello i dati pregressi
        sql='DELETE FROM %s WHERE ID_Diga=%d' % (TargetTabella,ID_Diga)
        cur.execute(sql)
        conn.commit()

    for i in range(len(progr_Points)):
        id_cur= ListaId[i]
        progr_cur=progr_Points[id_cur]

        pt=ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(listax[i],listay[i])

        pt.FlattenTo2D()
        pt_monte=pt.ExportToWkt()

        if id_cur in listaIdSecPrincipali:
            tipo_i=1
        else:
            tipo_i=0
        zcur=float("{0:.2f}".format(QuoteAlveo[i]))

        # salvo nel geodatabase
        sql='INSERT INTO %s (ID_Diga,id,type,elev,progr,geom) VALUES (%d'  %  (TargetTabella,ID_Diga)
        sql+=',%d' %  id_cur
        sql+=',%d' %  tipo_i
        sql+=',%.2f' % zcur
        sql+=',%.3f' % progr_cur

        GeomWKT="GeomFromText('%s',%d)" % (pt_monte,OriginEPSG)
        sql+=',%s' % GeomWKT
        sql+=');'
        cur.execute(sql)

        pt.Destroy()

    # salvo
    conn.commit()

    # aggiorno gli id delle sezioni principali
    # ----------------------------------------

    # tabella delle sezioni principali
    NomeTabella='MainCrossSec'

    for PKUID in PKUID_sec_id:
        ii=PKUID_sec_id[PKUID]
        sql='UPDATE %s SET id=%d WHERE ID_Diga=%d AND PKUID=%d'  %  (NomeTabella,ii,ID_Diga,PKUID)
        cur.execute(sql)

    conn.commit()
    # Close communication with the database
    cur.close()
    conn.close()


    return NotErr, errMsg

def SetCrossSec_2(mydb_path_user,PathFiles,ClipDEM,ID_Diga):

    """
    Crea le sezioni interpolate utilizzando le sezioni principali come direttrici

    Crea inoltre anche il grid StreamDH.tif delle altezze rispetto al fiume
    """

    NotErr=bool('True')
    errMsg='OK'


    # controllo esistenza dati di input
    # ---------------------------------
    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    CrossSec=PathFiles+os.sep+'CrossSec.shp'

    if not os.path.exists(ClipDEM):
        errMsg = "Manca per la diga num =%s il clip DEM\nEffettuare prima il ritaglio del modello digitale del terreno !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg


    # =====================================
    # apertura del database sqlite utente
    # ====================================

    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
   # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    # controlla esistenza linea tratto a valle
    NomeTabellaLinee='LineaValleDiga'

    # codice del sistema di riferimento della tabella
    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
    cur.execute(sql)

    record=cur.fetchone()
    if record!=None:
        DominioEPSG=record[0]
    else:
        DominioEPSG=32632

    dest_srs = ogr.osr.SpatialReference()
    dest_srs.ImportFromEPSG(DominioEPSG)

    sql='SELECT TotalLength,ST_AsText(geom) FROM %s WHERE ID_Diga=%d' % (NomeTabellaLinee,ID_Diga)
    cur.execute(sql)
    ChkDiga=cur.fetchone()

    if ChkDiga==None:
        errMsg = "Nella tabella= %s non ci sono dati per la diga num =%s \nEffettuare prima il calcolo della linea a valle !" % (NomeTabellaLinee,ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt_line=ChkDiga[1]
        TotalLength_1=ChkDiga[0]
        # crea la geometria della linea dell'asse del fiume
        StreamLine=ogr.CreateGeometryFromWkt(wkt_line)
        TotalLength=StreamLine.Length()

    # creo poly_stream di 50 m
    poly_stream=StreamLine.Buffer(50)

    TabellaPoligono='PoligonoValleDiga'
    sql='SELECT Area,ST_AsText(geom) FROM %s WHERE ID_Diga=%d' % (TabellaPoligono,ID_Diga)
    cur.execute(sql)
    ChkDigaPoly=cur.fetchone()

    if ChkDigaPoly==None:
        errMsg = "Nella tabella= %s non ci sono dati per la diga num =%s \nEffettuare prima il calcolo della zona di studio a valle !" % (TabellaPoligono,ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt_Dominio=ChkDigaPoly[1]
        TotalArea=ChkDigaPoly[0]

    CrossSecTot='%sTot.shp' % (CrossSec[:-4])
    CrossSecPoly='%sPoly.shp' % (CrossSec[:-4])

    # linea mediana
    LineaMediana='%sMedianPath.shp' % (CrossSec[:-4])

    # poligoni divisi in destra e sinistra
    CrossSecPoly2='%sPoly_2.shp' % (CrossSec[:-4])

    # traccia sezioni intermedie rappresentative del tratto
    CrossMedie='%sMean.shp' % (CrossSec[:-4])


    # traccia sezioni interpolate e sezioni intermedie rappresentative del tratto
    # sevono per creare StreamDH.tif meglio definito !!
    Cross_Tot_Medie='%sTot_and_Mean.shp' % (CrossSec[:-4])

    # Punto sul fiume della sezione rappresentativa del tratto
    CrossMediePunto='%sMeanPunto.shp' % (CrossSec[:-4])

    # poligoni destra e sinistra fiume
    PolySxDx=PathFiles + os.sep +'PolySxDx.shp'

    # carico il driver dello shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # creo il Multipolygon del dominio
    DominoPoly=ogr.CreateGeometryFromWkt(wkt_Dominio)

    Area=DominoPoly.Area()

    # prende il poligono piu' grande nei Multipoligons
    numpoly=DominoPoly.GetGeometryCount()

    if numpoly>1:
        Amax=0.0
        ii_max=0
        for ii in range(numpoly):
            g=DominoPoly.GetGeometryRef(ii)
            Aii=g.Area()
            if  Aii>Amax:
                Amax=Aii
                ii_max=ii
        DominoRing=DominoPoly.GetGeometryRef(ii_max)
        Area1=DominoRing.Area()

    else:
        # prende il primo
        DominoRing=DominoPoly.GetGeometryRef(0)
        Area1=DominoRing.Area()

    # leggo la geometria delle sezioni principali
    # ---------------------------------------------
    # tabella delle sezioni principali
    NomeTabella='MainCrossSec'


    sql='SELECT id,ST_AsText(geom) FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
    cur.execute(sql)
    ListaCrossSec=cur.fetchall()

    n=len(ListaCrossSec)

    if n <= 0:

        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni principali !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    else:

        # codice del sistema di riferimento della tabella
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        cur.execute(sql)
        record=cur.fetchone()
        if record!=None:
            OriginEPSG=record[0]
        else:
            OriginEPSG=32632

        # carico le geometrie in un dizionario
        CrossSecGeom={}

        for rec in ListaCrossSec:
            CrossSecGeom[rec[0]]=rec[1]

        # lettura punti di interpolazione
        # --------------------------------

        # tipo di sezione i-esima
        NumSezTipo={}

        # Caratteristiche della sezione
        NumSezElev={}
        NumSezProgr={}
        NumSezGeom={}

        # numero di sezione a monte di i-esima
        NumSezMonte={}
        # numero della sezione a valle di una sezione di monte
        NumSezValle={}

        # numero di sezione di monte corrente
        NumMonte=-1

        # creo il dizionario dei punti a monte
        dic_pt_Monte={}
        NumPt_Monte=0


        TabellaPoints='PathPoints'

        # codice del sistema di riferimento della tabella
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        cur.execute(sql)
        record=cur.fetchone()
        if record!=None:
            SourceEPSG=record[0]
        else:
            SourceEPSG=32632

       # controllo la esistenza di dati pregressi
        sql='SELECT id,type,elev,progr,ST_AsText(geom) FROM %s WHERE ID_Diga=%d ORDER BY id' % (TabellaPoints,ID_Diga)
        cur.execute(sql)
        ListaTratti=cur.fetchall()

        n=len(ListaTratti)

        if n <= 0:

            errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni principali !" % (ID_Diga)
            NotErr= bool()
            return NotErr, errMsg

        for rec in ListaTratti:

            NumSez=rec[0]
            tipo=rec[1]
            NumSezTipo[NumSez]=rec[1]
            NumSezElev[NumSez]=rec[2]
            NumSezProgr[NumSez]=rec[3]
            NumSezGeom[NumSez]=rec[4]

            # aggiorna NumMonte
            if tipo==1:
                if NumMonte>=0:
                    NumSezValle[NumMonte]=NumSez
                    NumSezMonte[NumSez]=NumMonte
                    NumMonte=NumSez
                else:
                    NumMonte=NumSez
            else:
                NumSezMonte[NumSez]=NumMonte

            dic_pt_Monte[NumSez]=NumPt_Monte
            NumPt_Monte=NumSez


    # Close communication with the database
    cur.close()
    conn.close()

    # ===============================
    # creazione sezioni interpolate
    # ===============================

    # -----------------------------------------
    # Crea lo shape file delle sezioni totali
    # -----------------------------------------

    # crea nuovo shapefile
    if os.path.exists(CrossSecTot):
        driver.DeleteDataSource(CrossSecTot)

    outDS2 = driver.CreateDataSource(CrossSecTot)
    if outDS2 is None:
        errMsg='Could not create file ' + CrossSecTot
        NotErr= bool()
        return NotErr, errMsg

    outLayer2 = outDS2.CreateLayer('CrossSecTot', dest_srs,geom_type=ogr.wkbLineString)

    # crea i nuovi campi id e Nome nello shapefile di output
    fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
    outLayer2.CreateField(fieldDefn2)
    fieldDefn3 = ogr.FieldDefn('type', ogr.OFTInteger)
    outLayer2.CreateField(fieldDefn3)
    fieldDefn4 = ogr.FieldDefn('zmin', ogr.OFTReal)
    outLayer2.CreateField(fieldDefn4)
    fieldDefn5 = ogr.FieldDefn('progr', ogr.OFTReal)
    outLayer2.CreateField(fieldDefn5)


    featureDefn2 = outLayer2.GetLayerDefn()

    bufferDistance = 50

    # lista delle sezioni salvate
    # viene utilizzate per la creazione congruente degli altri dati
    ListaNumSezTotOk=[]


     # tratto di fiume del punto della sezione
    dic_Sez_tratto_fiume={}

    for NumSez in NumSezTipo:

        AddFeature=bool()
        tipo=NumSezTipo[NumSez]
        z=NumSezElev[NumSez]
        Progr=NumSezProgr[NumSez]

        pt_curr=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

        pt_poly = pt_curr.Buffer(bufferDistance)

        poly_trattocorrente=pt_poly.Intersection(poly_stream)
        dic_Sez_tratto_fiume[NumSez]= poly_trattocorrente.ExportToWkt()


        feature = ogr.Feature(featureDefn2)
        feature.SetField('id', NumSez)
        feature.SetField('type',tipo)
        feature.SetField('zmin',z)
        feature.SetField('progr',Progr)


        try:
            if tipo==0:

##                txt='Sez interp %d : monte = %d - valle = %d' % (NumSez,NumSezMonte[NumSez],NumSezValle[NumSezMonte[NumSez]])
##                print (txt)
                numMonte=NumSezMonte[NumSez]
                SezMonteGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[numMonte])
                wktValle=CrossSecGeom[NumSezValle[numMonte]]
                SezValleGeom=ogr.CreateGeometryFromWkt(wktValle)
##                NewGeom = SezioneInterpolata(NumSez,pt_curr,SezMonteGeom,SezValleGeom)
                NewGeom = SezioneInterpolata(pt_curr,SezMonteGeom,SezValleGeom)

                try:
                    lengh1=NewGeom.Length()

##                    wkt_debug1=NewGeom.ExportToWkt()

                    # taglio lungo il dominio
                    NewGeom1=NewGeom.Intersection(DominoRing)

##                    wkt_debug2=NewGeom1.ExportToWkt()

                    lengh2=NewGeom1.Length()

                    if lengh2>0:
                        if NewGeom1.GetGeometryCount()>0:
                            for i in range(0, NewGeom1.GetGeometryCount()):
                                line=NewGeom1.GetGeometryRef(i)
                                if line.Intersect(pt_poly):
                                    feature.SetGeometry(line)
                                    lengh3=line.Length()
                                    if lengh3>0:
                                        AddFeature=bool('True')
                                    break
                        else:
                            AddFeature=bool('True')
                            feature.SetGeometry(NewGeom1)
                    else:
                        pass

                except:
                    pass

            else:

##                txt=' === Sezione principale  %s  ===' % NumSez
##                print (txt)
                NewGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[NumSez])

                try:
                    lengh1=NewGeom.Length()

                    NewGeom1=NewGeom.Intersection(DominoRing)

                    lengh2=NewGeom1.Length()

                    if lengh2>0:
                        if NewGeom1.GetGeometryCount()>0:
                            for i in range(0, NewGeom1.GetGeometryCount()):
                                line=NewGeom1.GetGeometryRef(i)
                                if line.Intersect(pt_poly):
                                    lengh3=line.Length()
                                    if lengh3>0:
                                        AddFeature=bool('True')
                                        feature.SetGeometry(line)
                                    break
                        else:
                            lengh3=NewGeom1.Length()
                            if lengh3>0:
                                AddFeature=bool('True')
                                feature.SetGeometry(NewGeom1)
                    else:
                        pass

                except:
                    pass
            if AddFeature:
                outLayer2.CreateFeature(feature)
                ListaNumSezTotOk.append(NumSez)
            else:
                txt=' === Error Sezione  %s  ===' % NumSez
                errMsg+='\n%s'% txt
            NewGeom.Destroy()

        except:
            txt=' === Error Sezione  %s  ===' % NumSez
            errMsg+='\n%s'% txt


    outDS2.Destroy()

    # carica le geometrie delle sezioni da clippare
    # ---------------------------------------------
    # sorgente dei dati in lettura
    ds = driver.Open(CrossSecTot, 0)
    if ds is None:
        errMsg='Could not open file %s' % CrossSecTot
        NotErr= bool()
        return NotErr, errMsg


    # leggo il layer dalla sorgente dei dati
    layer = ds.GetLayer()

    # conto il numero di Feature
    n = layer.GetFeatureCount()

    # dizionario geometrie
    dic_Sez_ToClip={}

    feat = layer.GetNextFeature()

    while feat:

        NumSez=feat.GetField('id')
        line_sez=feat.GetGeometryRef()
        wkt=line_sez.ExportToWkt()
        dic_Sez_ToClip[NumSez]=wkt

        feat = layer.GetNextFeature()

    # chiudo la connessione
    ds.Destroy()


    # controllo intersezione sezioni totali
    # --------------------------------------

    bufferClip=50.0

    distBufferClip=bufferClip*1.1

    #sorgente dei dati in scrittura
    ds = driver.Open(CrossSecTot, 1)
    if ds is None:
        errMsg='Could not open file %s' % CrossSecTot
        NotErr= bool()
        return NotErr, errMsg


    # leggo il layer dalla sorgente dei dati
    layer = ds.GetLayer()

    # conto il numero di Feature
    n = layer.GetFeatureCount()

    # dizionario dove salvato il punto di intersezione
    dic_Sez_PointIntersect={}

    # dizionario dove salvato il buffer
    dic_Sez_Clip={}
    # dizionario se il buffer e' a destra
    dic_Sez_lato_destro={}

    # Lista sezioni clippate al primo ciclo  while
    lista_clipped=[]

    # Lista sezioni di valle intersecate al primo ciclo  while
    lista_intersects=[]

    # leggo la prima feature
    feat = layer.GetNextFeature()

    num_sez=len(ListaNumSezTotOk)

    while feat:

        NumSez_Monte=feat.GetField('id')
        line_monte_geom=feat.GetGeometryRef()
        ii_ini=ListaNumSezTotOk.index(NumSez_Monte)+1
        # tratto corrente di linea di fiume
        line_trattocorrente = ogr.CreateGeometryFromWkt(dic_Sez_tratto_fiume[NumSez_Monte])

        update_sez=bool()

        for ii in range(ii_ini,num_sez):

            NumSez_Valle=ListaNumSezTotOk[ii]
            wkt=dic_Sez_ToClip[NumSez_Valle]
            line_valle= ogr.CreateGeometryFromWkt(wkt)


            # controllo se le due linee non si intersecano
            # ---------------------------------------------
            if line_valle.Crosses(line_monte_geom):

                update_sez=bool('True')

                # trova il punto di intersezione
                pt_interes_sez=line_monte_geom.Intersection(line_valle)

                # controllo se a sinistra o destra
                # --------------------------------

                # controllo se a destra o sinistra del fiume
                orientamento=ProdottoVettSezioni(line_monte_geom,line_valle)
                if  orientamento>0:
                    # si presume una intersezione in sponda sinistra
                    inters_destra=bool()
                else:
                    inters_destra=bool('True')

                # crea un buffer intorno al punto
                pt_poly = pt_interes_sez.Buffer(bufferClip)

                if NumSez_Valle not in  lista_intersects:
                    lista_intersects.append(NumSez_Valle)
                # salva il buffer e eventualmente sovrascrive
                dic_Sez_Clip[NumSez_Valle]= pt_poly.ExportToWkt()
                dic_Sez_PointIntersect[NumSez_Valle]= pt_interes_sez.ExportToWkt()
                dic_Sez_lato_destro[NumSez_Valle]= inters_destra

                # esegue il clip di solo quella di monte
                # la linea viene accorciata in modo da non toccare la sezione di valle
                lines=line_monte_geom.Difference(pt_poly)

                # conta il numero di geometrie: vi possono essere casi in cui eliminando
                # dalla linea la parte intersecata dal buffer si abbiano come risultato
                # due linee di cui una  a destra ed una a sinistra
                nn=lines.GetGeometryCount()
                if nn>1:
                    for i in range(nn):
                        g = lines.GetGeometryRef(i)
                        # sceglie la parte di linea che interseca l'asse del fiume
                        if g.Crosses(line_trattocorrente):
                            line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            # controllare il verso
                            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                            break
                else:
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

                wkt_monte_new=line.ExportToWkt()
                line_monte_geom= ogr.CreateGeometryFromWkt(wkt_monte_new)

        # salva le modifiche nello shapefile
        if update_sez:
            feat.SetGeometry(line)
            layer.SetFeature(feat)

        # leggo la feature successiva
        feat = layer.GetNextFeature()

    # esegue l'eventuale clip dell'ultima sezione
    if NumSez_Monte in lista_intersects:

        wkt=dic_Sez_ToClip[NumSez_Monte]
        line_valle= ogr.CreateGeometryFromWkt(wkt)

        # tratto corrente di linea di fiume
        line_trattocorrente = ogr.CreateGeometryFromWkt(dic_Sez_tratto_fiume[NumSez])

        pt_interes_sez= ogr.CreateGeometryFromWkt(dic_Sez_PointIntersect[NumSez])

        # controllo se a sinistra o destra
        # --------------------------------
        inters_destra= dic_Sez_lato_destro[NumSez]

        pt_poly=ogr.CreateGeometryFromWkt(dic_Sez_Clip[NumSez])

        # esegue il clip della sezione corrente
        lines=line_valle.Difference(pt_poly)

        # conta il numero di geometrie
        nn=lines.GetGeometryCount()
        if nn>1:

            for i in range(nn):
                g = lines.GetGeometryRef(i)
                if g.Crosses(line_trattocorrente):
                    line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                    # controllare il verso
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                    break
        else:
            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

        feat.SetGeometry(line)

        layer.SetFeature(feat)

    # chiudo la connessione
    ds.Destroy()


    # salvo la linea dell'ultima sezione
    # ----------------------------------
    #sorgente dei dati in lettura
    ds = driver.Open(CrossSecTot, 0)
    if ds is None:
        errMsg='Could not open file %s' % CrossSecTot
        NotErr= bool()
        return NotErr, errMsg

    # leggo il layer dalla sorgente dei dati
    layer = ds.GetLayer()

    # conto il numero di Feature
    n = layer.GetFeatureCount()
    feat_end = layer.GetFeature(n-1)
    NumSez_end=feat_end.GetField('id')
    line_valle=feat_end.GetGeometryRef()
    wkt_sez_end=line_valle.ExportToWkt()

    # chiudo la connessione
    ds.Destroy()

    # Crea anche il dizionario delle coordinate dei punti
    dic_PointPath={}

    for NumPoint in NumSezTipo:

        pt_curr=ogr.CreateGeometryFromWkt(NumSezGeom[NumPoint])

        punto=pt_curr.GetPoint()

         # aggiunge al dizionario
        dic_PointPath[NumPoint]=punto

    # ---------------------------------------------------------
    # crea i poligoni dei tratti  di fiume
    # -------------------------------------
    # e creazione del poligono globale verificando l'intersezione
    # ---------------------------------------------------------

    CreaPolygoni=1

    if CreaPolygoni>0:

        # salva Dizionario con la lista dei punti poligono
        DicPuntiPoligono={}

        # salva Dizionario con la geometria wkt del poligono finale clippato
        DicPoligonoClip={}

        # salva Dizionario con la lista dei punti dell'asse del fiume
        DicPuntiAsse={}

        # salva la lista dei numeri dei poligoni
        NumPolyList=[]

        # Creo il poligono unione dei tratti da utilizzare per il clip
        # --------------------------------------------------------------
        # inizializzo il poligono dei tratti
        PoligonoUnioneTratti_0=ogr.Geometry(ogr.wkbPolygon)

        #sorgente dei dati in lettura
        ds = driver.Open(CrossSecTot, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossSecTot
            NotErr= bool()
            return NotErr, errMsg


        # leggo il layer dalla sorgente dei dati
        Inlayer = ds.GetLayer()
        Spatialref = Inlayer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        SourceEPSG=int(Spatialref.GetAuthorityCode(None))


        # crea nuovo shapefile
        shpnew3=CrossSecPoly
        if os.path.exists(shpnew3):
            driver.DeleteDataSource(shpnew3)

        outDS3 = driver.CreateDataSource(shpnew3)
        if outDS3 is None:
            errMsg='Could not create file %s' % shpnew3
            NotErr= bool()
            return NotErr, errMsg

        outLayer3 = outDS3.CreateLayer('CrossSecPoly', Spatialref,geom_type=ogr.wkbPolygon)

        # crea i nuovi campi id e Nome nello shapefile di output
        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer3.CreateField(fieldDefn2)
        fieldDefn4 = ogr.FieldDefn('zmin', ogr.OFTReal)
        outLayer3.CreateField(fieldDefn4)

        featureDefn2 = outLayer3.GetLayerDefn()

        feat = Inlayer.GetNextFeature()
        NumMonte=feat.GetField('id')
        line_monte=feat.GetGeometryRef()
        line_monte_wkt=line_monte.ExportToWkt()

        while feat:

            feat = Inlayer.GetNextFeature()

            if feat:

                NumSValle=feat.GetField('id')
                Zvalle=float(feat.GetField('zmin'))

                line_valle=feat.GetGeometryRef()
                wkt=line_valle.ExportToWkt()

                feature = ogr.Feature(featureDefn2)
                feature.SetField('id', NumSValle)
                feature.SetField('zmin', Zvalle)

                # crea il poligono
                ring = ogr.Geometry(ogr.wkbLinearRing)
                line_monte=ogr.CreateGeometryFromWkt(line_monte_wkt)

                # controllo se le due linee non si intersecano
                # ---------------------------------------------
                if not line_monte.Crosses(line_valle):


                    ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                    ring.AddPoint(line_monte.GetX(1),line_monte.GetY(1))
                    ring.AddPoint(line_valle.GetX(1),line_valle.GetY(1))
                    ring.AddPoint(line_valle.GetX(0),line_valle.GetY(0))
                    ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))

                    # salvo i punti
                    PuntiPoligono=[]
                    PuntiPoligono.append((line_monte.GetX(0),line_monte.GetY(0)))
                    PuntiPoligono.append((line_monte.GetX(1),line_monte.GetY(1)))
                    PuntiPoligono.append((line_valle.GetX(1),line_valle.GetY(1)))
                    PuntiPoligono.append((line_valle.GetX(0),line_valle.GetY(0)))
                    DicPuntiPoligono[NumSValle]=PuntiPoligono

                else:

                    pt_interes_sez=line_monte.Intersection(line_valle)

                    # controllo distanza a sinistra
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                    dist_sx=pt_interes_sez.Distance(point)
                    # controllo distanza a destra
                    point.Destroy()
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(line_monte.GetX(1),line_monte.GetY(1))
                    dist_dx=pt_interes_sez.Distance(point)
                    point.Destroy()

                    if dist_dx>dist_sx:
                        # caso intersezione a sinistra
                        ring.AddPoint(pt_interes_sez.GetX(0),pt_interes_sez.GetY(0))
                        ring.AddPoint(line_monte.GetX(1),line_monte.GetY(1))
                        ring.AddPoint(line_valle.GetX(1),line_valle.GetY(1))
                        ring.AddPoint(pt_interes_sez.GetX(0),pt_interes_sez.GetY(0))
                        # salvo i punti
                        PuntiPoligono=[]
                        PuntiPoligono.append((pt_interes_sez.GetX(0),pt_interes_sez.GetY(0)))
                        PuntiPoligono.append((line_monte.GetX(1),line_monte.GetY(1)))
                        PuntiPoligono.append((line_valle.GetX(1),line_valle.GetY(1)))
                        PuntiPoligono.append((pt_interes_sez.GetX(0),pt_interes_sez.GetY(0)))
                        DicPuntiPoligono[NumSValle]=PuntiPoligono

                    else:
                        # caso intersezione a destra
                        ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                        ring.AddPoint(pt_interes_sez.GetX(0),pt_interes_sez.GetY(0))
                        ring.AddPoint(line_valle.GetX(0),line_valle.GetY(0))
                        ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                        # salvo i punti
                        PuntiPoligono=[]
                        PuntiPoligono.append((line_monte.GetX(0),line_monte.GetY(0)))
                        PuntiPoligono.append((pt_interes_sez.GetX(0),pt_interes_sez.GetY(0)))
                        PuntiPoligono.append((line_valle.GetX(0),line_valle.GetY(0)))
                        PuntiPoligono.append((line_monte.GetX(0),line_monte.GetY(0)))
                        DicPuntiPoligono[NumSValle]=PuntiPoligono

                    pt_interes_sez.Destroy()

                # salvo il numero nella lista
                NumPolyList.append(NumSValle)

                poly_new=ogr.Geometry(ogr.wkbPolygon)

                poly_new.AddGeometry(ring)

                numpt1=poly_new.GetGeometryRef(0).GetPointCount()

                # interseca il poligono del dominio
                NewPoly=poly_new.Intersection(DominoRing)

                if NewPoly!=None:
                    npoly=NewPoly.GetGeometryCount()
                    if  npoly>0:
                        numpt=NewPoly.GetGeometryRef(0).GetPointCount()
                        # intersezione con l'asse del fiume
                        TrattoStreamLine= StreamLine.Intersection(NewPoly)
                        nntratti=TrattoStreamLine.GetGeometryCount()
                        if nntratti>0:
                            # per evitare problemi di tratti con anse uso la linea retta che unisce
                            # il punto di monte con quello di valle
                            PuntiTratto=[]
                            PuntiTratto.append(dic_PointPath[NumMonte])
                            PuntiTratto.append(dic_PointPath[NumSValle])
                        else:
                            TrattoStream_max=TrattoStreamLine
                            # estraggo i punti
                            PuntiTratto=[]
                            for iii in range(0, TrattoStream_max.GetPointCount()):
                                # GetPoint returns a tuple not a Geometry
                                pt = TrattoStream_max.GetPoint(iii)
                                PuntiTratto.append(pt)

                        DicPuntiAsse[NumSValle]=PuntiTratto

                        poly_wkt=NewPoly.ExportToWkt()
                    else:
                        # geometria assente !!
                        poly_wkt=NewPoly.ExportToWkt()
                        # uso quella prima della intersezione
                        NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                        npoly=poly_new.GetGeometryCount()
                        geom_poly=poly_new.GetGeometryRef(0)
                        numpt=geom_poly.GetPointCount()

                        poly_wkt=NewPoly.ExportToWkt()
                        # salvo i dati
                        nome_prova=PathFiles +os.sep+ 'File_poly_%d.csv' % NumSValle
                        fout=open(nome_prova,'w')
                        txt='X;Y\n'
                        fout.write(txt)
                        for iii in range(0, numpt):
                            # GetPoint returns a tuple not a Geometry
                            pt = geom_poly.GetPoint(iii)
                            txt='%s;%s\n' %(pt[0],pt[1])
                            fout.write(txt)
                        fout.close()


                        # per evitare problemi di tratti con anse uso la linea retta che unisce
                        # il punto di monte con quello di valle
                        PuntiTratto=[]
                        PuntiTratto.append(dic_PointPath[NumMonte])
                        PuntiTratto.append(dic_PointPath[NumSValle])

                        DicPuntiAsse[NumSValle]=PuntiTratto

                else:

                    # casi di mancata intersezione ??????
                    NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                    npoly=poly_new.GetGeometryCount()
                    geom_poly=poly_new.GetGeometryRef(0)
                    numpt=geom_poly.GetPointCount()

                    poly_wkt=NewPoly.ExportToWkt()
                    # salvo i dati
                    nome_prova=PathFiles +os.sep+ 'File_poly_%d.csv' % NumSValle
                    fout=open(nome_prova,'w')
                    txt='X;Y\n'
                    fout.write(txt)
                    for iii in range(0, numpt):
                        # GetPoint returns a tuple not a Geometry
                        pt = geom_poly.GetPoint(iii)
                        txt='%s;%s\n' %(pt[0],pt[1])
                        fout.write(txt)
                    fout.close()


                    # per evitare problemi di tratti con anse uso la linea retta che unisce
                    # il punto di monte con quello di valle
                    PuntiTratto=[]
                    PuntiTratto.append(dic_PointPath[NumMonte])
                    PuntiTratto.append(dic_PointPath[NumSValle])

                    DicPuntiAsse[NumSValle]=PuntiTratto

                # controllo intersezione con le parti a monte
                NewPolyDiff=NewPoly.Difference(PoligonoUnioneTratti_0)

                # conta il numero di geometrie
                nn=NewPolyDiff.GetGeometryCount()
                if nn>1:
                    # tratto corrente di linea di fiume
                    poly_centro = ogr.CreateGeometryFromWkt(dic_Sez_tratto_fiume[NumMonte])
                    for i in range(nn):
                        g = NewPolyDiff.GetGeometryRef(i)
                        wkt_debug= g.ExportToWkt()

                        # sceglie la parte che interseca l'asse del fiume
                        if g.Intersect(poly_centro):
                            NewPoly_0=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            NewPoly=NewPoly_0.Buffer(0.1)
                            break
                else:
                    NewPoly=NewPolyDiff.Buffer(0.1)

                # unione
                nn0=PoligonoUnioneTratti_0.GetGeometryCount()

                PoligonoUnioneTratti_0=PoligonoUnioneTratti_0.Union(NewPoly)

                # controlla i casi con buchi
                nn1=PoligonoUnioneTratti_0.GetGeometryCount()
                if nn1>1:
                    type_name=PoligonoUnioneTratti_0.GetGeometryName()
                    if type_name=="POLYGON":
                        area_max=0.0
                        for k in range(nn1):
                            ring1=PoligonoUnioneTratti_0.GetGeometryRef(k)
                            are_curr=ring1.GetArea()
                            if are_curr>area_max:
                                k_max=k*1
                                area_max=are_curr
                        ring_max=PoligonoUnioneTratti_0.GetGeometryRef(k_max)
                        npt_max= ring_max.GetPointCount()
                        ring_new = ogr.Geometry(ogr.wkbLinearRing)
                        for i in range(npt_max):
                            pt_curr=ring_max.GetPoint(i)
                            ring_new.AddPoint(pt_curr[0],pt_curr[1])
                        PoligonoUnioneTratti_0=None
                        PoligonoUnioneTratti_0=ogr.Geometry(ogr.wkbPolygon)
                        PoligonoUnioneTratti_0.AddGeometry(ring_new)

                    elif type_name=="MULTIPOLYGON":
                        area_max=0.0
                        for k in range(nn1):
                            poli_cur=PoligonoUnioneTratti_0.GetGeometryRef(k)
                            are_curr=poli_cur.GetArea()
                            if are_curr>area_max:
                                k_max=k*1
                                area_max=are_curr
                                wkt_max=poli_cur.ExportToWkt()
                        PoligonoUnioneTratti_0=None
                        PoligonoUnioneTratti_0=ogr.CreateGeometryFromWkt(wkt_max)

                # salva il poligono clippato
                DicPoligonoClip[NumSValle]= NewPoly.ExportToWkt()

                feature.SetGeometry(NewPoly)
                outLayer3.CreateFeature(feature)
                ring.Destroy()
                poly_new.Destroy()
                NewPoly.Destroy()

                # aggiorno per lo step successivo

                line_monte_wkt=wkt

                NumMonte=NumSValle

        # chiudo la connessione
        ds.Destroy()

        outDS3.Destroy()

        PoligonoUnioneTratti_Wkt=PoligonoUnioneTratti_0.ExportToWkt()


    # Shapefile sezioni mediane
    CreaSezioniMediane=1

    if CreaSezioniMediane>0:

        # si presuppone di aver gia' memorizzato NumSezMonte e NumSezValle
        n = len(ListaTratti)

        # salva il buffer del buffer del punto di intersezione delle sezioni
        # mediane con l'asse del fiume
        dic_SezMediane_Clip={}


        # crea nuovo shapefile dei punti sul fiume delle sezioni
        # ------------------------------------------------------
        shpnew2=CrossMediePunto
        if os.path.exists(shpnew2):
            driver.DeleteDataSource(shpnew2)

        outDS2 = driver.CreateDataSource(shpnew2)
        if outDS2 is None:
            errMsg='Could not create file %s' % shpnew2
            NotErr= bool()
            return NotErr, errMsg

        outLayer2 = outDS2.CreateLayer('PuntoMedio', dest_srs,geom_type=ogr.wkbPoint)

        # crea nuovo shapefile delle sezioni
        # ----------------------------------
        shpnew3=CrossMedie
        if os.path.exists(shpnew3):
            driver.DeleteDataSource(shpnew3)

        outDS3 = driver.CreateDataSource(shpnew3)
        if outDS3 is None:
            errMsg='Could not create file %s' % shpnew3
            NotErr= bool()
            return NotErr, errMsg

        outLayer3 = outDS3.CreateLayer('SezioneMedia', dest_srs,geom_type=ogr.wkbLineString)

        # crea i nuovi campi id e Nome nello shapefile di output
        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer2.CreateField(fieldDefn2)
        outLayer3.CreateField(fieldDefn2)
        # crea campo distanza punto medio
        fieldDefnDist = ogr.FieldDefn('dist1', ogr.OFTReal)
        outLayer3.CreateField(fieldDefnDist)
        fieldDefnProgr = ogr.FieldDefn('progr', ogr.OFTReal)
        outLayer3.CreateField(fieldDefnProgr)

        fieldDefnZmin = ogr.FieldDefn('zmin', ogr.OFTReal)
        outLayer3.CreateField(fieldDefnZmin)

        fieldDefnWL = ogr.FieldDefn('hmax', ogr.OFTReal)
        outLayer3.CreateField(fieldDefnWL)

        featureDefn2 = outLayer2.GetLayerDefn()
        featureDefn3 = outLayer3.GetLayerDefn()

        # leggo la prima feature
        NumSez=ListaTratti[0][0]
        ProgrMonte=NumSezProgr[NumSez]
        Z_monte= NumSezElev[NumSez]

        pt_monte=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

        # creo il primo punto
        feature2 = ogr.Feature(featureDefn2)
        feature2.SetField('id', NumSez)
        feature2.SetGeometry(pt_monte)
        outLayer2.CreateFeature(feature2)

        # creo la prima sezione
        # ---------------------
        NewGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[NumSez])
        feature3 = ogr.Feature(featureDefn3)
        feature3.SetField('id', NumSez)
        feature3.SetField('progr', ProgrMonte)

        feature3.SetField('zmin',Z_monte)

        pt_sinistra=NewGeom.GetPoint(0)
        Geom_pt_sinistra=ogr.Geometry(ogr.wkbPoint)
        Geom_pt_sinistra.AddPoint(pt_sinistra[0],pt_sinistra[1])
        distanza=Geom_pt_sinistra.Distance(pt_monte)
        feature3.SetField('dist1', distanza)

        # taglio lungo il dominio
        pt1=pt_monte.GetPoint()
        pt_poly= pt_monte.Buffer(bufferClip)

        # Taglio lungo il poligono del primo tratti
        Poly_curr=ogr.CreateGeometryFromWkt(DicPoligonoClip[ListaTratti[1][0]])
        Poly_curr=Poly_curr.Buffer(0.5)
        NewGeom=NewGeom.Intersection(Poly_curr)

        lengh2=NewGeom.Length()

        if lengh2>0:
            if NewGeom.GetGeometryCount()>0:
                for i in range(0, NewGeom.GetGeometryCount()):
                    line=NewGeom.GetGeometryRef(i)
                    if line.Intersect(pt_poly):
                        feature3.SetGeometry(line)
                        break
            else:
                feature3.SetGeometry(NewGeom)
        else:
            feature3.SetGeometry(NewGeom)

        outLayer3.CreateFeature(feature3)


        # salva il buffer
        dic_SezMediane_Clip[NumSez]= pt_poly.ExportToWkt()

        numtratti=len(ListaTratti)

        for itratto in range(1,numtratti):

            NumSez=ListaTratti[itratto][0]

            # salta le sezioni non salvate
            if NumSez in ListaNumSezTotOk:

                tipo=NumSezTipo[NumSez]
                ProgrValle=NumSezProgr[NumSez]
                Z_valle= NumSezElev[NumSez]

                pt_valle=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

                # crea il punto intermedio per la sezione
                pt2=pt_valle.GetPoint()
                pt_mean=PuntoMedio(pt1,pt2)
                pt_curr=ogr.Geometry(ogr.wkbPoint)
                pt_curr.AddPoint(pt_mean[0],pt_mean[1])

                pt_poly = pt_curr.Buffer(bufferClip)
                # salva il buffer
                dic_SezMediane_Clip[NumSez]= pt_poly.ExportToWkt()


                feature2 = ogr.Feature(featureDefn2)
                feature2.SetField('id', NumSez)

                feature3 = ogr.Feature(featureDefn3)
                feature3.SetField('id', NumSez)
                ProgrMedia=(ProgrMonte+ProgrValle)/2.0
                feature3.SetField('progr', ProgrMedia)

                Z_media= (Z_monte+Z_valle)/2.0
                feature3.SetField('zmin',Z_media)


                if tipo==1:
                    # caso a monte di una sezione principale
                    SezValleGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[NumSez])
                    # trova la sezione di monte
                    numMonte=NumSezMonte[NumSez]
                    SezMonteGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[numMonte])

                else:
                   # caso di una sezione intermedia
                    numMonte=NumSezMonte[NumSez]
                    SezMonteGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[numMonte])
                    wktValle=CrossSecGeom[NumSezValle[numMonte]]
                    SezValleGeom=ogr.CreateGeometryFromWkt(wktValle)

                # calcola la sezione interpolate
                try:
                    NewGeom = SezioneInterpolata(NumSez,pt_curr,SezMonteGeom,SezValleGeom)

                    try:
                        lengh1=NewGeom.Length()

                        # taglio lungo il poligono del tratto corrente
                        Poly_curr=ogr.CreateGeometryFromWkt(DicPoligonoClip[NumSez])
                        NewGeom=NewGeom.Intersection(Poly_curr)

                        lengh2=NewGeom.Length()

                        if lengh2>0:
                            if NewGeom.GetGeometryCount()>0:
                                for i in range(0, NewGeom.GetGeometryCount()):
                                    line=NewGeom.GetGeometryRef(i)
                                    if line.Intersect(pt_poly):
                                        feature3.SetGeometry(line)
                                        break
                            else:
                                feature3.SetGeometry(NewGeom)
                        else:
                            feature3.SetGeometry(NewGeom)

                    except:

                        feature3.SetGeometry(NewGeom)

                    # calcola distanza del punto medio della sezione rispetto
                    # al punto a sinistra (serve per la successiva interpolazione
                    # di punti lungo la linea della sezione stessa)
                    pt_sinistra=NewGeom.GetPoint(0)
                    Geom_pt_sinistra=ogr.Geometry(ogr.wkbPoint)
                    Geom_pt_sinistra.AddPoint(pt_sinistra[0],pt_sinistra[1])

                    # trova il punto di intersezione con l'asse del fiume
                    pt_intersect= NewGeom.Intersection(StreamLine)
                    if pt_intersect!=None:
                        # controlla il numero di intersezioni
                        n_inters=pt_intersect.GetGeometryCount()
                        # in caso di piu' punti sceglie quello piu' vicino al medio
                        if n_inters>0:
                            distgg=NewGeom.Length()
                            for iii in range(0, n_inters):
                                g = pt_intersect.GetGeometryRef(iii)
                                distgg1=g.Distance(pt_curr)
                                if distgg1<distgg:
                                    distgg=distgg1*1.0
                                    ptt=g
                        else:
                            ptt=pt_intersect
                    else:
                        ptt=pt_curr

                    # distanza dal punto di intersezione con l'asse del fiume
                    distanza=Geom_pt_sinistra.Distance(ptt)
                    feature3.SetField('dist1', distanza)

                    feature2.SetGeometry(pt_curr)
                    outLayer2.CreateFeature(feature2)

                    outLayer3.CreateFeature(feature3)

                    NewGeom.Destroy()
                    pt_curr.Destroy()

                except:
                    pass

                # passo al punto successivo
                pt1=pt2
                ProgrMonte=ProgrValle*1.0
                Z_monte=Z_valle*1.0

            else:
##                print ('Esclusa sez: %d' % NumSez )
                txt='Esclusa sez: %d' % NumSez
                errMsg+='\n%s'% txt

        # aggiungo l'ultima sezione
        # =========================

        # leggo l'ultima feature
        NumSez=ListaTratti[-1][0]
        # controllo  e' fra le sezioni salvate
        if NumSez in ListaNumSezTotOk:

            ProgrMonte=NumSezProgr[NumSez]
            Z_monte= NumSezElev[NumSez]

            pt_monte=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

            # creo l'ultimo punto
            feature2 = ogr.Feature(featureDefn2)
            feature2.SetField('id', NumSez)
            feature2.SetGeometry(pt_monte)
            outLayer2.CreateFeature(feature2)

            # creo l'ultima sezione
            # ---------------------
            NewGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[NumSez])
            feature3 = ogr.Feature(featureDefn3)
            feature3.SetField('id', NumSez)
            feature3.SetField('progr', ProgrMonte)

            feature3.SetField('zmin',Z_monte)

            pt_sinistra=NewGeom.GetPoint(0)
            Geom_pt_sinistra=ogr.Geometry(ogr.wkbPoint)
            Geom_pt_sinistra.AddPoint(pt_sinistra[0],pt_sinistra[1])
            distanza=Geom_pt_sinistra.Distance(pt_monte)
            feature3.SetField('dist1', distanza)

            pt1=pt_monte.GetPoint()
            pt_poly= pt_monte.Buffer(bufferClip)

            # taglio lungo l'ultimo tratto
            Poly_curr=ogr.CreateGeometryFromWkt(DicPoligonoClip[NumSez])
            Poly_curr=Poly_curr.Buffer(0.5)
            NewGeom=NewGeom.Intersection(Poly_curr)

            lengh2=NewGeom.Length()

            if lengh2>0:
                if NewGeom.GetGeometryCount()>0:
                    for i in range(0, NewGeom.GetGeometryCount()):
                        line=NewGeom.GetGeometryRef(i)
                        if line.Intersect(pt_poly):
                            feature3.SetGeometry(line)
                            break
                else:
                    feature3.SetGeometry(NewGeom)
            else:
                feature3.SetGeometry(NewGeom)

            outLayer3.CreateFeature(feature3)


            # salva il buffer
            dic_SezMediane_Clip[NumSez]= pt_poly.ExportToWkt()


        # =========================
        # chiudo la connessione
        outDS2.Destroy()
        outDS3.Destroy()

        # Controllo intersezione sezioni mediane
        # --------------------------------------
##        bufferClip=200.0

        #sorgente dei dati in scrittura
        ds = driver.Open(CrossMedie, 1)
        if ds is None:
            errMsg='Could not open file %s' % CrossMedie
            NotErr= bool()
            return NotErr, errMsg


        # leggo il layer dalla sorgente dei dati
        layer = ds.GetLayer()

        # conto il numero di Feature
        n = layer.GetFeatureCount()

        # dizionario sezioni da clippare al secondo ciclo while
        dic_Sez_Clip={}
        # Lista sezioni clippate al primo ciclo  while
        lista_clipped=[]

        # leggo la prima feature
        feat = layer.GetNextFeature()

        line_monte=feat.GetGeometryRef()
        line_monte_wkt=line_monte.ExportToWkt()

        NumSez_Monte=feat.GetField('id')

        feat = layer.GetNextFeature()

        while feat:

            NumSez_Valle=feat.GetField('id')

            line_valle=feat.GetGeometryRef()

            wkt=line_valle.ExportToWkt()

            # controllo se le due linee non si intersecano
            # ---------------------------------------------
            line_monte_geom=ogr.CreateGeometryFromWkt(line_monte_wkt)

            if line_valle.Crosses(line_monte_geom):

                # trova il punto di intersezione
                pt_interes_sez=line_monte_geom.Intersection(line_valle)
                # crea un buffer intorno al punto
                pt_poly = pt_interes_sez.Buffer(bufferClip)

                if  NumSez_Monte not in lista_clipped:
                    # salva il buffer
                    dic_Sez_Clip[NumSez_Monte]=pt_poly.ExportToWkt()

                # memorizza che la sezione  e' stata oggetto di clip
                lista_clipped.append(NumSez_Valle)
                # memorizza come prossima linea monte
                line_monte_wkt= line_valle.ExportToWkt()

                # esegue il clip di solo quella di valle
                # la linea viene accorciata in modo da non toccare la sezione di monte
                lines=line_valle.Difference(pt_poly)

                # conta il numero di geometrie: vi possono essere casi in cui eliminando
                # dalla linea la parte intersecata dal buffer si abbia come risultato
                # due linee di cui una  a destra ed una s sinistra
                nn=lines.GetGeometryCount()
                if nn>1:
                    poly_centro= ogr.CreateGeometryFromWkt(dic_SezMediane_Clip[NumSez_Valle])
                    for i in range(nn):
                        g = lines.GetGeometryRef(i)
                        # sceglie la parte di linea che interseca l'asse del fiume
                        if g.Crosses(poly_centro):
                            line=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            break
                else:
                    line=lines

                # salva le modifiche nello shapefile
                feat.SetGeometry(line)
                layer.SetFeature(feat)

            line_monte_wkt=wkt
            NumSez_Monte= NumSez_Valle

            # leggo la feature successiva
            feat = layer.GetNextFeature()

        # clip delle sezioni di monte
        # ---------------------------
        # riparte dall'inizio
        layer.ResetReading()

        feat = layer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')

            if NumSez in dic_Sez_Clip:

                line_valle=feat.GetGeometryRef()

                pt_poly=ogr.CreateGeometryFromWkt(dic_Sez_Clip[NumSez])

                # esegue il clip della sezione corrente
                lines=line_valle.Difference(pt_poly)

                # conta il numero di geometrie
                nn=lines.GetGeometryCount()
                if nn>1:
                    poly_centro= ogr.CreateGeometryFromWkt(dic_SezMediane_Clip[NumSez])
                    for i in range(nn):
                        g = lines.GetGeometryRef(i)
                        if g.Crosses(poly_centro):
                            line=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            break
                else:
                    line=lines

                feat.SetGeometry(line)

                layer.SetFeature(feat)

            # leggo la feature successiva
            feat = layer.GetNextFeature()


        # chiudo la connessione
        ds.Destroy()

        # fine controllo


    CreaMediana=1

    if CreaMediana>0:

        # crea nuovo shapefile
        shpnew4=LineaMediana
        if os.path.exists(shpnew4):
            driver.DeleteDataSource(shpnew4)

        outDS4 = driver.CreateDataSource(shpnew4)
        if outDS4 is None:

            errMsg='Could not create file %s' % shpnew4
            NotErr= bool()
            return NotErr, errMsg

        outLayer4 = outDS4.CreateLayer('LineaMediana', dest_srs,geom_type=ogr.wkbLineString)

        # crea il campo nello shapefile di output
        fieldDefn4 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer4.CreateField(fieldDefn4)

        featureDefne4 = outLayer4.GetLayerDefn()

        # creo la feature
        feature = ogr.Feature(featureDefne4)
        feature.SetField('id', 1)

        # inizializzo la sua geometria
        line4=ogr.Geometry(ogr.wkbLineString)

        # conto il numero di Feature in input
        n = len(ListaTratti)

        NumMonte=0

        for NumPoint in NumSezTipo:

            pt_curr=ogr.CreateGeometryFromWkt(NumSezGeom[NumPoint])

            punto=pt_curr.GetPoint()

            # aggiunge al dizionario
            dic_PointPath[NumPoint]=punto

            # salta le sezioni non salvate
            if NumPoint in ListaNumSezTotOk:
                line4.AddPoint(punto[0],punto[1])

        feature.SetGeometry(line4)
        outLayer4.CreateFeature(feature)

        # salva in formato testo
        lineaMedianaWkt=line4.ExportToWkt()
        # azzera la memoria
        line4.Destroy()

        outDS4.Destroy()


    # Crea : Cross_Tot_Medie
    # ---------------------------
    CreaCross_Tot_Medie=1

    if CreaCross_Tot_Medie>0:

        #sorgente dei dati in lettura
        ds = driver.Open(CrossSecTot, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossSecTot
            NotErr= bool()
            return NotErr, errMsg

        # leggo il layer dalla sorgente dei dati
        Inlayer = ds.GetLayer()
        Spatialref = Inlayer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        SourceEPSG=int(Spatialref.GetAuthorityCode(None))

        # crea nuovo shapefile delle sezioni
        # ----------------------------------
        shpnew3=Cross_Tot_Medie
        if os.path.exists(shpnew3):
            driver.DeleteDataSource(shpnew3)

        outDS3 = driver.CreateDataSource(shpnew3)
        if outDS3 is None:
            errMsg='Could not create file %s' % shpnew3
            NotErr= bool()
            return NotErr, errMsg

        outLayer3 = outDS3.CreateLayer('SezioniTotMedie', dest_srs,geom_type=ogr.wkbLineString)

        # crea i nuovi campi id e Nome nello shapefile di output
        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer3.CreateField(fieldDefn2)

        # crea il campo zmin
        fieldDefnZmin = ogr.FieldDefn('zmin', ogr.OFTReal)
        outLayer3.CreateField(fieldDefnZmin)

        featureDefn3 = outLayer3.GetLayerDefn()

        # leggo la prima feature
        feat = Inlayer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')
            z_curr=float(feat.GetField('zmin'))

            feature = ogr.Feature(featureDefn3)
            feature.SetField('id', NumSez)
            feature.SetField('zmin', z_curr)

            # leggo la geometria
            linea=feat.GetGeometryRef()
            feature.SetGeometry(linea)
            outLayer3.CreateFeature(feature)


            feat = Inlayer.GetNextFeature()

         # chiudo la connessione
        ds.Destroy()

        # aggiungo le sezioni intermedie
         #sorgente dei dati in lettura
        ds = driver.Open(CrossMedie, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossMedie
            NotErr= bool()
            return NotErr, errMsg

        # leggo il layer dalla sorgente dei dati
        Inlayer = ds.GetLayer()

        # leggo la prima feature e la salto
        feat = Inlayer.GetNextFeature()

        # leggo la seconda feature
        feat = Inlayer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')
            z_curr=float(feat.GetField('zmin'))

            feature = ogr.Feature(featureDefn3)
            feature.SetField('id', -NumSez)
            feature.SetField('zmin', z_curr)

            # leggo la geometria
            linea=feat.GetGeometryRef()
            feature.SetGeometry(linea)
            outLayer3.CreateFeature(feature)

            feat = Inlayer.GetNextFeature()


         # chiudo la connessione
        ds.Destroy()
        # salvo
        outDS3.Destroy()



    # crea i poligoni dei tratti  di fiume divisi in destra e sinistra
    # ----------------------------------------------------------------
    CreaPolygoni2=1

    if CreaPolygoni2>0:

        #sorgente dei dati in lettura
        ds = driver.Open(CrossSecPoly, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossSecPoly
            NotErr= bool()
            return NotErr, errMsg


        # leggo il layer dalla sorgente dei dati
        Inlayer = ds.GetLayer()
        Spatialref = Inlayer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        SourceEPSG=int(Spatialref.GetAuthorityCode(None))


        # crea nuovo shapefile
        shpnew5=CrossSecPoly2
        if os.path.exists(shpnew5):
            driver.DeleteDataSource(shpnew5)

        outDS5 = driver.CreateDataSource(shpnew5)
        if outDS5 is None:
            errMsg='Could not create file %s' % shpnew5
            NotErr= bool()
            return NotErr, errMsg

        outLayer5 = outDS5.CreateLayer('CrossSecPoly_2', Spatialref,geom_type=ogr.wkbPolygon)

        # crea i nuovi campi id e Nome nello shapefile di output
        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn2)
        fieldDefn5 = ogr.FieldDefn('lato', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn5)

        featureDefn5 = outLayer5.GetLayerDefn()

        feat = Inlayer.GetNextFeature()

        while feat:

            # legge il numero
            NumSez=feat.GetField('id')

            # legge la geomeria
            PuntiPoligono=DicPuntiPoligono[NumSez]
            PuntiTratto=DicPuntiAsse[NumSez]

            # legge il poligono del tratto clippato
            Poly_curr=ogr.CreateGeometryFromWkt(DicPoligonoClip[NumSez])

            # crea i due poligoni
            for lato in range(2):

                feature = ogr.Feature(featureDefn5)
                feature.SetField('id', NumSez)
                feature.SetField('lato', lato)

                # crea il poligono
                ring = ogr.Geometry(ogr.wkbLinearRing)

                if lato==0:

                    # punto a monte sinistro
                    pt=PuntiPoligono[0]
                    ring.AddPoint(pt[0],pt[1])

                    # si aggiungono i punti dell'asse del fiume
                    # ..........................................
                    for pt in PuntiTratto:
                        ring.AddPoint(pt[0],pt[1])

                    # punto a valle sinistro
                    pt=PuntiPoligono[3]
                    ring.AddPoint(pt[0],pt[1])

                    # punto a monte sinistro
                    pt=PuntiPoligono[0]
                    ring.AddPoint(pt[0],pt[1])

                    ring.CloseRings()
                    poly_new=ogr.Geometry(ogr.wkbPolygon)
                    poly_new.AddGeometry(ring)

                else:

                    # punto a monte centrale
                    # uso il primo punto del tratto
                    pt= PuntiTratto[0]
                    ring.AddPoint(pt[0],pt[1])

                    # punto a monte destro
                    pt=PuntiPoligono[1]
                    ring.AddPoint(pt[0],pt[1])

                    # punto a valle sinistro
                    pt=PuntiPoligono[2]
                    ring.AddPoint(pt[0],pt[1])

                    # aggiunta punti dell'asse del fiume da valle verso monte
                    for iii in range(len(PuntiTratto)-1,-1,-1):
                        pt=PuntiTratto[iii]
                        ring.AddPoint(pt[0],pt[1])

                    ring.CloseRings()

                    poly_new=ogr.Geometry(ogr.wkbPolygon)
                    poly_new.AddGeometry(ring)

                # interseca il poligono del tratto corrente
                NewPoly=poly_new.Intersection(Poly_curr)

                if NewPoly!=None:
                    npoly=NewPoly.GetGeometryCount()
                    if  npoly>0:
                        numpt=NewPoly.GetGeometryRef(0).GetPointCount()
                    else:
                        # geometria assente !!
                        poly_wkt=NewPoly.ExportToWkt()
                        # uso quella prima della intersezione
                        NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                        npoly=poly_new.GetGeometryCount()
                        geom_poly=poly_new.GetGeometryRef(0)
                        numpt=geom_poly.GetPointCount()

                        poly_wkt=NewPoly.ExportToWkt()
                        # salvo i dati
                        nome_prova=PathFiles +os.sep+ 'File_poly_%d.csv' % NumSezValle
                        fout=open(nome_prova,'w')
                        txt='X;Y\n'
                        fout.write(txt)
                        for iii in range(0, numpt):
                            # GetPoint returns a tuple not a Geometry
                            pt = geom_poly.GetPoint(iii)
                            txt='%s;%s\n' %(pt[0],pt[1])
                            fout.write(txt)
                        fout.close()
                else:

                    # casi di mancata intersezione ??????
                    NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                    npoly=poly_new.GetGeometryCount()
                    geom_poly=poly_new.GetGeometryRef(0)
                    numpt=geom_poly.GetPointCount()

                    poly_wkt=NewPoly.ExportToWkt()
                    # salvo i dati
                    nome_prova=PathFiles +os.sep+ 'File_poly_%d.csv' % NumSezValle
                    fout=open(nome_prova,'w')
                    txt='X;Y\n'
                    fout.write(txt)
                    for iii in range(0, numpt):
                        # GetPoint returns a tuple not a Geometry
                        pt = geom_poly.GetPoint(iii)
                        txt='%s;%s\n' %(pt[0],pt[1])
                        fout.write(txt)
                    fout.close()

                # =================================
                feature.SetGeometry(NewPoly)
                outLayer5.CreateFeature(feature)

                ring.Destroy()
                poly_new.Destroy()
                NewPoly.Destroy()

            feat = Inlayer.GetNextFeature()

        # chiudo la connessione
        outDS5.Destroy()
        ds.Destroy()

    # creazione poligoni sinistra e destra fiume
    # -------------------------------------------
    CreaPolygoni3=1

    if CreaPolygoni3>0:

        # Creo il poligono unione dei tratti da utilizzare per il clip del fiume
        # ----------------------------------------------------------------------
        PoligonoUnioneTratti_0=ogr.CreateGeometryFromWkt(PoligonoUnioneTratti_Wkt)


        # crea la lista dei punti della streamLine
        # ========================================
        End_cross_sec_line=ogr.CreateGeometryFromWkt(wkt_sez_end)
        # interseca il poligono unione dei tratti
        NewStreamLine=StreamLine.Intersection(PoligonoUnioneTratti_0)
        if  NewStreamLine!=None:
            ListaPuntiStreamLine= SetStreamLinePoints(NewStreamLine,End_cross_sec_line)
        else:
            ListaPuntiStreamLine= SetStreamLinePoints(StreamLine,End_cross_sec_line)

        # crea nuovo shapefile
        shpnew6=PolySxDx
        if os.path.exists(shpnew6):
            driver.DeleteDataSource(shpnew6)

        outDS5 = driver.CreateDataSource(shpnew6)
        if outDS5 is None:
            errMsg='Could not open file %s' % shpnew6
            NotErr= bool()
            return NotErr, errMsg

        outLayer5 = outDS5.CreateLayer('PolySxDx', Spatialref,geom_type=ogr.wkbPolygon)

        # crea i nuovi campi id e Nome nello shapefile di output
        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn2)
        fieldDefn5 = ogr.FieldDefn('lato', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn5)

        featureDefn5 = outLayer5.GetLayerDefn()


        # creo il poligono sinistro
        # ---------------------------
        feature = ogr.Feature(featureDefn5)
        feature.SetField('id', 0)
        feature.SetField('lato', 0)

        ring = ogr.Geometry(ogr.wkbLinearRing)

        numPunti=len(NumPolyList)

        # lato sinistro
        for ii in range(numPunti):
            PuntiPoligono=DicPuntiPoligono[NumPolyList[ii]]
            # punto a monte sinistro
            pt=PuntiPoligono[0]
            ring.AddPoint(pt[0],pt[1])

        # aggiungo il punto sinistro dell'ultima sezione indice 3
        PuntiPoligono=DicPuntiPoligono[NumPolyList[-1]]
        pt=PuntiPoligono[3]
        ring.AddPoint(pt[0],pt[1])

        # aggiungo la line del fiume da valle verso monte
        for iii in range(len(ListaPuntiStreamLine)-1,-1,-1):
            pt=ListaPuntiStreamLine[iii]
            ring.AddPoint(pt[0],pt[1])

        ring.CloseRings()

        poly_sx=ogr.Geometry(ogr.wkbPolygon)
        poly_sx.AddGeometry(ring)

        # interseca il poligono del dominio
        NewPoly=poly_sx.Intersection(DominoRing)

        if  NewPoly!=None:
            feature.SetGeometry(NewPoly)
            outLayer5.CreateFeature(feature)

            ring.Destroy()
            poly_sx.Destroy()
            NewPoly.Destroy()
        else:
            feature.SetGeometry(poly_sx)
            outLayer5.CreateFeature(feature)

            ring.Destroy()
            poly_sx.Destroy()

        # crea il poligono destro
        # -------------------------
        feature = ogr.Feature(featureDefn5)
        feature.SetField('id', 1)
        feature.SetField('lato', 1)

        ring = ogr.Geometry(ogr.wkbLinearRing)

        # aggiungo la line del fiume da monte verso valle
        for iii in range(len(ListaPuntiStreamLine)):
            pt=ListaPuntiStreamLine[iii]
            ring.AddPoint(pt[0],pt[1])

        # lato destro da valle verso monte : indice 2
        for ii in range(numPunti-1,-1,-1):
            PuntiPoligono=DicPuntiPoligono[NumPolyList[ii]]
            pt=PuntiPoligono[2]
            ring.AddPoint(pt[0],pt[1])

        # aggiunge l'ultimo punto
        PuntiPoligono=DicPuntiPoligono[NumPolyList[0]]
        pt=PuntiPoligono[1]
        ring.AddPoint(pt[0],pt[1])

        ring.CloseRings()

        poly_dx=ogr.Geometry(ogr.wkbPolygon)
        poly_dx.AddGeometry(ring)

        # interseca il poligono del dominio
        NewPoly=poly_dx.Intersection(DominoRing)

        if  NewPoly!=None:
            feature.SetGeometry(NewPoly)
            outLayer5.CreateFeature(feature)

            ring.Destroy()
            poly_dx.Destroy()
            NewPoly.Destroy()

        else:
            feature.SetGeometry(poly_dx)
            outLayer5.CreateFeature(feature)

            ring.Destroy()
            poly_dx.Destroy()

         # chiudo la connessione
        outDS5.Destroy()


    # Creazione StreamDH.tif
    # ------------------------
    CreaCrid2=1

    if CreaCrid2>0:

        if not os.path.exists(ClipDEM):
            errMsg='File ClipDEM %s does not exists' % os.path.realpath(ClipDEM)
            NotErr= bool()
            return NotErr, errMsg

        infile=ClipDEM

        indatasetElev = gdal.Open( infile, GA_ReadOnly )
        if indatasetElev is None:
            errMsg='Could not open ' + infile
            NotErr= bool()
            return NotErr, errMsg

        prj = indatasetElev.GetProjectionRef()

        geotransform = indatasetElev.GetGeoTransform()

        originXElev = geotransform[0]
        originYElev = geotransform[3]
        pixelWidthElev = geotransform[1]
        pixelHeightElev = geotransform[5]
        colsElev=indatasetElev.RasterXSize
        rowsElev=indatasetElev.RasterYSize
        bandsElev=indatasetElev.RasterCount
        iBand = 1
        inbandElev = indatasetElev.GetRasterBand(iBand)
        inNoDataElev= inbandElev.GetNoDataValue()
        # reading the entire file at once
        terreno = inbandElev.ReadAsArray(0, 0, colsElev, rowsElev).astype(np.float32)

        # apre file curve con le sezioni
        # ------------------------------
        inDS1 = driver.Open(Cross_Tot_Medie, 0)
        if inDS1 is None:
            errMsg='Could not open ' + CrossSecTot
            NotErr= bool()
            return NotErr, errMsg

        # apro il layer in lettura
        InlayerCurve = inDS1.GetLayer()

        # mi leggo il sistema di riferimento
        spatialRef_sez=InlayerCurve.GetSpatialRef()

        feat_defn = InlayerCurve.GetLayerDefn()
        NumFields=feat_defn.GetFieldCount()
        flag_WL=0
        nomecampoLivello='zmin'
        for i in range(NumFields):
            nome = feat_defn.GetFieldDefn(i).GetName()
            if nomecampoLivello==nome:
                flag_WL=1
                break

        if flag_WL>0:

            # genera un grid con la quota minima delle sezioni
            GridSez=np.zeros((rowsElev,colsElev),np.float32)

            format = 'MEM'
            type = GDT_Float32

            driver2 = gdal.GetDriverByName(format)
            driver2.Register()

            gt=indatasetElev.GetGeoTransform()

            ds = driver2.Create('GridSez', indatasetElev.RasterXSize, indatasetElev.RasterYSize, 1, type)
            if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
                ds.SetGeoTransform(gt)

            # imposta il sistema di riferimento uguale al modello del tiranti: se manca imposta il default
            if prj is not None and len(prj) > 0:
                ds.SetProjection(prj)
            else:
                prj= spatialRef.ExportToWkt()
                ds.SetProjection(prj)

            iBand=1
            #CampoValore=["ATTRIBUTE=OBJECTID"]
            testo="ATTRIBUTE=%s"  % (nomecampoLivello)
            # Rasterize
            outband = ds.GetRasterBand(iBand)

            # Rasterize
            # azzero
            outband.WriteArray(GridSez, 0, 0)
            CampoValore=[testo]

            # creo la mappa dei valori
            # -------------------------------------------------
            err = gdal.RasterizeLayer(ds, [iBand], InlayerCurve,
                    burn_values=[0],
                    options=CampoValore)
            if err != 0:
                raise Exception("error rasterizing layer: %s" % err)

            # Reading WL
            GridSezWL = outband.ReadAsArray().astype(np.float32)

            ds= None

            # INTERPOLATE Water Level Grid
            # ----------------------------

            #size of grid
            xmin=originXElev
            xmax=xmin+colsElev*pixelWidthElev
            ymax=originYElev
            ymin=originYElev+rowsElev*pixelHeightElev

            nx = int((xmax - xmin + 1)/pixelWidthElev)
            ny = int(-(ymax - ymin + 1)/pixelHeightElev)

            # Generate a regular grid to interpolate the data.
            xi = np.linspace(xmin, xmax, nx)
            yi = np.linspace(ymin, ymax, ny)
            xi, yi = np.meshgrid(xi, yi)

            # Reading x,y,z
            mask=GridSezWL>0
            x=xi[mask]
            y=yi[mask]
            z=GridSezWL[mask]


            # Interpolate the values of z for all points in the rectangular grid

            # Interpolate  using scipy interpolate griddata
            WLArray = il.griddata((x, y), z, (xi, yi),method='linear') #(may use 'nearest', 'linear' or 'cubic'  - although constant problems w linear)

            maskTerreno=terreno<-10.0
            Wdepth=terreno-WLArray
            Nodata=-9999
            Wdepth= np.choose(maskTerreno,(Wdepth,Nodata))

            checkMask=numpy.isnan(Wdepth)

            nnan=checkMask.sum()

            if nnan>0:
                Wdepth[checkMask]=Nodata
            # nome del file output
            PathFiles=os.path.dirname(ClipDEM)
            FileDEM_out=PathFiles+os.sep+'StreamDH.tif'

            format = 'GTiff'
            driver = gdal.GetDriverByName(format)
            type = GDT_Float32
            gt=indatasetElev.GetGeoTransform()

            ds = driver.Create(FileDEM_out, indatasetElev.RasterXSize, indatasetElev.RasterYSize, 1, type)
            if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
                ds.SetGeoTransform(gt)

            # sets the reference system equal to the depth map of water: if it lacks sets the default
            if prj is not None and len(prj) > 0:
                ds.SetProjection(prj)
            else:
                prj= spatialRef.ExportToWkt()
                ds.SetProjection(prj)

            # writing raster
            iBand=1
            outband = ds.GetRasterBand(iBand)
            outband.WriteArray(Wdepth, 0, 0)

            outband.FlushCache()
            outband.SetNoDataValue(Nodata)
            outband.GetStatistics(0,1)
            outband = None

            ds = None

        inDS1.Destroy()


    CreaCrid3=1

    if CreaCrid3>0:

        indataset = gdal.Open(ClipDEM, GA_ReadOnly )
        if indataset is None:
            errMsg='Could not open ' + ClipDEM
            NotErr= bool()
            return NotErr, errMsg

        geotransform = indataset.GetGeoTransform()

        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        cols=indataset.RasterXSize
        rows=indataset.RasterYSize
        bands=indataset.RasterCount
        iBand = 1
        inband = indataset.GetRasterBand(iBand)
        inNoData= inband.GetNoDataValue()
        prj = indataset.GetProjectionRef()

        spatialRef = osr.SpatialReference()
        try:
            spatialRef.ImportFromWkt(prj)
        except:
            pass
            # definisco un default per il riferimento
            # imposto WGS84 UTM 32 N
            spatialRef.ImportFromEPSG(32632)

        terreno = inband.ReadAsArray(0, 0, cols, rows).astype(numpy.float)
        mask_Nodata= terreno==inNoData


        # rasterizza i poligoni
        # ====================
        orig_data_source = ogr.Open(CrossSecPoly2)
        source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_source, "")
        source_layer = source_ds.GetLayer()


        format = 'Gtiff'
        type = GDT_Int16

        driver3 = gdal.GetDriverByName(format)
        driver3.Register()

        PathFiles=os.path.dirname(ClipDEM)
        TrattiRaster=PathFiles+os.sep+'DestraSinistra.tif'

        dsRaster = driver3.Create(TrattiRaster, cols, rows, 1, type)
        gt1=geotransform
        if gt1 is not None and gt1 != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            dsRaster.SetGeoTransform(gt1)

        # imposta il sistema di riferimento uguale al modello del Terreno: se manca imposta il default
        if prj is not None and len(prj) > 0:
            dsRaster.SetProjection(prj)
        else:
            prj= spatialRef.ExportToWkt()
            dsRaster.SetProjection(prj)

        # Rasterize
        iBand=1
        outband = dsRaster.GetRasterBand(iBand)

        outNodata=-9999

        # scrive -1 su tutta la matrice
        ClassTratti=numpy.zeros((rows,cols)).astype(numpy.int)
        ClassTratti=ClassTratti-1

        outband.WriteArray(ClassTratti, 0, 0)

        # Rasterize
        err = gdal.RasterizeLayer(dsRaster, [1], source_layer,
                burn_values=[0],
                options=["ATTRIBUTE=lato"])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)

        # scrive i Nodata
        MatriceDati=outband.ReadAsArray(0, 0, cols, rows)
        MatriceDati=numpy.choose(mask_Nodata,(MatriceDati,outNodata))
        outband.WriteArray(MatriceDati, 0, 0)

        outband.FlushCache()
        outband.SetNoDataValue(outNodata)
        outband.GetStatistics(0,1)

        outband=None

        dsRaster=None
        # chiude le sorgenti dei dati
        orig_data_source.Destroy()



    return NotErr, errMsg

if __name__ == '__main__':


    mydb_path_user='..'+ os.sep+'db'+os.sep+'USER_GeoDB.sqlite'

    # San Giuliano
    ID_Diga=449
    PathFiles='..'+ os.sep+ str(ID_Diga)
    fileDEM=  PathFiles+ os.sep+'DTM_clip.tif'

##    NotErr, errMsg= SetIntermediatePoints(mydb_path_user,ID_Diga,fileDEM)
##
##    print(NotErr,errMsg)

    NotErr, errMsg= SetCrossSec_2(mydb_path_user,PathFiles,fileDEM,ID_Diga)

    print(NotErr,errMsg)
