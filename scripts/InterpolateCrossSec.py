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

Thee script reads:
    - the shapefile of the intersection points of the cross sections with the river path
      and the points on the axis where the creation of the interpolated sections is required
    - the shapefile of the main cross-sections from which to interpolate
    - the shapefile of the domain outline
    - the grid of DTM_clip

creates:

    - the shapefile of the total cross sections (original and interpolated) clipped
      along the outline of the domain

    - a polygonal shapefile composed of trapezoids between a cross section and the next one

    - a grid with the classification of the cells according to the id of the polygons that
      corresponds to the number of cells along the river path distant from the start point

    - a grid Tratti.tif with the river reach classes starting from the dam
      the reach number represents in number of cells of distance counted along the river path

    - a polygonal shapefile composed of trapezoids between a cross section and the next one
      dividing the trapezoid on the left (side = 0) and right (side = 1)

    - a grid DestraSinistra.tif with an indication of whether the cell is located at
      left (= 0) or right (= 1) of the river

    - a grid StreamDH with in values of the heights of the ground with respect to the river
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
    from osgeo.osr import osr
except:
    import ogr
    import osr
import numpy as np
import math
import matplotlib.pyplot as plt
import time
import scipy.interpolate as il

def isLeft(a,b,c):
    # check if the point c (x, y) is to the left of the segment a-b
    position=((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))

    return 1*numpy.sign(position)


def CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select):

    # control the line direction

    if inters_destra:
        # distance control to the right
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(line_select.GetX(1),line_select.GetY(1))
        dist_dx=pt_interes_sez.Distance(point)
        point.Destroy()
        if dist_dx > distBufferClip:
            # to flip
            line= flip_line(line_select)
        else:
            line=line_select
    else:
        # distance control to the left
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
    graph debug of a polygon
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
    graph debug two polygons
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
    Calculate the cross product between
    two vectors a and b in the x, y plane

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
    Creates the trace of a section having a direction defined by two points
    and centered on a point
    p0 [x,y]    : central point
    Length      : half-length
    psx [x,y]   : left point of the direction
    pdx [x,y]   : right point of direction
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

def SezioneInterpolata(pt_curr,SezMonteGeom,SezValleGeom):
    """
    Find an intermediate section between
    - SezMonteGeom: geom of the segment of the upstream section
    - SezValleGeom: geom of the segment of the downstream section
    - pt_curr     : intermediate point between the two sections and from which it
                    passes the interpolated section
                    the second point of the interpolated section is
                    the one in which the straight lines
                    of the two sections of upstream and downstream
                    they intersect each other
    """

    # Half length upstream section
    Length=SezMonteGeom.Length()/2.0*1.5

    # calculates the interpolated geometry
    # ------------------------------------
    p1 = np.array( [SezMonteGeom.GetX(0), SezMonteGeom.GetY(0)] )
    p2 = np.array( [SezMonteGeom.GetX(1), SezMonteGeom.GetY(1)] )


    p3 = np.array( [SezValleGeom.GetX(0), SezValleGeom.GetY(0)] )
    p4 = np.array( [SezValleGeom.GetX(1), SezValleGeom.GetY(1)] )

    # intersection of the two sections of upstream and downstream
    # -------------------------------------------------
    #this is the point at which the interpolated section must converge
    # and it can be in the right or left depending on the curvature
    # of the axis of the river
    Intersect_0=seg_intersect( p1,p2, p3,p4)

    # check whether to the right or left of the river
    orientamento=ProdottoVettSezioni(SezMonteGeom,SezValleGeom)

    # control on the joining line of the direct sections to the inner part of the curve
    LineInterna=ogr.Geometry(ogr.wkbLineString)

    if  orientamento>0:
        # an intersection is assumed on the left bank
        inters_sx=bool('True')
        inters_destra=bool()
        LineInterna.AddPoint(SezMonteGeom.GetX(0), SezMonteGeom.GetY(0))
        LineInterna.AddPoint(SezValleGeom.GetX(0), SezValleGeom.GetY(0))

    else:
        inters_sx=bool()
        inters_destra=bool('True')
        LineInterna.AddPoint(SezMonteGeom.GetX(1), SezMonteGeom.GetY(1))
        LineInterna.AddPoint(SezValleGeom.GetX(1), SezValleGeom.GetY(1))


     # extends the length of the section
    # -----------------------------------
    if  inters_sx:
        GeomNew=SetCrossSecPointDiretion(pt_curr.GetPoint(0),Length,Intersect_0,pt_curr.GetPoint(0))
    else:
        GeomNew=SetCrossSecPointDiretion(pt_curr.GetPoint(0),Length,pt_curr.GetPoint(0),Intersect_0)



    # check if the section reaches the intersection point of the main ones
    # --------------------------------------------------------------------------------
    bufferClip=50.0
    distBufferClip=bufferClip*1.1

    # creates a buffer around the point
    pt_interes_sez=ogr.Geometry(ogr.wkbPoint)
    pt_interes_sez.AddPoint(Intersect_0[0],Intersect_0[1])
    pt_poly = pt_interes_sez.Buffer(bufferClip)


    if GeomNew.Crosses(pt_poly):

        pt_curr_poly= pt_curr.Buffer(bufferClip)

        # the line is shortened
        lines=GeomNew.Difference(pt_poly)

        # counts the number of geometries
        nn=lines.GetGeometryCount()
        if nn>1:
            for i in range(nn):
                g = lines.GetGeometryRef(i)
                # choose the part of the line that intersects the axis of the river
                if g.Crosses(pt_curr_poly):
                    line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                    # check the direction
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                    break
        else:
            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

        wkt_new=line.ExportToWkt()
        GeomNew=ogr.CreateGeometryFromWkt(wkt_new)

    # check the intersection with the union line of the main sections
    # ----------------------------------------------------------------------
    line_poly=LineInterna.Buffer(1.0)

    if GeomNew.Crosses(LineInterna):

        pt_curr_poly= pt_curr.Buffer(bufferClip)

        # find the intersection point
        pt_interes_sez=GeomNew.Intersection(LineInterna)


        # the line is shortened
        lines=GeomNew.Difference(line_poly)

        # counts the number of geometries
        nn=lines.GetGeometryCount()
        if nn>1:
            for i in range(nn):
                g = lines.GetGeometryRef(i)
                # choose the part of the line that intersects the axis of the river
                if g.Crosses(pt_curr_poly):
                    line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                    # check the direction
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                    break
        else:
            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

        wkt_new=line.ExportToWkt()
        GeomNew=ogr.CreateGeometryFromWkt(wkt_new)

    return GeomNew


def PuntoMedio(pt1,pt2):

    # calculates the midpoint between two
    pt_mean=[(pt1[0]+pt2[0])/2.0,(pt1[1]+pt2[1])/2.0]

    return pt_mean

def PuntoIntermedio(pt1,pt2,percent):

    # calculates the intermediate point at a certain percentage between two
    pt_intermedio=[pt1[0]+percent*(pt2[0]-pt1[0]),pt1[1]+percent*(pt2[1]-pt1[1])]

    return pt_intermedio

def flip_line(line):
    # reverses the points of the line
    GeomNew=ogr.Geometry(ogr.wkbLineString)
    for i in range(line.GetPointCount()-1,-1,-1):
        # GetPoint returns a tuple not a Geometry
        pt = line.GetPoint(i)
        GeomNew.AddPoint(pt[0], pt[1])

    return GeomNew


def SetStreamLinePoints(MultiStreamLine,End_cross_sec_line):
    """
    Create the list of river points from the first to the last section
    Input:
        -  StreamLine           : line of the river path
        -  End_cross_sec_line   : final cross-section line
    """
    # create the list of points in the streamLine
    ListaPuntiStreamLine=[]
    if MultiStreamLine.GetGeometryName() == 'LINESTRING':
        StreamLine=MultiStreamLine
    else:

        # check the number of geometries
        # takes the largest polygon in the Multipoligons
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
            # takes the first one
            StreamLine=MultiStreamLine.GetGeometryRef(0)
            Length=StreamLine.Length()

    # entering the first point
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

def SetIntermediatePoints(mydb_path_user,PathFiles,fileDEM,DamID,DistanzaSezInterp=2000.0,DeltaSezPrincipale=5):
    """
    Create the sequence of intermediate points between the main cross sections
    """
    NotErr=bool('True')
    errMsg='OK'

    if not os.path.exists(PathFiles):
        errMsg = "There is no data for the dam num =%s \nPerform the editing of the main cross sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    if not os.path.exists(fileDEM):
        errMsg = "There is no data for the dam num =%s \nMake the downstream DEM clip first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
   # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    # check downstream line existence
    NomeTabellaLinee='Downstreampath'

    # reference system code of the table
    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
    cur.execute(sql)

    record=cur.fetchone()
    if record!=None:
        DominioEPSG=record[0]
    else:
        DominioEPSG=32632

    dest_srs = ogr.osr.SpatialReference()
    dest_srs.ImportFromEPSG(DominioEPSG)


    sql='SELECT ST_AsText(geom) FROM %s WHERE DamID=%d' % (NomeTabellaLinee,DamID)
    cur.execute(sql)
    ChkDiga=cur.fetchone()

    if ChkDiga==None:
        errMsg = "In the table= %s there is no data for the dam num =%s \nPerform the downstream line calculation first !" % (NomeTabellaLinee,DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt=ChkDiga[0]
        linea_geom=ogr.CreateGeometryFromWkt(wkt)
        TotalLength=linea_geom.Length()


    # table of main cross-sections
    NomeTabella='MainCrossSec'

    sql='SELECT PKUID,ST_AsText(geom) FROM %s WHERE DamID=%d' % (NomeTabella,DamID)
    sql+=' ORDER BY id'
    sql+=';'
    cur.execute(sql)
    ListaCrossSec=cur.fetchall()

    n=len(ListaCrossSec)

    if n <= 0:

        errMsg = "There is no data for the dam num =%s \nFirst calculate the main cross-sections !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:

        # reference system code of the table
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        cur.execute(sql)
        record=cur.fetchone()
        if record!=None:
            OriginEPSG=record[0]
        else:
            OriginEPSG=32632


    # --------------
    # Reading the DTM
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

    # defining a minimum distance between secondary cross-sections
    DistanzaSezInterp_min=3*pixelWidth

    # Create the list of the points of the river line
    # ---------------------------------------------
    ListaPuntiFiume=[]
    for i in range(0, linea_geom.GetPointCount()):
        ListaPuntiFiume.append(linea_geom.GetPoint(i))

    num_pt_fiume=len(ListaPuntiFiume)

    # dictionary PKUID_sec_id: stores the id values to be assigned to the main cross-sections
    PKUID_sec_id={}
    listaIdSecPrincipali=[]

    # dictionary id_Points : dictionary of points
    id_Points={}
    # dictionary of progressive distances
    progr_Points={}


    # -----------------------
    point_stream_start=0
    num_point_curr=0
    progr_monte=0.0
    pt_start=ListaPuntiFiume[0]
    id_Points[num_point_curr]=pt_start
    progr_Points[num_point_curr]=progr_monte

    # looking for the first main cross-section
    cross_sec_wkt=ListaCrossSec[0][1]

    # set id = 0 at the point of the first section
    PKUID_sec_id[ListaCrossSec[0][0]]=num_point_curr
    listaIdSecPrincipali.append(num_point_curr)

    cross_sec_line=ogr.CreateGeometryFromWkt(cross_sec_wkt)

    # stores the line of the corresponding interpolated section
    UpstreaSez_wkt=cross_sec_line.ExportToWkt()
    # stores the points of the section
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

                # sequence of secondary progressive distances
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
                # creating the progressive array
                progr_along=xvals + progr_monte

                # finding points at predefined distances by linear interpolation
                ArrayXSec=numpy.interp(xvals, dist_along, x0)
                ArrayYSec=numpy.interp(xvals, dist_along, y0)


                # creation of intermediate points
                # -----------------------------
                npt=len(ArrayXSec)

                pt_valle=[ArrayXSec[-1],ArrayYSec[-1]]

                Tooll_Distanza=tratto_length/float(npt)/4.0

                for i in range(1,npt-1):

                    num_point_curr+=1

                    pt_curr=ogr.Geometry(ogr.wkbPoint)
                    pt_curr.AddPoint(ArrayXSec[i],ArrayYSec[i])

                    # check the position of the i-th point with respect to the upstream segment
                    c=pt_curr.GetPoint(0)
                    PosizPtCur=isLeft(Point_a,Point_b,c)

                    # check the downstream point position with respect to the downstream segment
                    PosizPtValle=isLeft(Point_a,Point_b,pt_valle)


                    # save the upstream segment data for checking the next point
                    NewGeom = SezioneInterpolata(pt_curr,cross_sec_line,next_cross_sec_line)
                    Point_a=NewGeom.GetPoint(0)
                    Point_b=NewGeom.GetPoint(1)

                    # check the current mutual position of the two points
                    Congruenza=PosizPtCur*PosizPtValle

                    if Congruenza>0:

                        # also check that there is a minimum distance from the upstream section
                        # and from the directional downstream cross-section
                        UpstreaSez=ogr.CreateGeometryFromWkt(UpstreaSez_wkt)
                        Dist1=pt_curr.Distance(UpstreaSez)
                        Dist2=pt_curr.Distance(next_cross_sec_line)

                        if Dist1>Tooll_Distanza and Dist2>Tooll_Distanza:

                            # saving the point coordinates
                            id_Points[num_point_curr]=pt_curr.GetPoint(0)
                            progr_Points[num_point_curr]=progr_along[i]

                            # saves the geometry of the section upstream of the point i + 1-th
                            UpstreaSez_wkt=NewGeom.ExportToWkt()

                        else:
                            txt='--  Discarded point n: %d of the intermediate sections downstream cross-section directional : %d-th due to distance less than the minimum' % (i,isec-1)
                            errMsg+='\n%s'% txt

                            pt_curr.Destroy()

                    else:
                        txt='--  Discarded point n: %d  of the intermediate sections downstream cross-section directional :%d-th due to a meander that goes upstream' % (i,isec-1)
                        errMsg+='\n%s'% txt
                        pt_curr.Destroy()

                # -----------------------------------------------------------
                # Creation of the point of the main downstream cross-section
                # ----------------------------------------------------------

                # increase the number of points id
                num_point_curr+=1

                # creation of the downstream point
                PKUID_sec_id[ListaCrossSec[isec][0]]=num_point_curr
                listaIdSecPrincipali.append(num_point_curr)

                # saving the point coordinates
                id_Points[num_point_curr]=pt_valle
                progr_Points[num_point_curr]=progr_along[-1]

                # updating the upstream progressive
                progr_monte+=tratto_length

                # saves the geometry of the section upstream of the point i + 1-th
                UpstreaSez_wkt=ListaCrossSec[isec][1]

                cross_sec_line=ogr.CreateGeometryFromWkt(UpstreaSez_wkt)

                # save the points of the section
                Point_a=cross_sec_line.GetPoint(0)
                Point_b=cross_sec_line.GetPoint(1)

                # initialize the new reach
                line_tratto.Destroy()
                line_tratto=ogr.Geometry(ogr.wkbLineString)
                line_tratto.AddPoint(pt_valle[0],pt_valle[1])

                # update and exit the cycle
                point_stream_start=i_fiume

                break

                # clear the memory
                line_tratto.Destroy()

            else:
                line_tratto.AddPoint(ListaPuntiFiume[i_fiume+1][0],ListaPuntiFiume[i_fiume+1][1])

    # check the case of the last section
    # -------------------------------------
    NumTratti= NumSecPrincipali-1

    if NumTratti not in ListaSecOk:

        # case in which the last section (line_tratto) has not been evaluated
        tratto_length=line_tratto.Length()

        DistanzaSezInterpCurr=tratto_length/float(DeltaSezPrincipale)

        if DistanzaSezInterp<DistanzaSezInterp_min:
            DeltaSezPrincipaleCur=int(tratto_length/DistanzaSezInterp_min)
            if DeltaSezPrincipaleCur<1:
                DeltaSezPrincipaleCur=1
        else:
            DeltaSezPrincipaleCur=DeltaSezPrincipale*1


        NumSezSecondarieTratto=DeltaSezPrincipaleCur+1

        xvals=numpy.linspace(0, tratto_length, NumSezSecondarieTratto)

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

        ArrayXSec=numpy.interp(xvals, dist_along, x0)
        ArrayYSec=numpy.interp(xvals, dist_along, y0)


        npt=len(ArrayXSec)

        pt_valle=[ArrayXSec[-1],ArrayYSec[-1]]

        Tooll_Distanza=tratto_length/float(npt)/4.0

        for i in range(1,npt-1):

            num_point_curr+=1

            pt_curr=ogr.Geometry(ogr.wkbPoint)
            pt_curr.AddPoint(ArrayXSec[i],ArrayYSec[i])

            c=pt_curr.GetPoint(0)
            PosizPtCur=isLeft(Point_a,Point_b,c)

            PosizPtValle=isLeft(Point_a,Point_b,pt_valle)

            NewGeom = SezioneInterpolata(pt_curr,cross_sec_line,next_cross_sec_line)
            Point_a=NewGeom.GetPoint(0)
            Point_b=NewGeom.GetPoint(1)

            Congruenza=PosizPtCur*PosizPtValle

            if Congruenza>0:

                UpstreaSez=ogr.CreateGeometryFromWkt(UpstreaSez_wkt)
                Dist1=pt_curr.Distance(UpstreaSez)
                Dist2=pt_curr.Distance(next_cross_sec_line)

                if Dist1>Tooll_Distanza and Dist2>Tooll_Distanza:

                    id_Points[num_point_curr]=pt_curr.GetPoint(0)
                    progr_Points[num_point_curr]=progr_along[i]

                    UpstreaSez_wkt=NewGeom.ExportToWkt()

                else:
                    txt='--  Discarded point n: %d of the intermediate sections downstream cross-section directional : %d-th due to distance less than the minimum' % (i,isec-1)
                    errMsg+='\n%s'% txt

                    pt_curr.Destroy()

            else:
                txt='--  Discarded point n: %d  of the intermediate sections downstream cross-section directional :%d-th due to a meander that goes upstream' % (i,isec-1)
                errMsg+='\n%s'% txt
                pt_curr.Destroy()

        # -----------------------------------------------------------
        # Creation of the point of the main downstream cross-section
        # ----------------------------------------------------------

        num_point_curr+=1

        PKUID_sec_id[ListaCrossSec[isec][0]]=num_point_curr
        listaIdSecPrincipali.append(num_point_curr)

        id_Points[num_point_curr]=pt_valle
        progr_Points[num_point_curr]=progr_along[-1]

        # azzero la memoria
        line_tratto.Destroy()

    # --------------------------------------
    # find the elevation values of the points
    # --------------------------------------

    # finding the progressive distances of the secondary sections
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

    dist = numpy.sqrt((ArrayXSec[:-1] - ArrayXSec[1:])**2 + (ArrayYSec[:-1] - ArrayYSec[1:])**2)
    dist_along_2 = numpy.concatenate(([0], dist.cumsum()))
    nn_along_2=len(dist_along_2)


    ArrayIdSecPrinc=numpy.array(listaIdSecPrincipali)


    QuoteAlveo=[]
    listax=[]
    listay=[]

    # searching for the minimum elevation
    # ...........................................................
    # also checking a certain number of cells around the point
    deltapixel=1
    steps=deltapixel*2+1
    finale=len(ArrayXSec)-1

    # dictionary upstream cross-section
    dic_SezMonte={}
    # dictionary downstream cross-section
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

        # update dictionaries
        if i in  listaIdSecPrincipali:
            dic_SezMonte[i]=-1
            dic_SezValle[i]=-1
        else:
            id_valle=numpy.searchsorted(ArrayIdSecPrinc, i)
            dic_SezMonte[i]=ArrayIdSecPrinc[id_valle-1]
            dic_SezValle[i]=ArrayIdSecPrinc[id_valle]

        # also look for elevations of a certain number of points around!
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
            # taking the minimum value
            z=ZetaArray.min()
            # calculate the average value with the point of the section
            z=(zpunto+z)/2.0
        else:
            z=0.0
            txt='err sez=%d' % i
            errMsg+='\n%s'% txt
        QuoteAlveo.append(z)

    # building a sloping riverbed
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


    # =================================
    # save the data in the geodatabase
    # =================================

    TargetTabella='PathPoints'

    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (TargetTabella.lower())
    cur.execute(sql)
    record=cur.fetchone()
    if record!=None:
        OriginEPSG=record[0]
    else:
        OriginEPSG=32632

    sql='SELECT PKUID,id FROM %s WHERE DamID=%d' % (TargetTabella,DamID)
    cur.execute(sql)
    ListaTratti=cur.fetchall()
    if len(ListaTratti)>0:
        # cancello i dati pregressi
        sql='DELETE FROM %s WHERE DamID=%d' % (TargetTabella,DamID)
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

        # save
        sql='INSERT INTO %s (DamID,id,type,elev,progr,geom) VALUES (%d'  %  (TargetTabella,DamID)
        sql+=',%d' %  id_cur
        sql+=',%d' %  tipo_i
        sql+=',%.2f' % zcur
        sql+=',%.3f' % progr_cur

        GeomWKT="GeomFromText('%s',%d)" % (pt_monte,OriginEPSG)
        sql+=',%s' % GeomWKT
        sql+=');'
        cur.execute(sql)

        pt.Destroy()

    conn.commit()

    # update the id of the main cross-sections
    # ----------------------------------------

    NomeTabella='MainCrossSec'

    for PKUID in PKUID_sec_id:
        ii=PKUID_sec_id[PKUID]
        sql='UPDATE %s SET id=%d WHERE DamID=%d AND PKUID=%d'  %  (NomeTabella,ii,DamID,PKUID)
        cur.execute(sql)

    conn.commit()
    # Close communication with the database
    cur.close()
    conn.close()


    return NotErr, errMsg

def SetCrossSec_2(mydb_path_user,PathFiles,ClipDEM,DamID):

    """
    Create the interpolated sections using the main cross-sections as guidelines

    It also creates the StreamDH.tif grid of heights from the river bed
    """

    NotErr=bool('True')
    errMsg='OK'


    # check for existence of input data
    # ---------------------------------
    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "There is no data for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    CrossSec=PathFiles+os.sep+'CrossSec.shp'

    if not os.path.exists(ClipDEM):
        errMsg = "Missing for the dam num =%s the clip DEM\nFirst make the clip of the digital terrain model !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    # =====================================
    # connection to the geodatabase
    # ====================================

    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
   # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    # check if the downstrema line exists
    NomeTabellaLinee='Downstreampath'

    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
    cur.execute(sql)

    record=cur.fetchone()
    if record!=None:
        DominioEPSG=record[0]
    else:
        DominioEPSG=32632

    dest_srs = ogr.osr.SpatialReference()
    dest_srs.ImportFromEPSG(DominioEPSG)

    sql='SELECT TotalLength,ST_AsText(geom) FROM %s WHERE DamID=%d' % (NomeTabellaLinee,DamID)
    cur.execute(sql)
    ChkDiga=cur.fetchone()

    if ChkDiga==None:
        errMsg = "In the table= %s there is no data for the dam num =%s \nPerform the downstream line calculation first !" % (NomeTabellaLinee,DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt_line=ChkDiga[1]
        TotalLength_1=ChkDiga[0]
        # creates the geometry of the river axis line
        StreamLine=ogr.CreateGeometryFromWkt(wkt_line)
        TotalLength=StreamLine.Length()

    # create a poly_stream of 50 m
    poly_stream=StreamLine.Buffer(50)

    TabellaPoligono='StudyArea'
    sql='SELECT Area,ST_AsText(geom) FROM %s WHERE DamID=%d' % (TabellaPoligono,DamID)
    cur.execute(sql)
    ChkDigaPoly=cur.fetchone()

    if ChkDigaPoly==None:
        errMsg = "In the table= %s there is no data for the dam num =%s \nFirst calculate the study area downstream !" % (TabellaPoligono,DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt_Dominio=ChkDigaPoly[1]
        TotalArea=ChkDigaPoly[0]

    CrossSecTot='%sTot.shp' % (CrossSec[:-4])
    CrossSecPoly='%sPoly.shp' % (CrossSec[:-4])

    # midline
    LineaMediana='%sMedianPath.shp' % (CrossSec[:-4])

    # polygons divided into right and left
    CrossSecPoly2='%sPoly_2.shp' % (CrossSec[:-4])

    # trace of the intermediate cross-sections representative of the reach
    CrossMedie='%sMean.shp' % (CrossSec[:-4])


    # trace of the interpolated sections and intermediate sections representative of the reach
    # need to create StreamDH.tif better defined !!
    Cross_Tot_Medie='%sTot_and_Mean.shp' % (CrossSec[:-4])

    # Point on the river of the representative section of the reach
    CrossMediePunto='%sMeanPunto.shp' % (CrossSec[:-4])

    # polygons right and left river
    PolySxDx=PathFiles + os.sep +'PolySxDx.shp'

    driver = ogr.GetDriverByName('ESRI Shapefile')

    # creating the domain Multipolygon
    DominoPoly=ogr.CreateGeometryFromWkt(wkt_Dominio)

    Area=DominoPoly.Area()

    # takes the largest polygon in the Multipoligons
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
        # takes the first one
        DominoRing=DominoPoly.GetGeometryRef(0)
        Area1=DominoRing.Area()

    # reading the geometry of the main cross-sections
    # ---------------------------------------------
    NomeTabella='MainCrossSec'


    sql='SELECT id,ST_AsText(geom) FROM %s WHERE DamID=%d' % (NomeTabella,DamID)
    cur.execute(sql)
    ListaCrossSec=cur.fetchall()

    n=len(ListaCrossSec)

    if n <= 0:

        errMsg = "There is no data for the dam num =%s \nFirst calculate the main cross-sections !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:

        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        cur.execute(sql)
        record=cur.fetchone()
        if record!=None:
            OriginEPSG=record[0]
        else:
            OriginEPSG=32632

        # load the geometries in a dictionary
        CrossSecGeom={}

        for rec in ListaCrossSec:
            CrossSecGeom[rec[0]]=rec[1]

        # reading of interpolation points
        # --------------------------------

        # type of section i-th
        NumSezTipo={}

        # Section characteristics
        NumSezElev={}
        NumSezProgr={}
        NumSezGeom={}

        # cross-section number upstream of i-th
        NumSezMonte={}
        # cross-section number downstream of i-th
        NumSezValle={}

        # current upstream section number
        NumMonte=-1

        # creo il dizionario dei punti a monte
        dic_pt_Monte={}
        NumPt_Monte=0


        TabellaPoints='PathPoints'

        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        cur.execute(sql)
        record=cur.fetchone()
        if record!=None:
            SourceEPSG=record[0]
        else:
            SourceEPSG=32632

        sql='SELECT id,type,elev,progr,ST_AsText(geom) FROM %s WHERE DamID=%d ORDER BY id' % (TabellaPoints,DamID)
        cur.execute(sql)
        ListaTratti=cur.fetchall()

        n=len(ListaTratti)

        if n <= 0:

            errMsg = "There is no data for the dam num =%s \nFirst calculate the main cross-sections !" % (DamID)
            NotErr= bool()
            return NotErr, errMsg

        for rec in ListaTratti:

            NumSez=rec[0]
            tipo=rec[1]
            NumSezTipo[NumSez]=rec[1]
            NumSezElev[NumSez]=rec[2]
            NumSezProgr[NumSez]=rec[3]
            NumSezGeom[NumSez]=rec[4]

            # update NumMonte
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
    # creating interpolated sections
    # ===============================

    # -----------------------------------------
    # Create the shapefile of the total sections
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

    # list of saved sections
    ListaNumSezTotOk=[]


     # river reach of the section point
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

                numMonte=NumSezMonte[NumSez]
                SezMonteGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[numMonte])
                wktValle=CrossSecGeom[NumSezValle[numMonte]]
                SezValleGeom=ogr.CreateGeometryFromWkt(wktValle)
                NewGeom = SezioneInterpolata(pt_curr,SezMonteGeom,SezValleGeom)

                try:
                    lengh1=NewGeom.Length()


                    # clip by domain
                    NewGeom1=NewGeom.Intersection(DominoRing)

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

    # loads the geometry of the sections to be clipped
    # ---------------------------------------------
    ds = driver.Open(CrossSecTot, 0)
    if ds is None:
        errMsg='Could not open file %s' % CrossSecTot
        NotErr= bool()
        return NotErr, errMsg


    layer = ds.GetLayer()

    n = layer.GetFeatureCount()

    # dictionary geometries
    dic_Sez_ToClip={}

    feat = layer.GetNextFeature()

    while feat:

        NumSez=feat.GetField('id')
        line_sez=feat.GetGeometryRef()
        wkt=line_sez.ExportToWkt()
        dic_Sez_ToClip[NumSez]=wkt

        feat = layer.GetNextFeature()

    ds.Destroy()


    # intersection control of total cross-sections
    # --------------------------------------------

    bufferClip=50.0

    distBufferClip=bufferClip*1.1

    ds = driver.Open(CrossSecTot, 1)
    if ds is None:
        errMsg='Could not open file %s' % CrossSecTot
        NotErr= bool()
        return NotErr, errMsg


    layer = ds.GetLayer()

    n = layer.GetFeatureCount()

    # dictionary where to save the intersection point
    dic_Sez_PointIntersect={}

    # dictionary where to save the buffer
    dic_Sez_Clip={}

    # dictionary of check if the buffer is on the right
    dic_Sez_lato_destro={}

    # List sections clipped at the first while loop
    lista_clipped=[]

    # List of downstream sections intersected at the first while loop
    lista_intersects=[]

    feat = layer.GetNextFeature()

    num_sez=len(ListaNumSezTotOk)

    while feat:

        NumSez_Monte=feat.GetField('id')
        line_monte_geom=feat.GetGeometryRef()
        ii_ini=ListaNumSezTotOk.index(NumSez_Monte)+1
        # current reach of the river line
        line_trattocorrente = ogr.CreateGeometryFromWkt(dic_Sez_tratto_fiume[NumSez_Monte])

        update_sez=bool()

        for ii in range(ii_ini,num_sez):

            NumSez_Valle=ListaNumSezTotOk[ii]
            wkt=dic_Sez_ToClip[NumSez_Valle]
            line_valle= ogr.CreateGeometryFromWkt(wkt)


            # check if the two lines do not intersect
            # ---------------------------------------------
            if line_valle.Crosses(line_monte_geom):

                update_sez=bool('True')

                # find the intersection point
                pt_interes_sez=line_monte_geom.Intersection(line_valle)

                # check if left or right
                # --------------------------------
                orientamento=ProdottoVettSezioni(line_monte_geom,line_valle)
                if  orientamento>0:
                    # an intersection is assumed on the left bank
                    inters_destra=bool()
                else:
                    inters_destra=bool('True')

                # creates a buffer around the point
                pt_poly = pt_interes_sez.Buffer(bufferClip)

                if NumSez_Valle not in  lista_intersects:
                    lista_intersects.append(NumSez_Valle)
                # salva il buffer e eventualmente sovrascrive
                dic_Sez_Clip[NumSez_Valle]= pt_poly.ExportToWkt()
                dic_Sez_PointIntersect[NumSez_Valle]= pt_interes_sez.ExportToWkt()
                dic_Sez_lato_destro[NumSez_Valle]= inters_destra

                # executes the clip of only the one upstream
                # the line is shortened so as not to touch the downstream section
                lines=line_monte_geom.Difference(pt_poly)

                # counts the number of geometries: there may be cases where eliminating
                # from the line the part intersected by the buffer is obtained as a result
                # two lines, one on the right and one on the left
                nn=lines.GetGeometryCount()
                if nn>1:
                    for i in range(nn):
                        g = lines.GetGeometryRef(i)
                        # choose the part of the line that intersects the axis of the river
                        if g.Crosses(line_trattocorrente):
                            line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            # check the direction
                            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                            break
                else:
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

                wkt_monte_new=line.ExportToWkt()
                line_monte_geom= ogr.CreateGeometryFromWkt(wkt_monte_new)

        # save changes in the shapefile
        if update_sez:
            feat.SetGeometry(line)
            layer.SetFeature(feat)

        # reading the next feature
        feat = layer.GetNextFeature()

    # execute any clip in the last section
    if NumSez_Monte in lista_intersects:

        wkt=dic_Sez_ToClip[NumSez_Monte]
        line_valle= ogr.CreateGeometryFromWkt(wkt)

        line_trattocorrente = ogr.CreateGeometryFromWkt(dic_Sez_tratto_fiume[NumSez])

        pt_interes_sez= ogr.CreateGeometryFromWkt(dic_Sez_PointIntersect[NumSez])

        # --------------------------------
        inters_destra= dic_Sez_lato_destro[NumSez]

        pt_poly=ogr.CreateGeometryFromWkt(dic_Sez_Clip[NumSez])

        lines=line_valle.Difference(pt_poly)

        nn=lines.GetGeometryCount()
        if nn>1:

            for i in range(nn):
                g = lines.GetGeometryRef(i)
                if g.Crosses(line_trattocorrente):
                    line_select=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                    line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,line_select)
                    break
        else:
            line=CheckVerso(distBufferClip,inters_destra,pt_interes_sez,lines)

        feat.SetGeometry(line)

        layer.SetFeature(feat)

    ds.Destroy()


    # saving the line of the last section
    # ----------------------------------

    ds = driver.Open(CrossSecTot, 0)
    if ds is None:
        errMsg='Could not open file %s' % CrossSecTot
        NotErr= bool()
        return NotErr, errMsg

    layer = ds.GetLayer()

    n = layer.GetFeatureCount()
    feat_end = layer.GetFeature(n-1)
    NumSez_end=feat_end.GetField('id')
    line_valle=feat_end.GetGeometryRef()
    wkt_sez_end=line_valle.ExportToWkt()

    ds.Destroy()

    # Also creates the dictionary of point coordinates
    dic_PointPath={}

    for NumPoint in NumSezTipo:

        pt_curr=ogr.CreateGeometryFromWkt(NumSezGeom[NumPoint])

        punto=pt_curr.GetPoint()

         # adds to the dictionary
        dic_PointPath[NumPoint]=punto

    # ---------------------------------------------------------
    # creates the river reach polygons
    # -------------------------------------
    # and creation of the global polygon verifying the intersection
    # ---------------------------------------------------------

    CreaPolygoni=1

    if CreaPolygoni>0:

        # save Dictionary with the list of polygon points
        DicPuntiPoligono={}

        # save Dictionary with the wkt geometry of the clipped final polygon
        DicPoligonoClip={}

        # save Dictionary with the list of the points of the river axis
        DicPuntiAsse={}

        # save the list of polygon numbers
        NumPolyList=[]

        # Creating the polygon merge of reachs to be used for the clip
        # --------------------------------------------------------------
        # initializing the reach polygon
        PoligonoUnioneTratti_0=ogr.Geometry(ogr.wkbPolygon)

        ds = driver.Open(CrossSecTot, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossSecTot
            NotErr= bool()
            return NotErr, errMsg

        Inlayer = ds.GetLayer()
        Spatialref = Inlayer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        SourceEPSG=int(Spatialref.GetAuthorityCode(None))

        shpnew3=CrossSecPoly
        if os.path.exists(shpnew3):
            driver.DeleteDataSource(shpnew3)

        outDS3 = driver.CreateDataSource(shpnew3)
        if outDS3 is None:
            errMsg='Could not create file %s' % shpnew3
            NotErr= bool()
            return NotErr, errMsg

        outLayer3 = outDS3.CreateLayer('CrossSecPoly', Spatialref,geom_type=ogr.wkbPolygon)

        # create the new id and Name fields in the output shapefile
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

                # create the polygon
                ring = ogr.Geometry(ogr.wkbLinearRing)
                line_monte=ogr.CreateGeometryFromWkt(line_monte_wkt)

                # check if the two lines do not intersect
                # ---------------------------------------------
                if not line_monte.Crosses(line_valle):


                    ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                    ring.AddPoint(line_monte.GetX(1),line_monte.GetY(1))
                    ring.AddPoint(line_valle.GetX(1),line_valle.GetY(1))
                    ring.AddPoint(line_valle.GetX(0),line_valle.GetY(0))
                    ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))

                    # saving the points
                    PuntiPoligono=[]
                    PuntiPoligono.append((line_monte.GetX(0),line_monte.GetY(0)))
                    PuntiPoligono.append((line_monte.GetX(1),line_monte.GetY(1)))
                    PuntiPoligono.append((line_valle.GetX(1),line_valle.GetY(1)))
                    PuntiPoligono.append((line_valle.GetX(0),line_valle.GetY(0)))
                    DicPuntiPoligono[NumSValle]=PuntiPoligono

                else:

                    pt_interes_sez=line_monte.Intersection(line_valle)

                    # check the distance to the left
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                    dist_sx=pt_interes_sez.Distance(point)
                    # check the distance to the right
                    point.Destroy()
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(line_monte.GetX(1),line_monte.GetY(1))
                    dist_dx=pt_interes_sez.Distance(point)
                    point.Destroy()

                    if dist_dx>dist_sx:
                        # case of intersection on the left
                        ring.AddPoint(pt_interes_sez.GetX(0),pt_interes_sez.GetY(0))
                        ring.AddPoint(line_monte.GetX(1),line_monte.GetY(1))
                        ring.AddPoint(line_valle.GetX(1),line_valle.GetY(1))
                        ring.AddPoint(pt_interes_sez.GetX(0),pt_interes_sez.GetY(0))
                        # saving the points
                        PuntiPoligono=[]
                        PuntiPoligono.append((pt_interes_sez.GetX(0),pt_interes_sez.GetY(0)))
                        PuntiPoligono.append((line_monte.GetX(1),line_monte.GetY(1)))
                        PuntiPoligono.append((line_valle.GetX(1),line_valle.GetY(1)))
                        PuntiPoligono.append((pt_interes_sez.GetX(0),pt_interes_sez.GetY(0)))
                        DicPuntiPoligono[NumSValle]=PuntiPoligono

                    else:
                        # case of intersection on the right
                        ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                        ring.AddPoint(pt_interes_sez.GetX(0),pt_interes_sez.GetY(0))
                        ring.AddPoint(line_valle.GetX(0),line_valle.GetY(0))
                        ring.AddPoint(line_monte.GetX(0),line_monte.GetY(0))
                        # saving the points
                        PuntiPoligono=[]
                        PuntiPoligono.append((line_monte.GetX(0),line_monte.GetY(0)))
                        PuntiPoligono.append((pt_interes_sez.GetX(0),pt_interes_sez.GetY(0)))
                        PuntiPoligono.append((line_valle.GetX(0),line_valle.GetY(0)))
                        PuntiPoligono.append((line_monte.GetX(0),line_monte.GetY(0)))
                        DicPuntiPoligono[NumSValle]=PuntiPoligono

                    pt_interes_sez.Destroy()

                # saving the number in the list
                NumPolyList.append(NumSValle)

                poly_new=ogr.Geometry(ogr.wkbPolygon)

                poly_new.AddGeometry(ring)

                numpt1=poly_new.GetGeometryRef(0).GetPointCount()

                # intersects the polygon of the domain
                NewPoly=poly_new.Intersection(DominoRing)

                if NewPoly!=None:
                    npoly=NewPoly.GetGeometryCount()
                    if  npoly>0:
                        numpt=NewPoly.GetGeometryRef(0).GetPointCount()
                        # intersection with the axis of the river
                        TrattoStreamLine= StreamLine.Intersection(NewPoly)
                        nntratti=TrattoStreamLine.GetGeometryCount()
                        if nntratti>0:
                            # to avoid problems of reach with handles, use the straight line that joins
                            # the upstream point with the downstream point
                            PuntiTratto=[]
                            PuntiTratto.append(dic_PointPath[NumMonte])
                            PuntiTratto.append(dic_PointPath[NumSValle])
                        else:
                            TrattoStream_max=TrattoStreamLine
                            # extracting the points
                            PuntiTratto=[]
                            for iii in range(0, TrattoStream_max.GetPointCount()):
                                # GetPoint returns a tuple not a Geometry
                                pt = TrattoStream_max.GetPoint(iii)
                                PuntiTratto.append(pt)

                        DicPuntiAsse[NumSValle]=PuntiTratto

                        poly_wkt=NewPoly.ExportToWkt()
                    else:
                        # geometry does not exist !!
                        poly_wkt=NewPoly.ExportToWkt()
                        # using the one before the intersection
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


                        # to avoid problems of reach with handles, use the straight line that joins
                        # the upstream point with the downstream point
                        PuntiTratto=[]
                        PuntiTratto.append(dic_PointPath[NumMonte])
                        PuntiTratto.append(dic_PointPath[NumSValle])

                        DicPuntiAsse[NumSValle]=PuntiTratto

                else:

                    # cases of missed intersection ??????
                    NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                    npoly=poly_new.GetGeometryCount()
                    geom_poly=poly_new.GetGeometryRef(0)
                    numpt=geom_poly.GetPointCount()

                    poly_wkt=NewPoly.ExportToWkt()
                    # save the data
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


                    PuntiTratto=[]
                    PuntiTratto.append(dic_PointPath[NumMonte])
                    PuntiTratto.append(dic_PointPath[NumSValle])

                    DicPuntiAsse[NumSValle]=PuntiTratto

                # control of the intersection with the upstream parts
                NewPolyDiff=NewPoly.Difference(PoligonoUnioneTratti_0)

                # counts the number of geometries
                nn=NewPolyDiff.GetGeometryCount()
                if nn>1:
                    poly_centro = ogr.CreateGeometryFromWkt(dic_Sez_tratto_fiume[NumMonte])
                    for i in range(nn):
                        g = NewPolyDiff.GetGeometryRef(i)
                        wkt_debug= g.ExportToWkt()

                        if g.Intersect(poly_centro):
                            NewPoly_0=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            NewPoly=NewPoly_0.Buffer(0.1)
                            break
                else:
                    NewPoly=NewPolyDiff.Buffer(0.1)

                # union
                nn0=PoligonoUnioneTratti_0.GetGeometryCount()

                PoligonoUnioneTratti_0=PoligonoUnioneTratti_0.Union(NewPoly)

                # check the cases with holes
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

                # save the clipped polygon
                DicPoligonoClip[NumSValle]= NewPoly.ExportToWkt()

                feature.SetGeometry(NewPoly)
                outLayer3.CreateFeature(feature)
                ring.Destroy()
                poly_new.Destroy()
                NewPoly.Destroy()

                # update for the next step

                line_monte_wkt=wkt

                NumMonte=NumSValle

        # close the connection
        ds.Destroy()

        outDS3.Destroy()

        PoligonoUnioneTratti_Wkt=PoligonoUnioneTratti_0.ExportToWkt()


    # Shapefile middle cross-section
    CreaSezioniMediane=1

    if CreaSezioniMediane>0:

        # it is assumed that you have already stored NumSezMonte and NumSezValle
        n = len(ListaTratti)

        # save the buffer of the intersection point of the middle cross-sections
        # with the axis of the river
        dic_SezMediane_Clip={}


        # creates new shapefile of section points on the river
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

        # creates a new section shapefile
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

        # create the new id and Name fields in the output shapefile
        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer2.CreateField(fieldDefn2)
        outLayer3.CreateField(fieldDefn2)
        # create field midpoint distance
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

        # reading the first feature
        NumSez=ListaTratti[0][0]
        ProgrMonte=NumSezProgr[NumSez]
        Z_monte= NumSezElev[NumSez]

        pt_monte=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

        # creating the first point
        feature2 = ogr.Feature(featureDefn2)
        feature2.SetField('id', NumSez)
        feature2.SetGeometry(pt_monte)
        outLayer2.CreateFeature(feature2)

        # creating the first cross-section
        # --------------------------------
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

        # clipping along the domain
        pt1=pt_monte.GetPoint()
        pt_poly= pt_monte.Buffer(bufferClip)

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


        # save the buffer
        dic_SezMediane_Clip[NumSez]= pt_poly.ExportToWkt()

        numtratti=len(ListaTratti)

        for itratto in range(1,numtratti):

            NumSez=ListaTratti[itratto][0]

            # skip unsaved sections
            if NumSez in ListaNumSezTotOk:

                tipo=NumSezTipo[NumSez]
                ProgrValle=NumSezProgr[NumSez]
                Z_valle= NumSezElev[NumSez]

                pt_valle=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

                pt2=pt_valle.GetPoint()
                pt_mean=PuntoMedio(pt1,pt2)
                pt_curr=ogr.Geometry(ogr.wkbPoint)
                pt_curr.AddPoint(pt_mean[0],pt_mean[1])

                pt_poly = pt_curr.Buffer(bufferClip)
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
                    # case upstream a main cross-section
                    SezValleGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[NumSez])
                    # trova la sezione di monte
                    numMonte=NumSezMonte[NumSez]
                    SezMonteGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[numMonte])

                else:
                   # case of an intermediate section
                    numMonte=NumSezMonte[NumSez]
                    SezMonteGeom=ogr.CreateGeometryFromWkt(CrossSecGeom[numMonte])
                    wktValle=CrossSecGeom[NumSezValle[numMonte]]
                    SezValleGeom=ogr.CreateGeometryFromWkt(wktValle)

                # calculate the interpolate section
                try:
                    NewGeom = SezioneInterpolata(pt_curr,SezMonteGeom,SezValleGeom)

                    try:
                        lengh1=NewGeom.Length()

                        # cut along the polygon of the current reach
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

                    pt_sinistra=NewGeom.GetPoint(0)
                    Geom_pt_sinistra=ogr.Geometry(ogr.wkbPoint)
                    Geom_pt_sinistra.AddPoint(pt_sinistra[0],pt_sinistra[1])

                    pt_intersect= NewGeom.Intersection(StreamLine)
                    if pt_intersect!=None:
                        n_inters=pt_intersect.GetGeometryCount()
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

                    # distance from the point of intersection with the axis of the river
                    distanza=Geom_pt_sinistra.Distance(ptt)
                    feature3.SetField('dist1', distanza)

                    feature2.SetGeometry(pt_curr)
                    outLayer2.CreateFeature(feature2)

                    outLayer3.CreateFeature(feature3)

                    NewGeom.Destroy()
                    pt_curr.Destroy()

                except:
                    pass

                # next step
                pt1=pt2
                ProgrMonte=ProgrValle*1.0
                Z_monte=Z_valle*1.0

            else:
                txt='Excluding sez: %d' % NumSez
                errMsg+='\n%s'% txt

        # adding the last section
        # =========================

        NumSez=ListaTratti[-1][0]

        if NumSez in ListaNumSezTotOk:

            ProgrMonte=NumSezProgr[NumSez]
            Z_monte= NumSezElev[NumSez]

            pt_monte=ogr.CreateGeometryFromWkt(NumSezGeom[NumSez])

            feature2 = ogr.Feature(featureDefn2)
            feature2.SetField('id', NumSez)
            feature2.SetGeometry(pt_monte)
            outLayer2.CreateFeature(feature2)

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

            dic_SezMediane_Clip[NumSez]= pt_poly.ExportToWkt()


        # =========================
        outDS2.Destroy()
        outDS3.Destroy()

        # Check
        # --------------------------------------
        ds = driver.Open(CrossMedie, 1)
        if ds is None:
            errMsg='Could not open file %s' % CrossMedie
            NotErr= bool()
            return NotErr, errMsg


        layer = ds.GetLayer()

        n = layer.GetFeatureCount()

        dic_Sez_Clip={}
        lista_clipped=[]

        feat = layer.GetNextFeature()

        line_monte=feat.GetGeometryRef()
        line_monte_wkt=line_monte.ExportToWkt()

        NumSez_Monte=feat.GetField('id')

        feat = layer.GetNextFeature()

        while feat:

            NumSez_Valle=feat.GetField('id')

            line_valle=feat.GetGeometryRef()

            wkt=line_valle.ExportToWkt()

            line_monte_geom=ogr.CreateGeometryFromWkt(line_monte_wkt)

            if line_valle.Crosses(line_monte_geom):

                pt_interes_sez=line_monte_geom.Intersection(line_valle)
                pt_poly = pt_interes_sez.Buffer(bufferClip)

                if  NumSez_Monte not in lista_clipped:
                    dic_Sez_Clip[NumSez_Monte]=pt_poly.ExportToWkt()

                lista_clipped.append(NumSez_Valle)
                line_monte_wkt= line_valle.ExportToWkt()

                lines=line_valle.Difference(pt_poly)

                nn=lines.GetGeometryCount()
                if nn>1:
                    poly_centro= ogr.CreateGeometryFromWkt(dic_SezMediane_Clip[NumSez_Valle])
                    for i in range(nn):
                        g = lines.GetGeometryRef(i)
                        if g.Crosses(poly_centro):
                            line=ogr.CreateGeometryFromWkt(g.ExportToWkt())
                            break
                else:
                    line=lines

                feat.SetGeometry(line)
                layer.SetFeature(feat)

            line_monte_wkt=wkt
            NumSez_Monte= NumSez_Valle

            feat = layer.GetNextFeature()

        layer.ResetReading()

        feat = layer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')

            if NumSez in dic_Sez_Clip:

                line_valle=feat.GetGeometryRef()

                pt_poly=ogr.CreateGeometryFromWkt(dic_Sez_Clip[NumSez])

                lines=line_valle.Difference(pt_poly)

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

            feat = layer.GetNextFeature()

        ds.Destroy()

        # end check

    CreaMediana=1

    if CreaMediana>0:

        shpnew4=LineaMediana
        if os.path.exists(shpnew4):
            driver.DeleteDataSource(shpnew4)

        outDS4 = driver.CreateDataSource(shpnew4)
        if outDS4 is None:

            errMsg='Could not create file %s' % shpnew4
            NotErr= bool()
            return NotErr, errMsg

        outLayer4 = outDS4.CreateLayer('LineaMediana', dest_srs,geom_type=ogr.wkbLineString)

        fieldDefn4 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer4.CreateField(fieldDefn4)

        featureDefne4 = outLayer4.GetLayerDefn()

        feature = ogr.Feature(featureDefne4)
        feature.SetField('id', 1)

        line4=ogr.Geometry(ogr.wkbLineString)

        n = len(ListaTratti)

        NumMonte=0

        for NumPoint in NumSezTipo:

            pt_curr=ogr.CreateGeometryFromWkt(NumSezGeom[NumPoint])

            punto=pt_curr.GetPoint()

            dic_PointPath[NumPoint]=punto

            if NumPoint in ListaNumSezTotOk:
                line4.AddPoint(punto[0],punto[1])

        feature.SetGeometry(line4)
        outLayer4.CreateFeature(feature)

        lineaMedianaWkt=line4.ExportToWkt()
        line4.Destroy()

        outDS4.Destroy()


    # Create : Cross_Tot_Medie
    # ---------------------------
    CreaCross_Tot_Medie=1

    if CreaCross_Tot_Medie>0:

        ds = driver.Open(CrossSecTot, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossSecTot
            NotErr= bool()
            return NotErr, errMsg

        Inlayer = ds.GetLayer()
        Spatialref = Inlayer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        SourceEPSG=int(Spatialref.GetAuthorityCode(None))

        shpnew3=Cross_Tot_Medie
        if os.path.exists(shpnew3):
            driver.DeleteDataSource(shpnew3)

        outDS3 = driver.CreateDataSource(shpnew3)
        if outDS3 is None:
            errMsg='Could not create file %s' % shpnew3
            NotErr= bool()
            return NotErr, errMsg

        outLayer3 = outDS3.CreateLayer('SezioniTotMedie', dest_srs,geom_type=ogr.wkbLineString)

        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer3.CreateField(fieldDefn2)

        # create zmin field
        fieldDefnZmin = ogr.FieldDefn('zmin', ogr.OFTReal)
        outLayer3.CreateField(fieldDefnZmin)

        featureDefn3 = outLayer3.GetLayerDefn()

        feat = Inlayer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')
            z_curr=float(feat.GetField('zmin'))

            feature = ogr.Feature(featureDefn3)
            feature.SetField('id', NumSez)
            feature.SetField('zmin', z_curr)

            linea=feat.GetGeometryRef()
            feature.SetGeometry(linea)
            outLayer3.CreateFeature(feature)


            feat = Inlayer.GetNextFeature()

        ds.Destroy()

        # adding the intermediate sections
        ds = driver.Open(CrossMedie, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossMedie
            NotErr= bool()
            return NotErr, errMsg

        Inlayer = ds.GetLayer()

        feat = Inlayer.GetNextFeature()

        feat = Inlayer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')
            z_curr=float(feat.GetField('zmin'))

            feature = ogr.Feature(featureDefn3)
            feature.SetField('id', -NumSez)
            feature.SetField('zmin', z_curr)

            linea=feat.GetGeometryRef()
            feature.SetGeometry(linea)
            outLayer3.CreateFeature(feature)

            feat = Inlayer.GetNextFeature()


        ds.Destroy()
        # save
        outDS3.Destroy()



    # creates the polygons of the river reachs divided into right and left
    # ----------------------------------------------------------------
    CreaPolygoni2=1

    if CreaPolygoni2>0:

        ds = driver.Open(CrossSecPoly, 0)
        if ds is None:
            errMsg='Could not open file %s' % CrossSecPoly
            NotErr= bool()
            return NotErr, errMsg


        Inlayer = ds.GetLayer()
        Spatialref = Inlayer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        SourceEPSG=int(Spatialref.GetAuthorityCode(None))


        # create new shapefile
        shpnew5=CrossSecPoly2
        if os.path.exists(shpnew5):
            driver.DeleteDataSource(shpnew5)

        outDS5 = driver.CreateDataSource(shpnew5)
        if outDS5 is None:
            errMsg='Could not create file %s' % shpnew5
            NotErr= bool()
            return NotErr, errMsg

        outLayer5 = outDS5.CreateLayer('CrossSecPoly_2', Spatialref,geom_type=ogr.wkbPolygon)

        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn2)
        fieldDefn5 = ogr.FieldDefn('lato', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn5)

        featureDefn5 = outLayer5.GetLayerDefn()

        feat = Inlayer.GetNextFeature()

        while feat:

            NumSez=feat.GetField('id')

            PuntiPoligono=DicPuntiPoligono[NumSez]
            PuntiTratto=DicPuntiAsse[NumSez]

            Poly_curr=ogr.CreateGeometryFromWkt(DicPoligonoClip[NumSez])

            # create the two polygons
            for lato in range(2):

                feature = ogr.Feature(featureDefn5)
                feature.SetField('id', NumSez)
                feature.SetField('lato', lato)

                ring = ogr.Geometry(ogr.wkbLinearRing)

                if lato==0:

                    # left upstream point
                    pt=PuntiPoligono[0]
                    ring.AddPoint(pt[0],pt[1])

                    # the points of the river axis are added
                    # ..........................................
                    for pt in PuntiTratto:
                        ring.AddPoint(pt[0],pt[1])

                    # left downstream point
                    pt=PuntiPoligono[3]
                    ring.AddPoint(pt[0],pt[1])

                    # left upstream point
                    pt=PuntiPoligono[0]
                    ring.AddPoint(pt[0],pt[1])

                    ring.CloseRings()
                    poly_new=ogr.Geometry(ogr.wkbPolygon)
                    poly_new.AddGeometry(ring)

                else:

                    # central upstream point
                    # using the first point of reach
                    pt= PuntiTratto[0]
                    ring.AddPoint(pt[0],pt[1])

                    # right upstream point
                    pt=PuntiPoligono[1]
                    ring.AddPoint(pt[0],pt[1])

                    # left downstream point
                    pt=PuntiPoligono[2]
                    ring.AddPoint(pt[0],pt[1])

                    # adding points of the river axis from downstream to upstream
                    for iii in range(len(PuntiTratto)-1,-1,-1):
                        pt=PuntiTratto[iii]
                        ring.AddPoint(pt[0],pt[1])

                    ring.CloseRings()

                    poly_new=ogr.Geometry(ogr.wkbPolygon)
                    poly_new.AddGeometry(ring)

                # intersects the polygon of the current reach
                NewPoly=poly_new.Intersection(Poly_curr)

                if NewPoly!=None:
                    npoly=NewPoly.GetGeometryCount()
                    if  npoly>0:
                        numpt=NewPoly.GetGeometryRef(0).GetPointCount()
                    else:
                        # No geometry !!
                        poly_wkt=NewPoly.ExportToWkt()
                        # using the one before the intersection
                        NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                        npoly=poly_new.GetGeometryCount()
                        geom_poly=poly_new.GetGeometryRef(0)
                        numpt=geom_poly.GetPointCount()

                        poly_wkt=NewPoly.ExportToWkt()
                        # save the data
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

                    # cases of missed intersection ??????
                    NewPoly=ogr.CreateGeometryFromWkt(poly_new.ExportToWkt())
                    npoly=poly_new.GetGeometryCount()
                    geom_poly=poly_new.GetGeometryRef(0)
                    numpt=geom_poly.GetPointCount()

                    poly_wkt=NewPoly.ExportToWkt()
                    # save
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

        outDS5.Destroy()
        ds.Destroy()

    # polygon creation left and right river
    # -------------------------------------------
    CreaPolygoni3=1

    if CreaPolygoni3>0:

        # Creating the polygon, joining the reachs to be used for the river clip
        # ----------------------------------------------------------------------
        PoligonoUnioneTratti_0=ogr.CreateGeometryFromWkt(PoligonoUnioneTratti_Wkt)


        # create the list of points in the streamLine
        # ============================================
        End_cross_sec_line=ogr.CreateGeometryFromWkt(wkt_sez_end)
        # intersects the polygon union of reachs
        NewStreamLine=StreamLine.Intersection(PoligonoUnioneTratti_0)
        if  NewStreamLine!=None:
            ListaPuntiStreamLine= SetStreamLinePoints(NewStreamLine,End_cross_sec_line)
        else:
            ListaPuntiStreamLine= SetStreamLinePoints(StreamLine,End_cross_sec_line)

        # create the new shapefile
        shpnew6=PolySxDx
        if os.path.exists(shpnew6):
            driver.DeleteDataSource(shpnew6)

        outDS5 = driver.CreateDataSource(shpnew6)
        if outDS5 is None:
            errMsg='Could not open file %s' % shpnew6
            NotErr= bool()
            return NotErr, errMsg

        outLayer5 = outDS5.CreateLayer('PolySxDx', Spatialref,geom_type=ogr.wkbPolygon)

        fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn2)
        fieldDefn5 = ogr.FieldDefn('lato', ogr.OFTInteger)
        outLayer5.CreateField(fieldDefn5)

        featureDefn5 = outLayer5.GetLayerDefn()


        # creating the left polygon
        # ---------------------------
        feature = ogr.Feature(featureDefn5)
        feature.SetField('id', 0)
        feature.SetField('lato', 0)

        ring = ogr.Geometry(ogr.wkbLinearRing)

        numPunti=len(NumPolyList)

        # left side
        for ii in range(numPunti):
            PuntiPoligono=DicPuntiPoligono[NumPolyList[ii]]
            # left upstream point
            pt=PuntiPoligono[0]
            ring.AddPoint(pt[0],pt[1])

        # adding the left point of the last section index 3
        PuntiPoligono=DicPuntiPoligono[NumPolyList[-1]]
        pt=PuntiPoligono[3]
        ring.AddPoint(pt[0],pt[1])

        # adding points of the river axis from downstream to upstream
        for iii in range(len(ListaPuntiStreamLine)-1,-1,-1):
            pt=ListaPuntiStreamLine[iii]
            ring.AddPoint(pt[0],pt[1])

        ring.CloseRings()

        poly_sx=ogr.Geometry(ogr.wkbPolygon)
        poly_sx.AddGeometry(ring)

        # intersects the polygon of the domain
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

        # create the right polygon
        # -------------------------
        feature = ogr.Feature(featureDefn5)
        feature.SetField('id', 1)
        feature.SetField('lato', 1)

        ring = ogr.Geometry(ogr.wkbLinearRing)

        # adding the line of the river from upstream to the downstream
        for iii in range(len(ListaPuntiStreamLine)):
            pt=ListaPuntiStreamLine[iii]
            ring.AddPoint(pt[0],pt[1])

        # right side from downstream to upstream: index 2
        for ii in range(numPunti-1,-1,-1):
            PuntiPoligono=DicPuntiPoligono[NumPolyList[ii]]
            pt=PuntiPoligono[2]
            ring.AddPoint(pt[0],pt[1])

        # adds the last point
        PuntiPoligono=DicPuntiPoligono[NumPolyList[0]]
        pt=PuntiPoligono[1]
        ring.AddPoint(pt[0],pt[1])

        ring.CloseRings()

        poly_dx=ogr.Geometry(ogr.wkbPolygon)
        poly_dx.AddGeometry(ring)

        # intersects the polygon of the domain
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


    # Creating the StreamDH.tif
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

        # opens  files with cross-sections
        # ------------------------------
        inDS1 = driver.Open(Cross_Tot_Medie, 0)
        if inDS1 is None:
            errMsg='Could not open ' + CrossSecTot
            NotErr= bool()
            return NotErr, errMsg

        InlayerCurve = inDS1.GetLayer()

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

            # generates a grid with the minimum elevation of the sections
            GridSez=np.zeros((rowsElev,colsElev),np.float32)

            format = 'MEM'
            type = GDT_Float32

            driver2 = gdal.GetDriverByName(format)
            driver2.Register()

            gt=indatasetElev.GetGeoTransform()

            ds = driver2.Create('GridSez', indatasetElev.RasterXSize, indatasetElev.RasterYSize, 1, type)
            if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
                ds.SetGeoTransform(gt)

            if prj is not None and len(prj) > 0:
                ds.SetProjection(prj)
            else:
                prj= spatialRef.ExportToWkt()
                ds.SetProjection(prj)

            iBand=1
            testo="ATTRIBUTE=%s"  % (nomecampoLivello)
            # Rasterize
            outband = ds.GetRasterBand(iBand)

            outband.WriteArray(GridSez, 0, 0)
            CampoValore=[testo]

            # creating the map of values
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

            # output file name
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
            spatialRef.ImportFromEPSG(32632)

        terreno = inband.ReadAsArray(0, 0, cols, rows).astype(numpy.float)
        mask_Nodata= terreno==inNoData


        # rasterize the polygons
        # ======================
        orig_data_source = ogr.Open(CrossSecPoly2)
        source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_source, "")
        source_layer = source_ds.GetLayer()


        format = 'Gtiff'
        type = GDT_Int16

        driver3 = gdal.GetDriverByName(format)
        driver3.Register()

        TrattiRaster=PathFiles+os.sep+'DestraSinistra.tif'

        dsRaster = driver3.Create(TrattiRaster, cols, rows, 1, type)
        gt1=geotransform
        if gt1 is not None and gt1 != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            dsRaster.SetGeoTransform(gt1)

        if prj is not None and len(prj) > 0:
            dsRaster.SetProjection(prj)
        else:
            prj= spatialRef.ExportToWkt()
            dsRaster.SetProjection(prj)

        # Rasterize
        iBand=1
        outband = dsRaster.GetRasterBand(iBand)

        outNodata=-9999

        # writes -1 to the whole array
        ClassTratti=numpy.zeros((rows,cols)).astype(numpy.int)
        ClassTratti=ClassTratti-1

        outband.WriteArray(ClassTratti, 0, 0)

        # Rasterize
        err = gdal.RasterizeLayer(dsRaster, [1], source_layer,
                burn_values=[0],
                options=["ATTRIBUTE=lato"])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)

        # writes Nodata
        MatriceDati=outband.ReadAsArray(0, 0, cols, rows)
        MatriceDati=numpy.choose(mask_Nodata,(MatriceDati,outNodata))
        outband.WriteArray(MatriceDati, 0, 0)

        outband.FlushCache()
        outband.SetNoDataValue(outNodata)
        outband.GetStatistics(0,1)

        outband=None

        dsRaster=None
        # closes the data sources
        orig_data_source.Destroy()



    return NotErr, errMsg

if __name__ == '__main__':


    mydb_path_user='..'+ os.sep+'db'+os.sep+'USER_GeoDB.sqlite'

    # San Giuliano
    DamID=449
    PathFiles='..'+ os.sep+ str(DamID)
    fileDEM=  PathFiles+ os.sep+'DTM_clip.tif'

    NotErr, errMsg= SetIntermediatePoints(mydb_path_user,PathFiles,fileDEM,DamID,)

    print(NotErr,errMsg)

    NotErr, errMsg= SetCrossSec_2(mydb_path_user,DamID,PathFiles,fileDEM,DamID)

    print(NotErr,errMsg)
