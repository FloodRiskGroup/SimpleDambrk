# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SimpleDambrk - Create Flooded Area Polygon

                                 A QGIS plugin
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

The script create:

    - the feature of Floading Areas
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

import math
import csv
import scipy.interpolate as il
from scipy import signal

from datetime import datetime

def FloodSeverity(DV):

    # guidance from DSO-99-06 on selecting the flood severity category
    # DV in m2/s

    if DV<=4.6:
        fldsev='Low'
    elif DV>4.6 and DV<15.0:
        fldsev='Medium'
    else:
        fldsev='High'

    return fldsev

def ConseqFactot(FloodSeverity,Time):

    # guidance from DSO-99-06
    # Fatality Rate : Fraction of people at risk projected to die

    if 'Low' in FloodSeverity:
        if Time<15:
            Factor=0.01
        elif Time>=15 and Time<=60:
            Factor=0.005
        else:
            Factor=0.0003
    elif 'Medium' in FloodSeverity:
        if Time<15:
            Factor=0.15
        elif Time>=15 and Time<=60:
            Factor=0.03
        else:
            Factor=0.02
    else:
        if Time<15:
            Factor=0.75
        elif Time>=15 and Time<=60:
            Factor=0.4
        else:
            Factor=0.2
    return Factor


def PuntoIntermedio(pt1,pt2,percent):

    # calculates the intermediate point at a certain percentage between two
    pt_intermedio=[pt1[0]+percent*(pt2[0]-pt1[0]),pt1[1]+percent*(pt2[1]-pt1[1])]

    return pt_intermedio

def GridFloodingAreas(mydb_path_user,PathFiles,DamID,UseEnergyHead):

    NotErr=bool('True')
    errMsg='OK'
    MatriceRisultati=[]


    # ---------------------------------
    PathFiles=os.path.realpath(PathFiles)

    mydb_path_user=os.path.realpath(mydb_path_user)

    # polygon floodable area
    AreaInondabile=PathFiles+os.sep+'AreaInondabile_tot.shp'

    # Polygon Areas1: first component of the floodable area
    # ------------------------------------------------------
    # area on the right and left of the evaluated river axis
    # based on the width in the right and left obtained
    # from model propagation calculation
    # one-dimensional
    AreaInondabile_1=PathFiles+os.sep+'AreaInondabile_1.shp'

    if not os.path.exists(PathFiles):
        errMsg = "There is no data for the dam num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg
    else:
        os.chdir(PathFiles)
        log_file=open('log.txt','w')
        timenow_our=datetime.now().strftime('%y-%m-%d %H:%M')
        log_file.write('Start %s\n' % timenow_our)
        log_file.close()

    # trace intermediate sections representative of the section
    CrossMedie=PathFiles+os.sep+'CrossSecMean.shp'
    if not os.path.exists(CrossMedie):
        errMsg = "Missing CrossSecMean.shp for the dam num =%s \nPerform the calculation of the downstream sections first!" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDHFilled.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Missing for the dam num =%s il StreamDH\nCarry out first ModificaDH !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    # poligoni tratti
    CrossSecPoly=PathFiles+os.sep+ 'CrossSecPoly.shp'
    if not os.path.exists(CrossSecPoly):
        errMsg = "Missing CrossSecPoly.shp for the dam num =%s \nPerform the calculation of the downstream sections first!" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    # poligoni tratti divisi in destra e sinistra
    CrossSecPoly_2=PathFiles+os.sep+ 'CrossSecPoly_2.shp'
    if not os.path.exists(CrossSecPoly_2):
        errMsg = "Missing CrossSecPoly_2.shp for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    # i due poligoni sinistra e destra idraulica
    PolySxDx=PathFiles + os.sep +'PolySxDx.shp'
    if not os.path.exists(PolySxDx):
        errMsg = "Missing PolySxDx.shp for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    FileMatricePixDestra=PathFiles+os.sep+'MatricePixDestra.csv'
    if not os.path.exists(FileMatricePixDestra):
        errMsg = "Missing for the dam num =%s MatricePixDestra.csv\nPerform CreaCurveAreaAltezza first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    # =======================================
    # Reading the characteristics of the GRID
    # =======================================

    gdal.AllRegister()

    indataset = gdal.Open(StreamDH, GA_ReadOnly )
    if indataset is None:
        errMsg='Could not open file %s' % StreamDH
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

    inband=None
    indataset=None

    # --------------------
    # reading from the database
    # --------------------
    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)

    cur= conn.cursor()

    # reading the data for the calculation of flooding widths
    NomeTabella='Q_H_max'

    sql='SELECT '
    sql+=' PixDist'
    sql+=', Progr_fiume'
    sql+=', Qmax'
    sql+=', Hmax'
    sql+=', Bmax'
    sql+=', Vmax'
    sql+=', Time'
    sql+=' FROM %s' % NomeTabella
    sql+=' WHERE DamID=%d' % (DamID)
    sql+=' ORDER BY PixDist;'
    cur.execute(sql)
    MatriceDati=cur.fetchall()

    if len(MatriceDati)==0:
        errMsg = "Missing for the dam num =%s data Q_H_max\nCarry out first Calculation of propagation !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    ListaTratti=[]
    Progr_fiume=[]
    Qmax=[]
    Hmax=[]
    Bmax=[]
    Vmax=[]
    Time=[]

    for row in MatriceDati:
        ListaTratti.append(int(row[0]))
        Progr_fiume.append(float(row[1]))
        Qmax.append(float(row[2]))
        Hmax.append(float(row[3]))
        Bmax.append(float(row[4]))
        Vmax.append(float(row[5]))
        Time.append(float(row[6]))

    #  array
    Progr_fiume_array=numpy.array(Progr_fiume,dtype =numpy.float)
    Qmax_array=numpy.array(Qmax,dtype =numpy.float)
    Hmax_array=numpy.array(Hmax,dtype =numpy.float)
    Bmax_array=numpy.array(Bmax,dtype =numpy.float)
    Vmax_array=numpy.array(Vmax,dtype =numpy.float)
    Time_array=numpy.array(Time,dtype =numpy.float)

    # finding the maximum depth value
    Hmax_tot=Hmax_array.max()

    # reading of the curves necessary to evaluate the shift to the right of the flood area
    fin=open(FileMatricePixDestra,'r')
    reader = csv.reader(fin,delimiter=';')
    try:
        # python 2.7
        headers = reader.next()
    except:
        # python 3.4
        headers = reader.__next__()

    nn=len(headers)
    hvals=[]
    Vettoreh=[]
    for i in range(1,nn):
        hvals.append(headers[i])
        pp=headers[i].split('=')
        Vettoreh.append(float(pp[1]))

    # water depth array
    H_Array=numpy.array(Vettoreh,dtype =numpy.float)

    # dictionary of section numbers
    dic_PixDist={}
    # matrix of the quantities
    MatricePix=[]
    ii=-1
    for row in reader:
        ii+=1
        dic_PixDist[int(row[0])]=ii
        MatricePix.append(row[1:])

    fin.close()

    # matrix of the percentage of area on the right bank of the river for each height
    MatriceArray=numpy.array(MatricePix,dtype =numpy.float)


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

    if len(MatriceDati)==0:
        errMsg = "Missing for the dam num =%s data MatriceAexp\nCarry out first calculation of geometric quantitiese !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    # coefficient matrix: ka;ma;kq;mq;kcel;mcel for each section
    dic_MatriceCoeff={}
    # list of section numbers
    ListaNumSez=[]
    for row in MatriceDati:
        tratto_cur=int(row[0])
        ListaNumSez.append(tratto_cur)
        dic_MatriceCoeff[tratto_cur]=row[3:]

    # Close communication with the database
    cur.close()
    conn.close()

    # reading of the sections
    # ---------------------

    nomecampoAltezza='hmax'

    driver = ogr.GetDriverByName('ESRI Shapefile')

    ds = driver.Open(CrossMedie, 1)
    if ds is None:
        errMsg = 'Could not open ' + CrossMedie
        NotErr= bool()
        return NotErr, errMsg

    layer = ds.GetLayer()

    feat = layer.GetNextFeature()

    Spatialref = layer.GetSpatialRef()
    Spatialref.AutoIdentifyEPSG()
    SourceEPSG=int(Spatialref.GetAuthorityCode(None))


    # list of points in left and right
    ListaPtSx=[]
    ListaPtDx=[]

    # dictionary of flood limit distances in left and right
    dic_DistSx={}
    dic_DistDx={}

    DV_sez={}
    Time_min_sez={}

    while feat:

        NumSez=feat.GetField('id')

        if NumSez==0:

            NumSez=ListaNumSez[0]
        # midpoint distance
        dist1 = feat.GetField('dist1')
        # progressive along the river path
        progr = feat.GetField('progr')

        linea=feat.GetGeometryRef()
        Length=linea.Length()

        pt1=linea.GetPoint(0)
        pt2=linea.GetPoint(1)

        Qsez=numpy.interp(progr, Progr_fiume_array, Qmax_array)
        Hsez=numpy.interp(progr, Progr_fiume_array, Hmax_array)
        Bsez=numpy.interp(progr, Progr_fiume_array, Bmax_array)
        Vsez=numpy.interp(progr, Progr_fiume_array, Vmax_array)
        Timesez=numpy.interp(progr, Progr_fiume_array, Time_array)

        # check if use energy elevation
        if UseEnergyHead:
            # instead of the depth of water use energy elevation
            hcinetica=Vsez**2/2.0/9.81
            Htot=Hsez+hcinetica
        else:
            Htot=Hsez

        # load the dictionary
        DV_sez[NumSez]=Qsez/Bsez
        Time_min_sez[NumSez]=int(Timesez/60.0)

        feat.SetField(nomecampoAltezza, Htot)

        layer.SetFeature(feat)

        # reading the widths of the wet area on the right and left
        # ..........................................................
        try:
            MatriceCoeff=dic_MatriceCoeff[NumSez]
        except:
            pass

        ka=float(MatriceCoeff[2])
        ma=float(MatriceCoeff[3])
        mb=ma-1.0
        # wet width for water level
        Bsez_tot=ka*ma*math.pow(Htot,mb)

        PercDx= numpy.interp(Htot,H_Array,MatriceArray[dic_PixDist[NumSez]])
        Bdx=Bsez_tot*PercDx
        Bsx=Bsez_tot-Bdx

        dic_DistSx[NumSez]=Bsx
        dic_DistDx[NumSez]=Bdx

        PercAscSx=(dist1-Bsx)/Length
        PercAscDx=(dist1+Bdx)/Length

        Pt_Sx= PuntoIntermedio(pt1,pt2,PercAscSx)
        Pt_Dx= PuntoIntermedio(pt1,pt2,PercAscDx)

        ListaPtSx.append(Pt_Sx)
        ListaPtDx.append(Pt_Dx)

        feat = layer.GetNextFeature()

    ds.Destroy()

    log_file=open('log.txt','a')
    log_file.write('End scrittura hmax\n')
    log_file.close()


    # making the polygon based on the river path
    # ......................................................


    try:
        # creating/connecting the db
        conn = db.connect(mydb_path_user)

    except:

        conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
       # import extention
        conn.enable_load_extension(True)
        conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    NomeTabellaLinee='Downstreampath'

    sql='SELECT TotalLength,ST_AsText(geom) FROM %s WHERE DamID=%d' % (NomeTabellaLinee,DamID)
    cur.execute(sql)
    ChkDiga=cur.fetchone()

    if ChkDiga==None:
        errMsg = "Nella tabella= %s non ci sono dati per la diga num =%s \nEffettuare prima il calcolo della linea a valle !" % (NomeTabellaLinee,DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt_line=ChkDiga[1]
        TotalLength=ChkDiga[0]
        StreamLine=ogr.CreateGeometryFromWkt(wkt_line)
        StreamLine.FlattenTo2D()

        dic_StreamTratti={}

        inDS1 = driver.Open(CrossSecPoly, 0)
        if inDS1 is None:
            errMsg = 'Could not open ' + CrossSecPoly
            NotErr= bool()
            return NotErr, errMsg

        InlayerCurve = inDS1.GetLayer()

        num_tratti = InlayerCurve.GetFeatureCount()

        feat = InlayerCurve.GetNextFeature()

        dic_NumTratto={}

        ii=-1
        while feat:

            NumSez=feat.GetField('id')

            ii+=1
            dic_NumTratto[ii]=NumSez

            poly=feat.GetGeometryRef()

            line_curr=poly.Intersection(StreamLine)
            if line_curr!=None:
                dic_StreamTratti[NumSez]=line_curr.ExportToWkt()
            else:
                txt='No intersection cross-sec num=%d' % NumSez
                print(txt)

            feat = InlayerCurve.GetNextFeature()

        inDS1.Destroy()

    # Close communication with the database
    cur.close()
    conn.close()

    ds = driver.Open(PolySxDx, 0)

    if ds is None:
        errMsg = 'Could not open ' + PolySxDx
        NotErr= bool()
        return NotErr, errMsg

    layer = ds.GetLayer()

    filtro="lato = %d" % 0
    layer.SetAttributeFilter(filtro)

    feat=layer.GetNextFeature()

    PoligonoSx=feat.GetGeometryRef()
    PoligonoSx_wkt=PoligonoSx.ExportToWkt()

    layer.SetAttributeFilter(None)
    layer.ResetReading()
    filtro="lato = %d" % 1
    layer.SetAttributeFilter(filtro)

    feat=layer.GetNextFeature()

    PoligonoDx=feat.GetGeometryRef()
    PoligonoDx_wkt=PoligonoDx.ExportToWkt()

    ds.Destroy()

    # initializing the polygon of the floodable area
    PoligonoAree1=ogr.Geometry(ogr.wkbPolygon)

    PolySx=ogr.CreateGeometryFromWkt(PoligonoSx_wkt)
    PolyDx=ogr.CreateGeometryFromWkt(PoligonoDx_wkt)

    dist_min_pixel=pixelWidth

    for i in range(num_tratti):
        ii=dic_NumTratto[i]
        linea_curr_wkt=dic_StreamTratti[ii]
        linea_curr=ogr.CreateGeometryFromWkt(linea_curr_wkt)
        for lato in range(2):
            # check left
            if lato==0:
                if dic_DistSx[ii]>dist_min_pixel:
                    polytratto=linea_curr.Buffer(dic_DistSx[ii])
                else:
                    polytratto=linea_curr.Buffer(dist_min_pixel)
                NewGeom=polytratto.Intersection(PolySx)
                if NewGeom!=None:
                    PoligonoAree1=PoligonoAree1.Union(NewGeom)
                    polytratto.Destroy()
                    NewGeom.Destroy()
                else:
                    PoligonoAree1=PoligonoAree1.Union(polytratto)
            # check right
            elif lato==1:
                if dic_DistDx[ii]>dist_min_pixel:
                    polytratto=linea_curr.Buffer(dic_DistDx[ii])
                else:
                    polytratto=linea_curr.Buffer(dist_min_pixel)
                NewGeom=polytratto.Intersection(PolyDx)
                if NewGeom!=None:
                    PoligonoAree1=PoligonoAree1.Union(NewGeom)
                    polytratto.Destroy()
                    NewGeom.Destroy()
                else:
                    PoligonoAree1=PoligonoAree1.Union(polytratto)

    log_file=open('log.txt','a')
    log_file.write('End PoligonoAree1\n')
    log_file.close()

    # making a shapefile with the PolygonAree1
    # ----------------------------------------
    shpnew_1=AreaInondabile_1
    if os.path.exists(shpnew_1):
        driver.DeleteDataSource(shpnew_1)

    outDS_1 = driver.CreateDataSource(shpnew_1)
    if outDS_1 is None:
        errMsg='Could not create file %s' % shpnew_1
        NotErr= bool()
        return NotErr, errMsg

    outLayer_1 = outDS_1.CreateLayer('AreaInondabile_1', Spatialref,geom_type=ogr.wkbMultiPolygon)

    fieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
    outLayer_1.CreateField(fieldDefn2)

    featureDefn_1 = outLayer_1.GetLayerDefn()

    feature = ogr.Feature(featureDefn_1)
    feature.SetField('id', 1)

    feature.SetGeometry(PoligonoAree1)
    outLayer_1.CreateFeature(feature)

    outDS_1.Destroy()



    # making the polygon # 2 based on the digital terrain model
    # ---------------------------------------------------------------------

    if not os.path.exists(StreamDH):
        errMsg = 'File StreamDHFilled %s does not exist' % os.path.realpath(StreamDH)
        NotErr= bool()
        return NotErr, errMsg


    infile=StreamDH

    indatasetElev = gdal.Open( infile, GA_ReadOnly )
    if indatasetElev is None:
        errMsg = 'Could not open ' + infile
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
    DH = inbandElev.ReadAsArray(0, 0, colsElev, rowsElev).astype(numpy.float32)

    mask_Nodata= DH==inNoDataElev

    inDS1 = driver.Open(CrossMedie, 0)
    if inDS1 is None:
        errMsg = 'Could not open ' + CrossMedie
        NotErr= bool()
        return NotErr, errMsg

    InlayerCurve = inDS1.GetLayer()

    spatialRef_sez=InlayerCurve.GetSpatialRef()

    feat_defn = InlayerCurve.GetLayerDefn()
    NumFields=feat_defn.GetFieldCount()

    # creates a grid with depth to cross sections
    GridSez=numpy.zeros((rowsElev,colsElev),numpy.float32)

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
    testo="ATTRIBUTE=%s"  % (nomecampoAltezza)
    # Rasterize
    outband = ds.GetRasterBand(iBand)

    outband.WriteArray(GridSez, 0, 0)
    CampoValore=[testo]

    err = gdal.RasterizeLayer(ds, [iBand], InlayerCurve,
            burn_values=[0],
            options=CampoValore)
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    # Reading WL
    GridSezWL = outband.ReadAsArray().astype(numpy.float32)

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
    xi = numpy.linspace(xmin, xmax, nx)
    yi = numpy.linspace(ymin, ymax, ny)
    xi, yi = numpy.meshgrid(xi, yi)

    # Reading x,y,z
    mask=GridSezWL>0
    x=xi[mask]
    y=yi[mask]
    z=GridSezWL[mask]


    # Otherwise, try Method 2 - Interpolate  using scipy interpolate griddata
    WLArray = il.griddata((x, y), z, (xi, yi),method='linear') #(may use 'nearest', 'linear' or 'cubic'  - although constant problems w linear)

    checkMask=numpy.isnan(WLArray)

    nnan=checkMask.sum()

    Nodata=-9999
    if nnan>0:
        WLArray= numpy.choose(checkMask,(WLArray,Nodata))

    # WaterDepth calculation by difference between water and ground level
    Wdepth=WLArray-DH

    # filtering of isolated points and internal empty points
    Wdepth= signal.medfilt2d(Wdepth,kernel_size=7)

    # eliminates negative values
    maskWd=Wdepth<=0.0
    Wdepth= numpy.choose(maskWd,(Wdepth,Nodata))

    # eliminate external anomalous values due to the filtering algorithm
    maskWd=Wdepth>9999
    Wdepth= numpy.choose(maskWd,(Wdepth,Nodata))

    # adds the nodata of the terrain model
    Wdepth= numpy.choose(mask_Nodata,(Wdepth,Nodata))


    # output file
    FileDEM_out=PathFiles+os.sep+'Hmax.tif'

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


    log_file=open('log.txt','a')
    log_file.write('End Hmax.tif\n')
    log_file.close()

    # ----------------------------
    # Rasterize PoligonoAree1
    # ------------------------
    PoligonoAree1_Raster=PathFiles+os.sep+'PoligonoAree1.tif'

    orig_data_source = ogr.Open(shpnew_1)
    source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_source, "")
    source_layer = source_ds.GetLayer()


    format = 'Gtiff'
    type = GDT_Int16

    driver3 = gdal.GetDriverByName(format)
    driver3.Register()


    dsRaster = driver3.Create(PoligonoAree1_Raster, cols, rows, 1, type)
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

    ClassTratti=numpy.zeros((rows,cols)).astype(numpy.int)

    outband.WriteArray(ClassTratti, 0, 0)

    # Rasterize
    err = gdal.RasterizeLayer(dsRaster, [1], source_layer,
            burn_values=[0],
            options=["ATTRIBUTE=id"])
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    # Reading from the raster of the matrix with value 1 in a flooded area
    MatriceDatiArea1=outband.ReadAsArray(0, 0, cols, rows)

    # eliminates any points with H greater than Hmax
    DH_MatriceDatiArea1= DH*MatriceDatiArea1
    mask_greatHmax =  DH_MatriceDatiArea1 > Hmax_tot
    nnn=mask_greatHmax.sum()
    MatriceDatiArea1=numpy.choose(mask_greatHmax,(MatriceDatiArea1,0))


    # writing Nodata
    mask_Nodata=MatriceDatiArea1==0
    MatriceDati=numpy.choose(mask_Nodata,(MatriceDatiArea1,outNodata))
    outband.WriteArray(MatriceDati, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(outNodata)
    outband.GetStatistics(0,1)

    outband=None

    dsRaster=None
    orig_data_source.Destroy()


    # ----------------------------

    # name of the output file with 1 in the wet cells
    FileDEM_out_1=PathFiles+os.sep+'HH.tif'

    format = 'GTiff'
    driver = gdal.GetDriverByName(format)
    type = GDT_Int16
    gt=indatasetElev.GetGeoTransform()

    ds = driver.Create(FileDEM_out_1, indatasetElev.RasterXSize, indatasetElev.RasterYSize, 1, type)
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
    WW=Wdepth>0
    # adding polygon areas 1
    # ---------------------------
    mask_Data1=MatriceDatiArea1==1
    # saving in the raster
    WW= numpy.choose(mask_Data1,(WW,1))

    outband.WriteArray(WW, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(Nodata)
    outband.GetStatistics(0,1)
    outband = None

    ds = None

    log_file=open('log.txt','a')
    log_file.write('End HH.tif\n')
    log_file.close()

    # Raster to vector
    # -------------------------

    # this allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    log_file=open('log.txt','a')
    log_file.write('End gdal.UseExceptions()\n')
    log_file.close()


    fileName=FileDEM_out_1
    src_ds = gdal.Open(fileName)
    if src_ds is None:
        errMsg = 'Could not open ' + fileName
        NotErr= bool()
        return NotErr, errMsg

    srcband = src_ds.GetRasterBand(1)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjection())

    log_file=open('log.txt','a')
    log_file.write('End srs.ImportFromWkt(src_ds.GetProjection()\n')
    log_file.close()

    dst_layername = "PolyFtr"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_filename=PathFiles+os.sep+dst_layername + ".shp"
    if os.path.exists(dst_filename):
        drv.DeleteDataSource(dst_filename)

    dst_ds = drv.CreateDataSource(dst_filename)
    dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs)
    newField = ogr.FieldDefn('id', ogr.OFTInteger)
    dst_layer.CreateField(newField)

    log_file=open('log.txt','a')
    log_file.write('End dst_layer.CreateField(newField)\n')
    log_file.close()

    # con bandmask
    gdal.Polygonize(srcband, srcband, dst_layer, 0, [],
    callback=None )

    log_file=open('log.txt','a')
    log_file.write('End Polygonize\n')
    log_file.close()

    src_ds=None
    dst_ds.Destroy()

    # deleting the temporary grid
    os.remove(fileName)

    log_file=open('log.txt','a')
    log_file.write('End remove HH.tif\n')
    log_file.close()

    # performing the union of the polygons
    # ----------------------------------
    in_layername=PathFiles+os.sep+ "PolyFtr.shp"

    shpdriver = ogr.GetDriverByName('ESRI Shapefile')

    inDS1 = shpdriver.Open(in_layername, 0)
    if inDS1 is None:
        errMsg = 'Could not open ' + in_layername
        NotErr= bool()
        return NotErr, errMsg

    InlayerCurve = inDS1.GetLayer()

    feat = InlayerCurve.GetNextFeature()

    poly_tot=ogr.Geometry(ogr.wkbMultiPolygon)

    while feat:

        poly=feat.GetGeometryRef()
        # aggiungo geometria poligonale
        poly_tot.AddGeometry(poly)

        feat = InlayerCurve.GetNextFeature()

    inDS1.Destroy()

    log_file=open('log.txt','a')
    log_file.write('End PolyFtr.shp\n')
    log_file.close()

    # creating the final flood area
    # -----------------------------

    # saving in the geodatabase
    # ---------------------------------------

    try:
        # creating/connecting the db
        conn = db.connect(mydb_path_user)

    except:

        conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
       # import extention
        conn.enable_load_extension(True)
        conn.execute('SELECT load_extension("mod_spatialite")')

    cur= conn.cursor()

    TargetTabella='FloodExtent'

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
        # delete previous data
        sql='DELETE FROM %s WHERE DamID=%d' % (TargetTabella,DamID)
        cur.execute(sql)
        conn.commit()


    inDS1 = shpdriver.Open(CrossSecPoly, 0)
    if inDS1 is None:
        errMsg = 'Could not open ' + CrossSecPoly
        NotErr= bool()
        return NotErr, errMsg

    InlayerCurve = inDS1.GetLayer()

    feat = InlayerCurve.GetNextFeature()

    while feat:

        NumSez=feat.GetField('id')

        poly=feat.GetGeometryRef()

        FloodSeverityString=FloodSeverity(DV_sez[NumSez])

        Factor=ConseqFactot(FloodSeverityString,Time_min_sez[NumSez])

        # making the intersection to get the polygon
        poly_curr=poly.Intersection(poly_tot)

        if poly_curr!=None:

            sql='INSERT INTO %s (DamID,id,DV,FloodSeverity,WarningTimeMin,FatalityRate,geom) VALUES (%d'  %  (TargetTabella,DamID)
            sql+=',%d' % NumSez
            sql+=',%.2f' % DV_sez[NumSez]
            sql+=',"%s"' % FloodSeverityString
            sql+=',%d' % Time_min_sez[NumSez]
            sql+=',%.3f' % Factor

            poly_curr.FlattenTo2D()

            # check if MULTIPOLYGON
            TipoGeom=poly_curr.GetGeometryName()
            if TipoGeom=='POLYGON':
                multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
                multipolygon.AddGeometry(poly_curr)
                wkt2=multipolygon.ExportToWkt()
            elif poly_curr.GetGeometryName()=='MULTIPOLYGON':
                wkt2=poly_curr.ExportToWkt()
            poly2=ogr.CreateGeometryFromWkt(wkt2)

            GeomWKT="GeomFromText('%s',%d)" % (wkt2,OriginEPSG)
            sql+=',%s' % GeomWKT
            sql+=');'
            cur.execute(sql)

        else:

            log_file=open('log.txt','a')
            log_file.write('Err tratto n=%d\n' % NumSez)
            log_file.close()


        feat = InlayerCurve.GetNextFeature()


    inDS1.Destroy()

    log_file=open('log.txt','a')
    log_file.write('End routine\n')
    log_file.close()

    conn.commit()
    # Close communication with the database
    cur.close()
    conn.close()


    return NotErr, errMsg

if __name__ == '__main__':

    myGDB_user='..' + os.sep + 'db'+os.sep+ 'USER_GeoDB.sqlite'

    # ==================
    # San Giuliano  dam
    # ==================
    DamID=449

    root_dir = os.path.dirname(__file__) + os.sep + '..'

    PathFiles=root_dir +os.sep+ str(DamID)
    PathFiles=os.path.realpath(PathFiles)

    UseEnergyHead= bool('True')
    NotErr, errMsg = GridFloodingAreas(myGDB_user,PathFiles,DamID,UseEnergyHead,UseEnergyHead)

    print(NotErr,errMsg)
