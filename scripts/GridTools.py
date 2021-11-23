# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SimpleDamBrk - GridTolls
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

"""
import os, sys
import csv

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
except:
    import ogr

try:
    from osgeo import osr
except:
    import osr

import numpy as np
import math
import time
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from matplotlib.font_manager import FontProperties

def ControlloCongruenzaGRID(OriginData,indataset,gt):
    toll=0.5
    toll2=0.001
    ok=bool('True')
    # OriginData=[originX,originY,pixelWidth,pixelHeight,cols,rows]
    if abs(OriginData[0] - gt[0])>toll:
        print ('originX non congruente con OriginData')
        ok=bool()
        #sys.exit(1)
    if abs(OriginData[1] - gt[3])>toll:
        print ('originY non congruente con OriginData')
        ok=bool()
        #sys.exit(1)
    if abs(OriginData[2] - gt[1])>toll2:
        print ('pixelWidth non congruente con OriginData')
        ok=bool()
        #sys.exit(1)
    if abs(OriginData[3] - gt[5])>toll2:
        print ('pixelHeight non congruente con OriginData')
        ok=bool()
        #sys.exit(1)
    if OriginData[4] != indataset.RasterXSize:
        print ('cols non congruente con OriginData')
        ok=bool()
        #sys.exit(1)
    if OriginData[5] !=indataset.RasterYSize:
        print ('rows non congruente con OriginData')
        ok=bool()
        #sys.exit(1)
    return ok


def GridDistanzaFiume(mydb_path_user,DamID,PathFiles,ClipDEM,DistanzaMax=4000.0):

    """
    The script load:
    - grid ClipDEM
    - shapefile total cross-sec
    - shapefile river-path
    create:
    - shapefile (in memory) of polygons at different distances from the line of the river's path
    - grid Distances.tif with in values of distances with respect to the river

    """

    NotErr=bool('True')
    errMsg='OK'

    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "Missing data for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    CrossSec=PathFiles+os.sep+'CrossSec.shp'
    CrossSecTot=PathFiles+os.sep+'CrossSecTot.shp'

    if not os.path.exists(CrossSecTot):
        errMsg = "Missing data for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    if not os.path.exists(ClipDEM):
        errMsg = "Missing for the dam num =%s il clip DEM\nFirst make the clip of the digital terrain model !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    CrossSecPoly='%sPoly.shp' % (CrossSec[:-4])
    if not os.path.exists(CrossSecPoly):
        errMsg = "Missing for the dam num =%s il Grid Tratti\nFirst made the CreaSezInterpolate !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    # =====================================
    # Opening user database sqlite
    # ====================================

    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
   # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    # check downstream line
    NomeTabellaLinee='Downstreampath'

    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
    cur.execute(sql)

    record=cur.fetchone()
    if record!=None:
        DominioEPSG=record[0]
    else:
        DominioEPSG=32632

    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(DominioEPSG)


    sql='SELECT TotalLength,ST_AsText(geom) FROM %s WHERE DamID=%d' % (NomeTabellaLinee,DamID)
    cur.execute(sql)
    ChkDiga=cur.fetchone()

    if ChkDiga==None:
        errMsg = "Into the table= %s there is no data for the dam num =%s \nPerform the downstream line calculation first !" % (NomeTabellaLinee,DamID)
        NotErr= bool()
        return NotErr, errMsg

    else:
        wkt_line=ChkDiga[1]
        TotalLength=ChkDiga[0]


    # Close communication with the database
    cur.close()
    conn.close()

    # output file name
    FileDEM_out=PathFiles+os.sep+'Distances.tif'

    # ==================
    # Legge i dati GRID
    # ==================

    gdal.AllRegister()

    indataset = gdal.Open(ClipDEM, GA_ReadOnly )
    if indataset is None:
        print ('Could not open ' + ClipDEM)
        sys.exit(1)

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
        spatialRef.ImportFromEPSG(32632)

    terreno = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)
    mask_Nodata= terreno==inNoData

    gt=indataset.GetGeoTransform()

    inband = None

    indataset = None

    driver = ogr.GetDriverByName('ESRI Shapefile')


    InDS = driver.Open(CrossSecTot, 0)
    if InDS is None:
    	print ('Could not open ' + CrossSec)
    	sys.exit(1)     #exit with an error code


    Inlayer = InDS.GetLayer()

    inFeature=Inlayer.GetNextFeature()

    Lengthtmax=0.0

    while inFeature:
        NumSez=inFeature.GetField('id')
        geom = inFeature.GetGeometryRef()
        length=geom.Length()
        if length>Lengthtmax:
            Lengthtmax=length
        inFeature=Inlayer.GetNextFeature()

    InDS.Destroy()

    # number of steps
    NumSteps=int(Lengthtmax/pixelWidth/1.5)

    ListaDist=[]
    for i in range(NumSteps):
        dist=(i+1)*pixelWidth
        ListaDist.append(dist)

    shpnew2=PathFiles+os.sep +"distances.shp"
    nomecampoDist='dist'


    PathGeom=ogr.CreateGeometryFromWkt(wkt_line)


    #create an output datasource in memory
    outdriver=ogr.GetDriverByName('MEMORY')
    outDS2=outdriver.CreateDataSource('memData')

    #open the memory datasource with write access
    tmp=outdriver.Open('memData',1)

    outLayer2 = outDS2.CreateLayer('distances', dest_srs,geom_type=ogr.wkbPolygon)

    fieldDefn1 = ogr.FieldDefn('id', ogr.OFTInteger)
    outLayer2.CreateField(fieldDefn1)
    fieldDefn2 = ogr.FieldDefn(nomecampoDist, ogr.OFTReal)
    outLayer2.CreateField(fieldDefn2)

    featureDefn2 = outLayer2.GetLayerDefn()

    bufferDistance = ListaDist[0]

    NewGeom = PathGeom.Buffer(bufferDistance)

    poly_old_wkt=NewGeom.ExportToWkt()

    feature = ogr.Feature(featureDefn2)
    feature.SetField('id', 0)
    feature.SetField(nomecampoDist,ListaDist[0])
    feature.SetGeometry(NewGeom)
    outLayer2.CreateFeature(feature)
    NewGeom.Destroy()

    for i in range(1,NumSteps):
        bufferDistance=ListaDist[i]
        NewGeom1 = PathGeom.Buffer(bufferDistance)

        poly_old=ogr.CreateGeometryFromWkt(poly_old_wkt)

        poly_old_wkt=NewGeom1.ExportToWkt()

        # symmetric difference
        simdiff = NewGeom1.SymmetricDifference(poly_old)
        area=simdiff.Area()

        feature = ogr.Feature(featureDefn2)
        feature.SetField('id', i)
        feature.SetField(nomecampoDist,ListaDist[i])
        feature.SetGeometry(simdiff)
        outLayer2.CreateFeature(feature)
        simdiff.Destroy()

    Inlayer =outLayer2


    # Rasterize
    # -------------------
    format = 'GTiff'
    type = GDT_Float32

    driver2 = gdal.GetDriverByName(format)
    driver2.Register()

    DistMax=np.zeros((rows,cols),np.float32)
    DistMax=DistMax+ListaDist[-1]

    ds = driver2.Create(FileDEM_out, cols,  rows, 1, type)
    if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
        ds.SetGeoTransform(gt)

    if prj is not None and len(prj) > 0:
        ds.SetProjection(prj)
    else:
        prj= spatialRef.ExportToWkt()
        ds.SetProjection(prj)

    iBand=1
    testo="ATTRIBUTE=%s"  % (nomecampoDist)
    # Rasterize
    outband = ds.GetRasterBand(iBand)

    # Rasterize
    outband.WriteArray(DistMax, 0, 0)
    CampoValore=[testo]

    # create the map of values
    # -------------------------------------------------
    err = gdal.RasterizeLayer(ds, [iBand], Inlayer,
            burn_values=[0],
            options=CampoValore)
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    # Reading
    dist = outband.ReadAsArray().astype(np.float32)
    Nodata=-9999
    dist= np.choose(mask_Nodata,(dist,Nodata))
    # save
    outband.WriteArray(dist, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(Nodata)
    outband.GetStatistics(0,1)
    outband = None

    ds= None

    # creation of the grid of river sections
    # ---------------------------------------
    CreaCrid=1

    if CreaCrid>0:

        mask_Dist= np.greater_equal(dist,DistanzaMax)


        # ==================
        # Reading GRID
        # ==================

        gdal.AllRegister()

        indataset = gdal.Open(ClipDEM, GA_ReadOnly )
        if indataset is None:
            errMsg='Could not open file %s' % ClipDEM
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
            # default
            # WGS84 UTM 32 N
            spatialRef.ImportFromEPSG(32632)

        terreno = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)
        mask_Nodata= terreno==inNoData


        # Rasterize polygon
        # ====================
        orig_data_source = ogr.Open(CrossSecPoly)
        source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_source, "")
        source_layer = source_ds.GetLayer()

        result=source_ds.ExecuteSQL('SELECT id,zmin FROM CrossSecPoly')
        resultFeat = result.GetNextFeature()
        ListaId=[]
        QuotaId={}
        while resultFeat:
            dd=resultFeat.GetField('id')
            ListaId.append(dd)
            zz=resultFeat.GetField('zmin')
            QuotaId[dd]=zz
            resultFeat = result.GetNextFeature()
        source_ds.ReleaseResultSet(result)

        format = 'Gtiff'
        type = GDT_Int16

        driver3 = gdal.GetDriverByName(format)
        driver3.Register()

        TrattiRaster=PathFiles+os.sep+'Tratti.tif'

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

        # writes -1 to the whole matrix
        ClassTratti=np.zeros((rows,cols)).astype(np.int)
        ClassTratti=ClassTratti-1

        outband.WriteArray(ClassTratti, 0, 0)

        # Rasterize
        err = gdal.RasterizeLayer(dsRaster, [1], source_layer,
                burn_values=[0],
                options=["ATTRIBUTE=id"])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)

        MatriceDati=outband.ReadAsArray(0, 0, cols, rows)
        MatriceDati=np.choose(mask_Nodata,(MatriceDati,outNodata))

        # discard points beyond the maximum distance by assigning Nodata
        MatriceDati=np.choose(mask_Dist,(MatriceDati,outNodata))

        outband.WriteArray(MatriceDati, 0, 0)

        outband.FlushCache()
        outband.SetNoDataValue(outNodata)
        outband.GetStatistics(0,1)

        outband=None

        dsRaster=None
        orig_data_source.Destroy()


    return NotErr, errMsg

def ModDTM(mydb_path_user,DamID,PathFiles,ClipDEM):

    """
    The script load:
        - StreamDH grid
        - Distances grid
        - shapefile river-path
    create:
        - StreamDHFilled grid:
          heights that we have at least one slope with respect to the river> = PendMin
          and that on the course of the river they have value = 0

    """

    NotErr=bool('True')
    errMsg='OK'

    # minimum slope to the river
    PendMin=float(0.002)

    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "Missing data for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDH.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Missing for the dam num =%s il StreamDH\nPerform first CreaSezInterpolate !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    Distances=PathFiles+os.sep+'Distances.tif'
    if not os.path.exists(Distances):
        errMsg = "Missing for the dam num =%s il Distances\nPerform first CreaGridDistanza !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    # import extention
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')

    # creating a Cursor
    cur = conn.cursor()

    # check downstream line
    NomeTabellaLinee='Downstreampath'

    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
    cur.execute(sql)

    record=cur.fetchone()
    if record!=None:
        DominioEPSG=record[0]
    else:
        DominioEPSG=32632

    dest_srs = osr.SpatialReference()
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
        TotalLength=ChkDiga[0]


    # Close communication with the database
    cur.close()
    conn.close()

    # output file name
    FileDEM_out=PathFiles+os.sep+'StreamDHFilled.tif'

    # ==================================
    # Reading GRID data of STREAM DH
    # ==================================

    gdal.AllRegister()

    indataset = gdal.Open(StreamDH, GA_ReadOnly )
    if indataset is None:
        errMsg = "Could not open " + StreamDH
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

    StreamDHA = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)

    # creates the mask of the area with Nodata
    mask_Nodata= StreamDHA==inNoData

    gt=indataset.GetGeoTransform()

    inband = None

    indataset = None

    # ==================================
    # Reading GRID data of DISTANCES
    # ==================================

    gdal.AllRegister()

    indataset = gdal.Open(Distances, GA_ReadOnly )
    if indataset is None:
        print ('Could not open ' + Distances)
        sys.exit(1)

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

    DistancesA = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)

    gt=indataset.GetGeoTransform()

    inband = None

    indataset = None

    # PT 2 --------------------

    MatriceHmin = DistancesA*PendMin

    # PT 3 e 4 ----------------

    # search for a mask of the StreamDH input points that do not satisfy
    # he mask = StreamDH <PendMin condition (operation between an array of numpy and the value PendMin)

    Inf=StreamDHA<MatriceHmin

    StreamDHFilled=np.choose(Inf,(StreamDHA,MatriceHmin))

    StreamDHFilled=np.choose(mask_Nodata,(StreamDHFilled,inNoData))

    # PT 5 - creation of the river axis grid by reading the shapefile of the river axis

    #create an output datasource in memory
    outdriver=ogr.GetDriverByName('MEMORY')
    source=outdriver.CreateDataSource('memData')

    #open the memory datasource with write access
    tmp=outdriver.Open('memData',1)

    outLayer = source.CreateLayer("stream", dest_srs, geom_type=ogr.wkbLineString)
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    outLayer.CreateField(idField)

    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    line=ogr.CreateGeometryFromWkt(wkt_line)
    feature.SetGeometry(line)
    feature.SetField("id", 1)
    outLayer.CreateFeature(feature)
    feature = None


    # Rasterize
    # -------------------
    format = 'MEM'
    type = GDT_Int16

    driver2 = gdal.GetDriverByName(format)
    driver2.Register()

    FileDEM_tmp=PathFiles+os.sep+'Stream.tif'

    ds = driver2.Create(FileDEM_tmp, cols,  rows, 1, type)
    if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
        ds.SetGeoTransform(gt)

    if prj is not None and len(prj) > 0:
        ds.SetProjection(prj)
    else:
        prj= spatialRef.ExportToWkt()
        ds.SetProjection(prj)

    iBand=1
    # Rasterize
    outband = ds.GetRasterBand(iBand)

    # burn_values=[1] indicates to rasterize the vector with a value of 1
    err = gdal.RasterizeLayer(ds, [iBand], outLayer,
            burn_values=[1])
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    # Reading
    stream = outband.ReadAsArray().astype(np.int)
    Nodata=-9999
    stream= np.choose(mask_Nodata,(stream,Nodata))
    outband.WriteArray(stream, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(Nodata)
    outband.GetStatistics(0,1)
    outband = None

    ds= None


    # PT 6 ------------- replacement in StreamDHFilled of the value 0 at the points where stream = 1

    maskstream=stream==1

    StreamDHFilled=np.choose(maskstream,(StreamDHFilled,0.0))

    # PT 7 ------------- saving the final StreamDHFilled matrix in the file tif = StreamDHFilled.tif

    # output file name
    FileDEM_out=PathFiles+os.sep+'StreamDHFilled.tif'

    format = 'GTiff'
    driver = gdal.GetDriverByName(format)
    type = GDT_Float32

    driver3 = gdal.GetDriverByName(format)
    driver3.Register()

    # writing raster
    ds = driver3.Create(FileDEM_out, cols,  rows, 1, type)
    if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
        ds.SetGeoTransform(gt)

    if prj is not None and len(prj) > 0:
        ds.SetProjection(prj)
    else:
        prj= spatialRef.ExportToWkt()
        ds.SetProjection(prj)

    iBand=1
    # Rasterize
    outband = ds.GetRasterBand(iBand)
    outband.WriteArray(StreamDHFilled, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(Nodata)
    outband.GetStatistics(0,1)
    outband = None

    ds = None

    return NotErr, errMsg

def CurveAreaAltezza(mydb_path_user,DamID,PathFiles):

    """
    The cript loads:
        - StreamDHFilled.tif
        - Tratti.tif         : grid with classes of stretches of river starting from the dam
                               the number represents in number of cells of
                               distance counted along the axis of the river

    counts for each section and for each dh with a 1-meter step in the number
    of cells underlying the difference in height dh from the river
    It also constructs the curves of the number of cumulated cells along the path of the
    river for every dh with 1 meter step

    Save counts in two csv files:
        - MatricePixel.csv
        - MatricePixCum.csv
    """


    NotErr=bool('True')
    errMsg='OK'

    PathFiles=os.path.realpath(PathFiles)

    # files output
    filecsv1=PathFiles+os.sep+'MatricePixel.csv'
    filecsv2=PathFiles+os.sep+'MatricePixCum.csv'

    # contains the percentage of area in the orographic right
    filecsvPixDestra=PathFiles+os.sep+'MatricePixDestra.csv'

    if not os.path.exists(PathFiles):
        errMsg = "Missing data for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDHFilled.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Missing for the dam num =%s the grid StreamDHFilled\nPerform first ModificaDH !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    Tratti=PathFiles+os.sep+'Tratti.tif'
    if not os.path.exists(Tratti):
        errMsg = "Missing for the dam num =%s the Grid Tratti\nPerform first CreaSezInterpolate !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    DestraSinistra=PathFiles+os.sep+'DestraSinistra.tif'
    if not os.path.exists(DestraSinistra):
        errMsg = "Missing for the dam num =%s the Grid DestraSinistra\nPerform first CreaSezInterpolate !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg



    # ==================
    # Reading GRID
    # ==================

    gdal.AllRegister()

    indataset = gdal.Open(Tratti, GA_ReadOnly )
    if indataset is None:
        print ('Could not open ' + Tratti)
        sys.exit(1)

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

    cellsize=pixelWidth

    OriginData=[originX,originY,pixelWidth,pixelHeight,cols,rows]

    prj = indataset.GetProjectionRef()

    spatialRef = osr.SpatialReference()
    try:
        spatialRef.ImportFromWkt(prj)
    except:
        pass
        spatialRef.ImportFromEPSG(32632)

    TrattiArray = inband.ReadAsArray(0, 0, cols, rows).astype(np.int)

    TrattiVector_ini=np.unique(TrattiArray)

    inband = None
    indataset = None

    # Reading StreamDH
    # -----------------

    if not os.path.exists(StreamDH):
        errMsg = "File StreamDH %s does not exists" % os.path.realpath(StreamDH)
        NotErr= bool()
        return NotErr, errMsg

    infile=StreamDH

    indatasetElev = gdal.Open( infile, GA_ReadOnly )
    if indatasetElev is None:
        print ('Could not open ' + infile)
        sys.exit(1)

    prj = indatasetElev.GetProjectionRef()

    gt = indatasetElev.GetGeoTransform()

    ok=ControlloCongruenzaGRID(OriginData,indatasetElev,gt)

    if not ok:
        errMsg = 'Grid error: %s' % infile
        NotErr= bool()
        return NotErr, errMsg


    originXElev = gt[0]
    originYElev = gt[3]
    pixelWidthElev = gt[1]
    pixelHeightElev = gt[5]
    colsElev=indatasetElev.RasterXSize
    rowsElev=indatasetElev.RasterYSize
    bandsElev=indatasetElev.RasterCount
    iBand = 1
    inbandElev = indatasetElev.GetRasterBand(iBand)
    inNoDataElev= inbandElev.GetNoDataValue()

    # reading the entire file at once
    DH = inbandElev.ReadAsArray(0, 0, colsElev, rowsElev).astype(np.float32)
    mask_Nodata= DH==inNoDataElev

    TrattiVector_ini=np.unique(TrattiArray)

    # assigning, for congruency, the Nodata map of the DH also to the TrattiArray
    TrattiArray=np.choose(mask_Nodata,(TrattiArray,inNoData))

    TrattiVector_ini2=np.unique(TrattiArray)

    # eliminates negative values
    DH=np.choose(np.less(DH,0.0),(DH,0.0))

    inbandElev=None

    indatasetElev= None

    # Reading grid DestraSinistra
    # --------------------------
    indatasetDx = gdal.Open(DestraSinistra, GA_ReadOnly )
    if indatasetDx is None:
        errMsg = 'Could not open ' + DestraSinistra
        NotErr= bool()
        return NotErr, errMsg

    gtDx = indatasetDx.GetGeoTransform()

    ok=ControlloCongruenzaGRID(OriginData,indatasetDx,gtDx)

    if not ok:
        errMsg = 'Grid non conguente: %s' % indatasetDx
        NotErr= bool()
        return NotErr, errMsg

    originX = gtDx[0]
    originY = gtDx[3]
    pixelWidth = gtDx[1]
    pixelHeight = gtDx[5]
    cols=indatasetDx.RasterXSize
    rows=indatasetDx.RasterYSize
    bands=indatasetDx.RasterCount
    iBand = 1
    inbandDx = indatasetDx.GetRasterBand(iBand)
    inNoDataDx= inbandDx.GetNoDataValue()

    prj = indatasetDx.GetProjectionRef()

    spatialRef = osr.SpatialReference()
    try:
        spatialRef.ImportFromWkt(prj)
    except:
        pass
        spatialRef.ImportFromEPSG(32632)

    # reading array of the right fluvial zone = 1
    TrattiArrayDx = inbandDx.ReadAsArray(0, 0, cols, rows).astype(np.int)

    # creating the mask of the right area
    mask_Dx=np.equal(TrattiArrayDx,1)
    numDx=mask_Dx.sum()

    inbandDx = None
    indatasetDx = None


    # creates the list of river sections
    # ------------------------
    TrattiVector1=np.unique(TrattiArray)
    # discard data <0 (Nodata and data outside the river sections)
    mask=np.where(TrattiVector1>0)[0]

    # final array
    TrattiVector=TrattiVector1[mask]

    # Maximum height for which the curves are created
    Hmax=51

    MatricePix=[]
    VettorePixHmax=[]
    VettoreVolHmax=[]

    # Matrix of the percentage on the right
    MatricePixDx=[]

    for tratto in TrattiVector:
        # mask
        mask_tratto=np.equal(TrattiArray,tratto)

        numeropixel=mask_tratto.sum()
        VettorePix=[]
        VettorePixDx=[]
        for h in range(1,Hmax):
            mask=np.less_equal(DH,h) & mask_tratto
            nn=mask.sum()
            if nn>0:
                # select the pixels of the river section <h
                DH_cur=DH[np.where(mask)]
                # sum the heights: equivalent to the volume in terms of pixels * h
                Vol_h=np.sum(np.absolute(DH_cur),dtype=np.float32)
                # the empty volume is obtained by the difference
                Vol_d_valle=float(nn)*h-Vol_h

                VettorePix.append(Vol_d_valle)

                # finding those on the right hydrographic
                # ---------------------------------
                mask2=mask & mask_Dx
                ndx=mask2.sum()
                # select the pixels of the river section <h
                DH_curDx=DH[np.where(mask2)]
                # sum the heights: equivalent to the volume in terms of pixels * h
                Vol_h_Dx=np.sum(np.absolute(DH_curDx),dtype=np.float32)
                # the empty volume is obtained by the difference
                Vol_h_valleDx=float(ndx)*h-Vol_h_Dx

                if Vol_d_valle>0:
                    PercDx=Vol_h_valleDx/Vol_d_valle
                else:
                    PercDx=0.0
                VettorePixDx.append(PercDx)
            else:
                VettorePix.append(0.0)
                VettorePixDx.append(0.0)

        MatricePix.append(VettorePix)
        MatricePixDx.append(VettorePixDx)

    # creates the matrix of the areas accumulated at equal height from the riverbed
    MatricePixCum=[]
    nh=len(MatricePix[0])
    ndist=len(MatricePix)
    for j in range(nh):
        curva=[]
        Acum=0
        for i in range(ndist):
            Acum+=MatricePix[i][j]
            curva.append(Acum)
        MatricePixCum.append(curva)


    # save MatricePixel & MatricePixDx matrices
    fout=open(filecsv1,'w')
    foutDx=open(filecsvPixDestra,'w')

    txt='PixDist'
    hvals=[]
    for j in range(1,Hmax):
        nome='h=%.1f' % (j)
        hvals.append(nome)
        txt+=';h=%.1f' % (j)
    txt+='\n'
    fout.write(txt)
    foutDx.write(txt)

    nn=len(TrattiVector)
    for i in range(nn):
        txt='%d' % TrattiVector[i]
        row = MatricePix[i]
        for rec in row:
            txt+=';%.2f' % rec
        txt+='\n'
        fout.write(txt)

        rowDx = MatricePixDx[i]
        txt='%d' % TrattiVector[i]
        for rec in rowDx:
            txt+=';%.4f' % rec
        txt+='\n'
        foutDx.write(txt)

    fout.close()
    foutDx.close()

    # save MatricePixCum
    fout=open(filecsv2,'w')
    txt='PixDist'
    for j in range(1,Hmax):
        txt+=';h=%.1f' % (j)
    txt+='\n'
    fout.write(txt)
    nn=len(TrattiVector)
    for i in range(nn):
        txt='%d' % TrattiVector[i]
        for j in range(nh):
           txt+=';%.2f' % MatricePixCum[j][i]
        txt+='\n'
        fout.write(txt)
    fout.close()



    grafici=0
    if grafici>0:
        # ----------
        fontP = FontProperties()
        fontP.set_size('small')

        fig, ax = plt.subplots()
        nomi=[]
        x=np.array(TrattiVector,dtype =np.float)
        x=x*cellsize/1000.0
        A1=np.array(MatricePix,dtype =np.float)
        A1T=A1.T

        i=-1
        for curva in A1T:
            i=i+1
            ax.plot(x, curva, '-o')
            nomi.append(hvals[i])

        ax.legend( (nomi),loc='upper left',prop = fontP )
        ax.set_ylabel('Num. Cells',color='blue')
        ax.set_xlabel('Distance (km)',color='blue')
        txt='Area of the Valley for different heights from the river, starting from the dam\n'
        plt.title(txt)

        plt.grid()
        plt.show()

        # ---------
        # GRAPH 2
        fig, ax = plt.subplots()
        nomi=[]
        A=np.array(MatricePixCum,dtype =np.float)

        i=-1
        for curva in A:
            i=i+1
            ax.plot(x, curva, '-')
            nomi.append(hvals[i])

        ax.legend( (nomi),loc='upper left',prop = fontP )
        ax.set_ylabel('Num. Cells',color='blue')
        ax.set_xlabel('Distance (km)',color='blue')
        txt='Cumulated Valley Area for different heights from the river, starting from the dam\n'
        plt.title(txt)

        plt.grid()
        plt.show()

    return NotErr, errMsg

def ParamGeomHydro(mydb_path_user,DamID,PathFiles):
    """
    The script calculates the curve of the conveyance for a cross section representative of the river's reach
    Loads distances and differences in height from the point file along the axis of the river from
     - PathPoints.shp

    Inputs data:
     - cellsize : cell size of the starting raster
     - nMann    : n of Manning assumed constant
     - pendMin  : minimum acceptable slope, in case of lower slopes it assumes pendMin

    Loads the curves of the river reach areas from the file
     - MatricePixel.csv

    Calculate the distances and slopes of each river reach from the shapefile PathPoints.shp
    Starting from the points of the reach areas, find a monomy formula using the least squares
    to approximate the area variation with water depth

    Area=ka*h^ma

    from this, applying the uniform motion formula calculates a monomy formula
    of the rating curves, assimilating the hydraulic radius to the height h

    Q= 1/nMann*ka**radq(pend)*h^(ma+2/3) = kq*h^mq

    Calculate the celerity c = dQ/dA always as a monomial formula

    c=kcel*h^mcel

    Calculate the width of the free surface B deriving dA / dh

     B=ka*ma*h^(ma-1) = kB*h^mB

    Save the coefficients and exponents found of the monomule formulas in the file:

        - MatriceAexp.csv

    """

    pythonver_log=sys.version
    ppp=pythonver_log.split(' ')
    py_ver=ppp[0][0:1]


    NotErr=bool('True')
    errMsg='OK'

    PathFiles=os.path.realpath(PathFiles)

    # An acceptable minimum slope of one per thousand is defined
    pendMin=0.001
    # A constant Manning value is defined
    nMann=0.06

    if not os.path.exists(PathFiles):
        errMsg = "Missing data for the dam num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDHFilled.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Missing for the dam  num =%s il StreamDH\nPerform first ModificaDH !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg


    filecsv1=PathFiles+os.sep+'MatricePixel.csv'
    if not os.path.exists(filecsv1):
        errMsg = "Missing for the dam  num =%s il MatricePixel.csv\nPerform first CreaCurveAreaAltezza !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    # =====================================
    # Connecting the database sqlite
    # ====================================

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

    TabellaPoints='PathPoints'

    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (TabellaPoints.lower())
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

        errMsg = "Missing data for the dam  num =%s \nPerform the calculation of the downstream sections first !" % (DamID)
        NotErr= bool()
        return NotErr, errMsg

    NumSezTipo={}
    NumSezElev={}
    NumSezProgr={}
    NumSezGeom={}

    for rec in ListaTratti:

        NumSez=rec[0]
        tipo=rec[1]
        NumSezTipo[NumSez]=rec[1]
        NumSezElev[NumSez]=rec[2]
        NumSezProgr[NumSez]=rec[3]
        NumSezGeom[NumSez]=rec[4]


    # ==================================
    # Find cellsize
    # ==================================

    # Reading from STREAM DH
    gdal.AllRegister()

    indataset = gdal.Open(StreamDH, GA_ReadOnly )
    if indataset is None:
        print ('Could not open ' + StreamDH)
        sys.exit(1)

    geotransform = indataset.GetGeoTransform()

    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]

    cellsize=pixelWidth

    indataset = None


    # Number of Feature
    n = len(ListaTratti)

    id0=ListaTratti[0][0]
    geom=ogr.CreateGeometryFromWkt(NumSezGeom[id0])

    # reading the coordinates
    x0=geom.GetX()
    y0=geom.GetY()
    # reading river reach
    elev0=NumSezElev[id0]
    ProgrMonte=NumSezProgr[id0]

    # dictionary of the reach slope
    TrattoPend={}
    # dictionary of reach length in a straight line
    TrattoDistanza={}
    # dictionary of the length of the reach along the river
    TrattoDistanza_fiume={}
    # dictionary of the average value of the progressive reach along the rive
    TrattoProg_fiume={}

    numtratti=len(ListaTratti)

    for itratto in range(1,numtratti):

        # reading data
        id1=ListaTratti[itratto][0]
        elev1=NumSezElev[id1]
        ProgrValle=NumSezProgr[id1]


        # reading his geometry
        geom = ogr.CreateGeometryFromWkt(NumSezGeom[id1])

        x1=geom.GetX()
        y1=geom.GetY()
        # find the distance in a straight line between the two points
        distanza=math.sqrt(math.pow((x1-x0),2)+math.pow((y1-y0),2))
        # find the partial distance, as the difference of the progressive ones
        distanza_fiume=ProgrValle-ProgrMonte
        # find the intermediate progressive distance
        ProgrMedia=(ProgrMonte+ProgrValle)/2.0

        # the slope is the tangent of the angle
        pend=(elev0-elev1)/distanza

        # check if the slope <0
        if pend<=pendMin:
            pend=pendMin*1.0

        # save in dictionaries
        TrattoPend[id1]=pend
        TrattoDistanza[id1]=distanza
        TrattoDistanza_fiume[id1]=distanza_fiume
        TrattoProg_fiume[id1]=ProgrMedia

        # go to the next step
        # point 1 of a reach becomes the 0 point of the next reach
        x0=x1*1.0
        y0=y1*1.0
        elev0=elev1*1.0
        ProgrMonte=ProgrValle*1.0


    # reading Curves of reach volumes
    fin=open(filecsv1,'r')
    reader = csv.reader(fin,delimiter=';')
    if py_ver=='2':
        headers = reader.next()
    else:
        headers = reader.__next__()


    nn=len(headers)
    hvals=[]
    Vettoreh=[]
    for i in range(1,nn):
        hvals.append(headers[i])
        pp=headers[i].split('=')
        Vettoreh.append(float(pp[1]))

    # water depth array
    H_Array=np.array(Vettoreh,dtype =np.float)

    ListaAscisse=[]
    MatricePix=[]
    for row in reader:
        ListaAscisse.append(row[0])
        MatricePix.append(row[1:])

    fin.close()

    # matrix of the areas intended as numeropixel * h for each dh with respect to the river
    MatriceArray=np.array(MatricePix,dtype =np.float)

    # cell area
    cellArea=cellsize*cellsize

    # number of cells affecting the length of the first reach
    # this is also the constant step with which they are also built
    # the subsequent reachs
    deltacelle=int(ListaAscisse[0])

    # find areas where the area is less than the minimum:
    # it assumes at least a number of cells equal to those of the path
    # of the river's path in a reach
    zonamin=MatriceArray<deltacelle

    # clear areamin areas
    MatriceArray=np.choose(zonamin,(MatriceArray,deltacelle))


    # The volume in the Array Matrix is the number of cells x height
    # the actual volume in square meters is obtained by multiplying by the area of the cell
    MatriceAree=MatriceArray*cellArea


    # ------------------------------------------------------------------------------------------
    # calculates for each reach the monomic curve of the Area according to the depth of the water
    # -------------------------------------------------------------------------------------------

    # save the MatriceAexp

    # delete any previous data
    # -----------------------------------
    NomeTabella='MatriceAexp'
    sql='DELETE FROM %s WHERE DamID=%d' % (NomeTabella,DamID)
    cur.execute(sql)
    conn.commit()

    nn=len(MatriceAree)

    # Area chart
    graficoA=0
    # wet area width chart
    graficoB=0


    for i in range(nn):
        curva1= MatriceAree[i]

        # eliminates extreme points
        # as they can influence the first part of the curve too much
        valmax=curva1[-1]
        soglia=valmax*0.99
        mask=curva1<soglia
        curva= curva1[mask]
        # find the number of points left
        numpunti=len(curva)

        # creating an array even for the depths limited to the first few points
        xx=H_Array[:numpunti]
        x=np.log(xx)

        tratto=int(ListaAscisse[i])
        distanza=TrattoDistanza[tratto]
        distanza_fiume=TrattoDistanza_fiume[tratto]
        progressiva=TrattoProg_fiume[tratto]

        # reads the slope of the reach
        pend=TrattoPend[tratto]

        # dividing the volume curve by the length of the reach
        # we get the average transverse area of the reach
        curvaAreaMedia=curva/distanza

        # find the function monomy Area = const * h ^ m
        # we adopt the least squares technique
        # let us make linear the formula through the passage to logarithms
        # ln(A)=ln(c) + m*ln(h)
        # the angular coefficient m of the straight line between the logarithms is
        # equal to the exponent sought for the monomy formula
        y=np.log(curvaAreaMedia)

        # We can rewrite the line equation as y = Ap, where A = [[x 1]] and p = [[m], [c]]
        A = np.vstack([x, np.ones(len(x))]).T
        try:
            # Now use lstsq to solve for p:
            m, c = np.linalg.lstsq(A, y,rcond=None)[0]

            # check the case sections that do not widen with increasing depth
            if m<1.0:
##                m=1.5
                pass

            # doing the exponential of the logarithm of c we get the constant of the formula monomy
            const=math.exp(c)
        except:
            # in case of error use a fixed scale Area = const * h ^ m
            m=1.5
            const=10.0

        # enter a value> 0 if you want to view the graph
        graficoA=0
        if graficoA>0:
            yy=const*np.power(xx,m)
            plt.plot(xx, curvaAreaMedia, 'o', label='Original data', markersize=10)
            plt.plot(xx, yy, 'r', label='Fitted line')
            txt='tratto : %d distance: %d' % (i,tratto)
            plt.title(txt)
            txt='Area= %.2f *H^%.3f' % (const,m)
            asc=xx[0]
            ordin=plt.ylim()[1]*0.9
            plt.text(asc, ordin, txt, horizontalalignment='left',verticalalignment='top', fontsize=12, rotation=0 )

            plt.legend()
            plt.grid()
            plt.show()

        # calculate the rating curves by the Chezy formula v=X*radq(R*if) assuming that, according to Manning, X=1/n*R^(1/6)
        # -----------------------------------------------------------------------------------------------------------------
        # Q=1/nMann*A*R^2/3*radq(pend) =approximately 1/nMann*A*h^2/3*radq(pend)= 1/nMann*const**radq(pend)*h^(m+2/3)
        # Q=kq*h^mq
        kq=1.0/nMann*const*math.sqrt(pend)
        mq=(m+2.0/3.0)

        # calculate the celerity c=dQ/dA  - c=dQ/dh*dh/dA= kq*mq*h^(mq-1)*dh/dA
        # ----------------------------------------------------------------------
        # B: water width; h: water depth =approximately Hydraulic Radius
        # Q=B*h*1/nMann*h^(1/6)*radq(h*if) = B*1/nMann*h^(5/3)*radq(if)
        # B =dA/dh
        # celerity' c=dQ/dA = dQ/dh * dh/dA = 5/3*B*1/nMann*h^(2/3)*radq(if) * 1/B = 5/3*1/nMann*h^(2/3)*radq(if)
        # celerity=kcel*h^mcel
        kcel=5.0/3.0/nMann*math.sqrt(pend)
        mcel=2.0/3.0

        # calculate water width  B=dA/dh
        # ------------------------------------------------
        # B=const*m*h^(m-1) = kB*h^mB
        kB=const*m
        mB=m-1.0

        sql='INSERT INTO %s (DamID' % (NomeTabella)
        sql_value=') VALUES (%d' % DamID
        sql_value+=',%d,%.1f,%.1f,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f' %(tratto,progressiva,distanza_fiume,distanza,pend,const,m,kq,mq,kcel,mcel)
        sql+=', PixDist'
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
        sql+='%s' % sql_value
        sql+=');'
        cur.execute(sql)


        # calculation and graph of celerity
        # ...........................
        graficoc=0
        if graficoc>0:
            yy=kcel*np.power(xx,mcel)
            plt.plot(xx, yy, 'r', label='celerity')
            txt='tratto : %d distance: %d' % (i,tratto)
            plt.title(txt)
            txt='c= %.2f *H^%.3f' % (kcel,mcel)
            asc=xx[0]
            ordin=plt.ylim()[1]*0.9
            plt.text(asc, ordin, txt, horizontalalignment='left',verticalalignment='top', fontsize=12, rotation=0 )
            plt.xlabel('h (m)')
            plt.ylabel('c (m/s)')

            plt.legend()
            plt.grid()
            plt.show()

        # calculation and graph of water width
        # ............................................
        graficoB=0
        if graficoB>0:
            yy=kB*np.power(xx,mB)
            plt.plot(xx, yy, 'r', label='water width')
            txt='tratto : %d distance: %d' % (i,tratto)
            plt.title(txt)
            txt='B= %.2f *H^%.3f' % (kB,mB)
            asc=xx[0]
            ordin=plt.ylim()[1]*0.9
            plt.text(asc, ordin, txt, horizontalalignment='left',verticalalignment='top', fontsize=12, rotation=0 )
            plt.xlabel('h (m)')
            plt.ylabel('B (m)')

            plt.legend()
            plt.grid()
            plt.show()


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
    fileDEM=  PathFiles+ os.sep+'DTM_clip.tif'


    NotErr, errMsg= CurveAreaAltezza(mydb_path_user,DamID,PathFiles)

    print(NotErr,errMsg)

