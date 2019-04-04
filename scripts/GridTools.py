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
    # importa i sistemi di riferimento
    from osgeo.osr import osr
except:
    import ogr
    # importa i sistemi di riferimento
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


def GridDistanzaFiume(mydb_path_user,ID_Diga,PathFiles,ClipDEM,DistanzaMax=4000.0):

    """
    Lo script legge:
    - il grid del ClipDEM
    - lo shape file delle sezioni totali
    - lo shape la linea dell'asse del fiume verso valle
    crea:
    - lo shape file (in memoria) dei poligoni a diversa distanza dalla linea dell'asse
    - un grid Distances.tif con in valori delle distanze rispetto al fiume

    Dato input:
        - DistanzaMax: larghezza massima del flusso di piena per cui si ammette
                       il moto monodimensionale. Oltre tale distanza eventuali
                       allagamenti si ipotizza che siano zone in cui l'onda si
                       si espande ma il flusso è statico (es. vasca di espansione)
    Crea:

    - grid Distances.tif il grid delle distanze dal fiume

    - un grid Tratti.tif con le classi dei tratti di fiume a partire dalla diga
      il numero del tratto rappresenta in numero di celle di
      distanza contate lungo l'asse del fiume

    Utilizza il parametro DistanzaMax per eliminare da Tratti.tif il pixel oltre
    la  DistanzaMax dal fiume

    """

    NotErr=bool('True')
    errMsg='OK'

    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    CrossSec=PathFiles+os.sep+'CrossSec.shp'
    CrossSecTot=PathFiles+os.sep+'CrossSecTot.shp'

    if not os.path.exists(CrossSecTot):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    if not os.path.exists(ClipDEM):
        errMsg = "Manca per la diga num =%s il clip DEM\nEffettuare prima il ritaglio del modello digitale del terreno !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    CrossSecPoly='%sPoly.shp' % (CrossSec[:-4])
    if not os.path.exists(CrossSecPoly):
        errMsg = "Manca per la diga num =%s il Grid Tratti\nEffettuare prima il CreaSezInterpolate !" % (ID_Diga)
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
        TotalLength=ChkDiga[0]


    # Close communication with the database
    cur.close()
    conn.close()

    # nome del file output
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
        pass
        # definisco un default per il riferimento
        # imposto WGS84 UTM 32 N
        spatialRef.ImportFromEPSG(32632)

    terreno = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)
    mask_Nodata= terreno==inNoData

    gt=indataset.GetGeoTransform()

    inband = None

    indataset = None


    # carico il driver dello shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')


    # leggo lo shape file delle sezioni
    InDS = driver.Open(CrossSecTot, 0)
    if InDS is None:
    	print ('Could not open ' + CrossSec)
    	sys.exit(1)     #exit with an error code


    # leggo il layer dalla sorgente dei dati
    Inlayer = InDS.GetLayer()

    inFeature=Inlayer.GetNextFeature()

    Lengthtmax=0.0

    while inFeature:
        NumSez=inFeature.GetField('id')
        geom = inFeature.GetGeometryRef()
        length=geom.Length()
        if length>Lengthtmax:
            Lengthtmax=length
        # ogr.CreateGeometryFromWkt(wkt)
        inFeature=Inlayer.GetNextFeature()

    InDS.Destroy()

    # calcola il numero di steps
    NumSteps=int(Lengthtmax/pixelWidth/1.5)

    ListaDist=[]
    for i in range(NumSteps):
        dist=(i+1)*pixelWidth
        ListaDist.append(dist)

    shpnew2=PathFiles+os.sep +"distances.shp"
    nomecampoDist='dist'


    PathGeom=ogr.CreateGeometryFromWkt(wkt_line)


    # creo in memoria
    # ---------------

    #create an output datasource in memory
    outdriver=ogr.GetDriverByName('MEMORY')
    outDS2=outdriver.CreateDataSource('memData')

    #open the memory datasource with write access
    tmp=outdriver.Open('memData',1)

    outLayer2 = outDS2.CreateLayer('distances', dest_srs,geom_type=ogr.wkbPolygon)

    # crea i nuovi campi id e Nome nello shapefile di output
    fieldDefn1 = ogr.FieldDefn('id', ogr.OFTInteger)
    outLayer2.CreateField(fieldDefn1)
    fieldDefn2 = ogr.FieldDefn(nomecampoDist, ogr.OFTReal)
    outLayer2.CreateField(fieldDefn2)

    featureDefn2 = outLayer2.GetLayerDefn()

    # prima distanza
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
        # crea la geometria precedente
        poly_old=ogr.CreateGeometryFromWkt(poly_old_wkt)
        # salva la geometria attuan
        poly_old_wkt=NewGeom1.ExportToWkt()
        # effettua la differenza della geomatria attuane con la geometria precedente
        # symmetric difference
        simdiff = NewGeom1.SymmetricDifference(poly_old)
        area=simdiff.Area()

        feature = ogr.Feature(featureDefn2)
        feature.SetField('id', i)
        feature.SetField(nomecampoDist,ListaDist[i])
        feature.SetGeometry(simdiff)
        outLayer2.CreateFeature(feature)
        simdiff.Destroy()

    # leggo dalla memoria
    Inlayer =outLayer2


    # rasterizza
    # -------------------
    format = 'GTiff'
    type = GDT_Float32

    driver2 = gdal.GetDriverByName(format)
    driver2.Register()

    # matrice massime distanza
    DistMax=np.zeros((rows,cols),np.float32)
    DistMax=DistMax+ListaDist[-1]

    ds = driver2.Create(FileDEM_out, cols,  rows, 1, type)
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
    testo="ATTRIBUTE=%s"  % (nomecampoDist)
    # Rasterize
    outband = ds.GetRasterBand(iBand)

    # Rasterize
    # inizializzo con distanza massima
    outband.WriteArray(DistMax, 0, 0)
    CampoValore=[testo]

    # creo la mappa dei valori
    # -------------------------------------------------
    err = gdal.RasterizeLayer(ds, [iBand], Inlayer,
            burn_values=[0],
            options=CampoValore)
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    # Reading
    dist = outband.ReadAsArray().astype(np.float32)
    # aggiungo i nodata
    Nodata=-9999
    dist= np.choose(mask_Nodata,(dist,Nodata))
    # salvo
    outband.WriteArray(dist, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(Nodata)
    outband.GetStatistics(0,1)
    outband = None

    ds= None

    # creazione del grid dei tratti di fiume
    # ---------------------------------------
    CreaCrid=1

    if CreaCrid>0:

        # creo la maschera della distanza massima di analisi geometria flusso monodimensionale
        mask_Dist= np.greater_equal(dist,DistanzaMax)


        # ==================
        # Legge i dati GRID
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
            # definisco un default per il riferimento
            # imposto WGS84 UTM 32 N
            spatialRef.ImportFromEPSG(32632)

        terreno = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)
        mask_Nodata= terreno==inNoData


        # rasterizza i poligoni
        # ====================
        orig_data_source = ogr.Open(CrossSecPoly)
        source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_source, "")
        source_layer = source_ds.GetLayer()

        # legge la lista degli id
        # -----------------------
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

##        format = 'MEM'
        format = 'Gtiff'
##        type = GDT_Float32
        type = GDT_Int16

        driver3 = gdal.GetDriverByName(format)
        driver3.Register()

        PathFiles=os.path.dirname(ClipDEM)
        TrattiRaster=PathFiles+os.sep+'Tratti.tif'

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
        ClassTratti=np.zeros((rows,cols)).astype(np.int)
        ClassTratti=ClassTratti-1

        outband.WriteArray(ClassTratti, 0, 0)

        # Rasterize
        err = gdal.RasterizeLayer(dsRaster, [1], source_layer,
                burn_values=[0],
                options=["ATTRIBUTE=id"])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)

        # scrive i Nodata
        MatriceDati=outband.ReadAsArray(0, 0, cols, rows)
        MatriceDati=np.choose(mask_Nodata,(MatriceDati,outNodata))

        # scarta i punti oltre la distanza massima assegnando Nodata
        MatriceDati=np.choose(mask_Dist,(MatriceDati,outNodata))

        outband.WriteArray(MatriceDati, 0, 0)

        outband.FlushCache()
        outband.SetNoDataValue(outNodata)
        outband.GetStatistics(0,1)

        outband=None

        dsRaster=None
        # chiude le sorgenti dei dati
        orig_data_source.Destroy()



    return NotErr, errMsg

def ModDTM(mydb_path_user,ID_Diga,PathFiles,ClipDEM):

    """
    Lo script legge:
        - il grid del StreamDH
        - il grid del Distances
        - lo shape la linea dell'asse del fiume verso valle
    crea:
        - il grid del StreamDHFilled: altezze che
          abbiamo almeno un pendenza rispetto al fiume >=PendMin
          e che lungo l'asse del fiume abbiano valore =0

    #-------------------------
    #------- medodologia -----
    # ------------------------

    1) lettura dei due grid di input
    2) creazione di una matrice delle altezze minime che rispettano
       la condizione pendenza_rispetto_al_fiume=PendMin

       Hmin=DistancesArray*PendMin  (operazione fra array di numpy)

    3) ricerca di una maschera dei punti di StreamDH in input che non soddisfano la
       condizione

        mask= StreamDH < PendMin (operazione fra un array di numpy ed il valore PendMin)

    4) sostituzione nei punti della maschera con i valori Hmin

        StreamDHFilled[mask] = Hmin[mask]  (operazione fra array di numpy)

    5) creazione del grid "stream" dell'asse del fiume avente le stesse dimensioni e gridsize
       di StreamDH ed avente valore 1 sull'asse

    6) sostituzione in StreamDHFilled del valore 0 nei punti in cui stream=1

        maskstream=stream==1

        StreamDHFilled=np.choose(maskstream,(StreamDHFilled,0.0))

    7) salvo la matrice StreamDHFilled finale nel file tif = StreamDHFilled.tif

    """

    NotErr=bool('True')
    errMsg='OK'

    # assegno una pendenza minima es. 1/mille
    PendMin=float(0.002)

    PathFiles=os.path.realpath(PathFiles)

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDH.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Manca per la diga num =%s il StreamDH\nEffettuare prima CreaSezInterpolate !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    Distances=PathFiles+os.sep+'Distances.tif'
    if not os.path.exists(Distances):
        errMsg = "Manca per la diga num =%s il Distances\nEffettuare prima il CreaGridDistanza !" % (ID_Diga)
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
        TotalLength=ChkDiga[0]


    # Close communication with the database
    cur.close()
    conn.close()

    # nome del file output
    FileDEM_out=PathFiles+os.sep+'StreamDHFilled.tif'

    # ==================================
    # Legge i dati GRID dello STREAM DH
    # ==================================

    gdal.AllRegister()

    indataset = gdal.Open(StreamDH, GA_ReadOnly )
    if indataset is None:
##        print ('Could not open ' + StreamDH)
##        sys.exit(1)
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
        # definisco un default per il riferimento
        # imposto WGS84 UTM 32 N
        spatialRef.ImportFromEPSG(32632)

    StreamDHA = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)

    # crea la maschera della zona con Nodata
    mask_Nodata= StreamDHA==inNoData

    gt=indataset.GetGeoTransform()

    inband = None

    indataset = None

    # ==================================
    # Legge i dati GRID del DISTANCES
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
        # definisco un default per il riferimento
        # imposto WGS84 UTM 32 N
        spatialRef.ImportFromEPSG(32632)

    DistancesA = inband.ReadAsArray(0, 0, cols, rows).astype(np.float)

    gt=indataset.GetGeoTransform()

    inband = None

    indataset = None

    # PT 2 --------------------

    MatriceHmin = DistancesA*PendMin

    # PT 3 e 4 ----------------

    # ricerca di una maschera dei punti di StreamDH in input che non soddisfano
    # la condizione mask= StreamDH < PendMin (operazione fra un array di numpy ed il valore PendMin)

    Inf=StreamDHA<MatriceHmin

    StreamDHFilled=np.choose(Inf,(StreamDHA,MatriceHmin))

    StreamDHFilled=np.choose(mask_Nodata,(StreamDHFilled,inNoData))

    # punto 5 ----------------- creazione del grid dell'asse del fiume leggo lo shape file dell'asse del fiume


    # creo in memoria
    # ---------------

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


    # rasterizza
    # -------------------
##    format = 'GTiff'
    # sostituendo al format='GTiff' il format = 'MEM' la matrice
    # viene salvata solo in memoria e non sul disco
    format = 'MEM'
##    type = GDT_Float32
    type = GDT_Int16

    driver2 = gdal.GetDriverByName(format)
    driver2.Register()

    FileDEM_tmp=PathFiles+os.sep+'Stream.tif'

    ds = driver2.Create(FileDEM_tmp, cols,  rows, 1, type)
    if gt is not None and gt != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
        ds.SetGeoTransform(gt)

    # imposta il sistema di riferimento uguale al modello del tiranti: se manca imposta il default
    if prj is not None and len(prj) > 0:
        ds.SetProjection(prj)
    else:
        prj= spatialRef.ExportToWkt()
        ds.SetProjection(prj)

    iBand=1
    # Rasterize
    outband = ds.GetRasterBand(iBand)

    # creo la mappa dei valori
    # ---------------------------------------------
    # nota burn_values=[1] indica di rasterizzare con valore 1 il vettore
    err = gdal.RasterizeLayer(ds, [iBand], outLayer,
            burn_values=[1])
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)

    # Reading
    stream = outband.ReadAsArray().astype(np.int)
    # aggiungo i nodata
    Nodata=-9999
    stream= np.choose(mask_Nodata,(stream,Nodata))
    # salvo
    outband.WriteArray(stream, 0, 0)

    outband.FlushCache()
    outband.SetNoDataValue(Nodata)
    outband.GetStatistics(0,1)
    outband = None

    ds= None

    # implementare i punti 6 e 7 della metodologia

    # PT 6 ------------- sostituzione in StreamDHFilled del valore 0 nei punti in cui stream=1

    maskstream=stream==1

    StreamDHFilled=np.choose(maskstream,(StreamDHFilled,0.0))

    # PT 7 ------------- salvo la matrice StreamDHFilled finale nel file tif = StreamDHFilled.tif

    # nome del file output
    PathFiles=os.path.dirname(StreamDH)
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

    # imposta il sistema di riferimento uguale al modello precedente
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

def CurveAreaAltezza(mydb_path_user,ID_Diga,PathFiles):

    """
    Lo script legge:
        - StreamDHFilled.tif : grid con in valori delle altezze del terreno rispetto al fiume
        - Tratti.tif         : grid con le classi dei tratti di fiume a partire dalla diga
                              il numero del tratto rappresenta in numero di celle di
                               distanza contate lungo l'asse del fiume

    effettua il conteggio, per ogni tratto, e per ogni dh con passo 1 metro, nel numero
    di celle soggiacenti la differenza di altezza dh dal fiume
    Costruisce anche le curve del numero di celle cumulate lungo il percorso del
    fiume per ogni dh con passo 1 metro

    Salva i due conteggi in due file csv:
        - MatricePixel.csv
        - MatricePixCum.csv
    """


    NotErr=bool('True')
    errMsg='OK'

    PathFiles=os.path.realpath(PathFiles)

    # files output
    filecsv1=PathFiles+os.sep+'MatricePixel.csv'
    filecsv2=PathFiles+os.sep+'MatricePixCum.csv'

    # contiene la percentuale di area in destra orografica
    filecsvPixDestra=PathFiles+os.sep+'MatricePixDestra.csv'

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDHFilled.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Manca per la diga num =%s il grid StreamDHFilled\nEffettuare prima ModificaDH !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    Tratti=PathFiles+os.sep+'Tratti.tif'
    if not os.path.exists(Tratti):
        errMsg = "Manca per la diga num =%s il Grid Tratti\nEffettuare prima il CreaSezInterpolate !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    DestraSinistra=PathFiles+os.sep+'DestraSinistra.tif'
    if not os.path.exists(DestraSinistra):
        errMsg = "Manca per la diga num =%s il Grid DestraSinistra\nEffettuare prima il CreaSezInterpolate !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg



    # ==================
    # Legge i dati GRID
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
        # definisco un default per il riferimento
        # imposto WGS84 UTM 32 N
        spatialRef.ImportFromEPSG(32632)

    TrattiArray = inband.ReadAsArray(0, 0, cols, rows).astype(np.int)
##    mask_Nodata= TrattiArray==inNoData

    TrattiVector_ini=np.unique(TrattiArray)

    inband = None
    indataset = None

    # lettura StreamDH
    # -----------------

    if not os.path.exists(StreamDH):
        errMsg = "File StreamDH %s non esiste" % os.path.realpath(StreamDH)
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
        errMsg = 'Grid non conguente: %s' % infile
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

    # assegno, per congruenza, la mappa Nodata dei DH anche ai TrattiArray
    TrattiArray=np.choose(mask_Nodata,(TrattiArray,inNoData))

    TrattiVector_ini2=np.unique(TrattiArray)

    # elimina i  valori negativi
    DH=np.choose(np.less(DH,0.0),(DH,0.0))

    inbandElev=None

    indatasetElev= None

    # legge grid DestraSinistra
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
        # definisco un default per il riferimento
        # imposto WGS84 UTM 32 N
        spatialRef.ImportFromEPSG(32632)

    # lettura matrice zona destra fluviale =1
    TrattiArrayDx = inbandDx.ReadAsArray(0, 0, cols, rows).astype(np.int)

    # creo la maschera della zona destra
    mask_Dx=np.equal(TrattiArrayDx,1)
    numDx=mask_Dx.sum()

    inbandDx = None
    indatasetDx = None


    # crea la lista dei tratti
    # ------------------------
    TrattiVector1=np.unique(TrattiArray)
    # scarta i dati <0 (Nodata e dati esterni ai tratti)
    mask=np.where(TrattiVector1>0)[0]

    # vettore finale
    TrattiVector=TrattiVector1[mask]

    # Massima altezza per cui vengono create le curve
    Hmax=51

    MatricePix=[]
    VettorePixHmax=[]
    VettoreVolHmax=[]

    # Matrice percentuale dx
    MatricePixDx=[]

    # scorre i tratti
    for tratto in TrattiVector:
        # maschera del tratto
        mask_tratto=np.equal(TrattiArray,tratto)
        # numero di pichel del tratto
        numeropixel=mask_tratto.sum()
        VettorePix=[]
        VettorePixDx=[]
        for h in range(1,Hmax):
            mask=np.less_equal(DH,h) & mask_tratto
            nn=mask.sum()
            if nn>0:
                # seleziona i pixel del tratto < di h
                DH_cur=DH[np.where(mask)]
                # somma le altezze che equivale al volume in termini di pixel*h
                Vol_h=np.sum(np.absolute(DH_cur),dtype=np.float32)
                # il volume vuoto si ottiene per la differenza
                Vol_d_valle=float(nn)*h-Vol_h

                VettorePix.append(Vol_d_valle)

                # trovo quelli a destra idrografica
                # ---------------------------------
                mask2=mask & mask_Dx
                ndx=mask2.sum()
                # seleziona i pixel del tratto < di h
                DH_curDx=DH[np.where(mask2)]
                # somma le altezze che equivale al volume in termini di pixel*h
                Vol_h_Dx=np.sum(np.absolute(DH_curDx),dtype=np.float32)
                # il volume vuoto si ottiene per la differenza
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

    # crea la matrice delle aree cumulate a parita' di altezza dall'alveo
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


    # salva le matrici MatricePixel & MatricePixDx
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

    # salva la matrice MatricePixCum
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
        # GRAFICO 1
        # ----------
        #prepara i fonts per eventuali grafici di controllo
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
        ax.set_xlabel('Distanza (km)',color='blue')
        txt='Area della Valle per diverse altezze rispetto al fiume a partire dalla diga\n'
        plt.title(txt)

        plt.grid()
        plt.show()

        # ---------
        # GRAFICO 2
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
        ax.set_xlabel('Distanza (km)',color='blue')
        txt='Area Cumulata della Valle per diverse altezze rispetto al fiume a partire dalla diga\n'
        plt.title(txt)

        plt.grid()
        plt.show()

    return NotErr, errMsg

def ParamGeomHydro(mydb_path_user,ID_Diga,PathFiles):
    """
    Lo script calcola la curve delle conveyance per una sezione rappresentativa del tratto
    Legge le distanze ed i dislivelli dal file dei punti lungo l'asse del fiume
     - PathPoints.shp

    Legge dalla lista di parametri passati in input:
     - cellsize : dimensione della cella del raster di partenza
     - nMann    : numero di Manning assunto constante
     - pendMin  : pensenza minima accettabile, in caso di pendenze inferiori assume
                  pendMin

    Legge le curve delle aree dei tratti dal file
     - MatricePixel.csv

    Calcola le distanze e le pendenze di ciascun tratto dallo shape file PathPoints.shp
    Partendo dai punti delle aree dei tratti, trova con i minimi quadrati una formula
    monomia per approssimare la variazione di area con l'altezza d'acqua

    Area=ka*h^ma

    da questa, applicando la formula del moto uniforme calcola un formula monomia
    della scala di delfusso, assimilando il raggio idraulico all'altezza h

    Q= 1/nMann*ka**radq(pend)*h^(ma+2/3) = kq*h^mq

    Calcola la celerita' c=dQ/dA sempre come formula monomia

    c=kcel*h^mcel

    Calcola la larghezza del pelo libero B derivando dA/dh

     B=ka*ma*h^(ma-1) = kB*h^mB

    Salva i coefficienti e gli esponenti trovati delle formule monomia nel file:

        - MatriceAexp.csv

    """

##    dir_script= os.path.dirname(os.path.realpath(__file__))
##    os.chdir(dir_script)
    pythonver_log=sys.version
    ppp=pythonver_log.split(' ')
    py_ver=ppp[0][0:1]


    NotErr=bool('True')
    errMsg='OK'

    PathFiles=os.path.realpath(PathFiles)

    # assegno una pendenza minima accettabile pari all'uno per mille
    pendMin=0.001
    # assumo Manning costante
    nMann=0.06

##    filecsv_out=PathFiles+os.sep+'MatriceAexp.csv'

    if not os.path.exists(PathFiles):
        errMsg = "Non ci sono dati per la diga num =%s \nEffettuare prima il calcolo delle sezioni a valle !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    StreamDH=PathFiles+os.sep+'StreamDHFilled.tif'
    if not os.path.exists(StreamDH):
        errMsg = "Manca per la diga num =%s il StreamDH\nEffettuare prima ModificaDH !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg


    filecsv1=PathFiles+os.sep+'MatricePixel.csv'
    if not os.path.exists(filecsv1):
        errMsg = "Manca per la diga num =%s il MatricePixel.csv\nEffettuare prima CreaCurveAreaAltezza !" % (ID_Diga)
        NotErr= bool()
        return NotErr, errMsg

    # =====================================
    # apertura del database sqlite utente
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

    # codice del sistema di riferimento della tabella
    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (TabellaPoints.lower())
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
    # cerco cellsize
    # ==================================

    # Legge i dati GRID dello STREAM DH
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


    # legge shape file
    # carico il layer dello shapefile PathPoints.shp
##    driver = ogr.GetDriverByName('ESRI Shapefile')


    # conto il numero di Feature
    n = len(ListaTratti)

    # leggo la prima feature
    id0=ListaTratti[0][0]
    geom=ogr.CreateGeometryFromWkt(NumSezGeom[id0])

    # leggo le coordinate
    x0=geom.GetX()
    y0=geom.GetY()
    # leggo i dati
    elev0=NumSezElev[id0]
    ProgrMonte=NumSezProgr[id0]

    # dizionari
    TrattoPend={}
    # dizionario lunghezza del tratto in linea retta
    TrattoDistanza={}
    # dizionario lunghezza del tratto lungo il fiume
    TrattoDistanza_fiume={}
    # dizionario progressiva media del tratto lungo il fiume
    TrattoProg_fiume={}

    numtratti=len(ListaTratti)

    for itratto in range(1,numtratti):

        # leggo i dati
        id1=ListaTratti[itratto][0]
        elev1=NumSezElev[id1]
        ProgrValle=NumSezProgr[id1]


        # leggo la sua geometria
        geom = ogr.CreateGeometryFromWkt(NumSezGeom[id1])
        # leggo le coordinate
        x1=geom.GetX()
        y1=geom.GetY()
        # trovo la distanza in linea retta fra i due punti
        distanza=math.sqrt(math.pow((x1-x0),2)+math.pow((y1-y0),2))
        # trovo la distanza parziale come differenza delle progressive
        distanza_fiume=ProgrValle-ProgrMonte
        # trovo la distanza progressiva intermedia
        ProgrMedia=(ProgrMonte+ProgrValle)/2.0

        # la pendenza e' la tangente dell'angolo
        pend=(elev0-elev1)/distanza

        # controllo pendenza 0
        if pend<=pendMin:
            pend=pendMin*1.0

        # salvo nei dizionari
        TrattoPend[id1]=pend
        TrattoDistanza[id1]=distanza
        TrattoDistanza_fiume[id1]=distanza_fiume
        TrattoProg_fiume[id1]=ProgrMedia

##        # passo al punto successivo
        # il punto 1 di un tratto diventa il punto 0 del tratto seguente
        x0=x1*1.0
        y0=y1*1.0
        elev0=elev1*1.0
        ProgrMonte=ProgrValle*1.0


    # lettura Curve dei volumi del tratto
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

    # vettore delle altezze d'acqua
    H_Array=np.array(Vettoreh,dtype =np.float)

    ListaAscisse=[]
    MatricePix=[]
    for row in reader:
        ListaAscisse.append(row[0])
        MatricePix.append(row[1:])

    fin.close()

    # matrice delle aree intese come numeropixel*h per ogni dh rispetto al fiume
    MatriceArray=np.array(MatricePix,dtype =np.float)

    # area della cella
    cellArea=cellsize*cellsize

    # numero di celle che interessano la lunghezza del primo tratto
    # questo e' anche il passo costante con cui sono costruiti anche
    # i successivi tratti
    deltacelle=int(ListaAscisse[0])

    # trova le zone dove l'area e' inferiore alla minima:
    # si assume come minimo un numero di celle pari a quelle del percorso
    # dell'asse del fiume in un tratto
    zonamin=MatriceArray<deltacelle
##    n1=zonamin.sum()
    # elimina zone areamin
    MatriceArray=np.choose(zonamin,(MatriceArray,deltacelle))
##    # verifica
##    zonamin=MatriceArray<deltacelle
##    nn=zonamin.sum()
    # Il volume nella MatriceArray e' numero di celle x altezza
    # il volume effettivo in mq si ottiene moltiplicando per l'area della cella
    MatriceAree=MatriceArray*cellArea


    # -----------------------------------------------------------------------------
    # calcola per ogni tratto la curva monomia dell'Area in funzione dell'altezza
    # -----------------------------------------------------------------------------

    # salva la matrice MatriceAexp
##    fout=open(filecsv_out,'w')
##    txt='PixDist;progr_lungo_fiume;distanza_fiume;distanza_linea_retta;pend;ka;ma;kq;mq;kcel;mcel\n'
##    fout.write(txt)

##    # salva nel database
##    conn = sqlite3.connect(mydb_path_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
##
##    # conn.cursor will return a cursor object, you can use this cursor to perform queries
##    cur= conn.cursor()

    # cancella eventuali dati precedenti
    # -----------------------------------
    NomeTabella='MatriceAexp'
    sql='DELETE FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
    cur.execute(sql)
    conn.commit()

    nn=len(MatriceAree)

    # grafico Area
    graficoA=0
    # grafico larghezza pelo libero
    graficoB=0

    # debug
##    ListaDebug=[]
##    ListaDebug.append(6)
##    ListaDebug.append(7)
##    ListaDebug.append(8)
##    ListaDebug.append(9)
##    # grafico Area
##    graficoA=1
##    graficoB=1
##    for i in ListaDebug:
    # fine debug

    for i in range(nn):
        curva1= MatriceAree[i]

        # elemina i punti estremi
        # in quanto possono influenzare troppo la prima parte della curva
        valmax=curva1[-1]
        soglia=valmax*0.99
        mask=curva1<soglia
        curva= curva1[mask]
        # trovo il numero di punti rimasti
        numpunti=len(curva)

        # creo un array anche per le altezze limitato ai primi numpunti
        xx=H_Array[:numpunti]
        x=np.log(xx)

        tratto=int(ListaAscisse[i])
        distanza=TrattoDistanza[tratto]
        distanza_fiume=TrattoDistanza_fiume[tratto]
        progressiva=TrattoProg_fiume[tratto]

        # legge la pendenza del tratto
        pend=TrattoPend[tratto]

        # dividendo la curva dei volumi per la lunghezza del tratto
        # ottengo l'area trasversale media del tratto
        curvaAreaMedia=curva/distanza

        # trova la funzione monomia Area=const*h^m
        # adotto la tecnica dei minimi quadrati
        # rendo lineare la formula monomia mediante il passaggio ai logaritmi
        # ln(A)=ln(c) + m*ln(h)
        # il coefficiente angolare m della retta fra i logarirmi e'
        # pari all'esponente cercato della formula monomia
        y=np.log(curvaAreaMedia)

        # We can rewrite the line equation as y = Ap, where A = [[x 1]] and p = [[m], [c]]
        A = np.vstack([x, np.ones(len(x))]).T
        try:
            # Now use lstsq to solve for p:
            m, c = np.linalg.lstsq(A, y,rcond=None)[0]

            # controllo caso sezioni che non si allargano all'aumentare dell'altezza
            if m<1.0:
##                m=1.5
                pass

            # facendo l'esponenziale del logaritmo di c ottendo la costante della formula monomia
            const=math.exp(c)
        except:
            # in caso di errore usa una scala fissa Area=const*h^m
            m=1.5
            const=10.0

        # inserire un valore >0 se si vuole visionare il grafico
        graficoA=0
        if graficoA>0:
            yy=const*np.power(xx,m)
            plt.plot(xx, curvaAreaMedia, 'o', label='Original data', markersize=10)
            plt.plot(xx, yy, 'r', label='Fitted line')
            txt='tratto : %d distanza: %d' % (i,tratto)
            plt.title(txt)
            txt='Area= %.2f *H^%.3f' % (const,m)
            asc=xx[0]
            ordin=plt.ylim()[1]*0.9
            plt.text(asc, ordin, txt, horizontalalignment='left',verticalalignment='top', fontsize=12, rotation=0 )

            plt.legend()
            plt.grid()
            plt.show()

        # calcola la scala di deflusso con la  formula di Chezy v=X*radq(R*if) assumendo secondo Manning X=1/n*R^(1/6)
        # -----------------------------------------------------------------------------------------------------------
        # Q=1/nMann*A*R^2/3*radq(pend) =circa 1/nMann*A*h^2/3*radq(pend)= 1/nMann*const**radq(pend)*h^(m+2/3)
        # Q=kq*h^mq
        kq=1.0/nMann*const*math.sqrt(pend)
        mq=(m+2.0/3.0)

        # calcolo della celerita' c=dQ/dA  - c=dQ/dh*dh/dA= kq*mq*h^(mq-1)*dh/dA
        # ----------------------------------------------------------------------
        # B: larghezza pelo libero; h: altrezza =circa Raggio Idraulico
        # Q=B*h*1/nMann*h^(1/6)*radq(h*if) = B*1/nMann*h^(5/3)*radq(if)
        # B =dA/dh
        # celerita' c=dQ/dA = dQ/dh * dh/dA = 5/3*B*1/nMann*h^(2/3)*radq(if) * 1/B = 5/3*1/nMann*h^(2/3)*radq(if)
        # celerita=kcel*h^mcel
        kcel=5.0/3.0/nMann*math.sqrt(pend)
        mcel=2.0/3.0

        # calcolo della larghezza del pelo libero B=dA/dh
        # ------------------------------------------------
        # B=const*m*h^(m-1) = kB*h^mB
        kB=const*m
        mB=m-1.0

        sql='INSERT INTO %s (ID_Diga' % (NomeTabella)
        sql_value=') VALUES (%d' % ID_Diga
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


        # calcolo e grafico celerita'
        # ...........................
        graficoc=0
        if graficoc>0:
            yy=kcel*np.power(xx,mcel)
            plt.plot(xx, yy, 'r', label='celerita')
            txt='tratto : %d distanza: %d' % (i,tratto)
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

        # calcolo e grafico larghezza del pelo libero
        # ............................................
        graficoB=0
        if graficoB>0:
            yy=kB*np.power(xx,mB)
            plt.plot(xx, yy, 'r', label='larghezza pelo libero')
            txt='tratto : %d distanza: %d' % (i,tratto)
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


##    fout.close()

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
    fileDEM=  PathFiles+ os.sep+'DTM_clip.tif'


##    NotErr, errMsg= GridDistanzaFiume(mydb_path_user,ID_Diga,PathFiles,fileDEM)
##
##    print(NotErr,errMsg)

##    NotErr, errMsg= ModDTM(mydb_path_user,ID_Diga,PathFiles,fileDEM)

    NotErr, errMsg= CurveAreaAltezza(mydb_path_user,ID_Diga,PathFiles)

    print(NotErr,errMsg)

