# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SimpleDamBrk - main

        begin                : 2019-02-19
        git sha              : $Format:%H$
        copyright            : (C) 2019 by L. Mancusi /RSE S.p.A
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
import os
import shutil
import sqlite3
import time

try:
    from osgeo import gdal
    gdal.TermProgress = gdal.TermProgress_nocb
    from osgeo.gdalconst import *
except ImportError:
    import gdal
    from gdalconst import *

try:
    from osgeo import ogr
    # import reference systems module
    from osgeo import osr
except:
    import ogr
    # import reference systems module
    import osr

from InterpolateCrossSec import SetIntermediatePoints, SetCrossSec_2
from GridTools import GridDistanzaFiume, ModDTM, CurveAreaAltezza, ParamGeomHydro
from RoutingCinemat import RunCinemat
from CalFloodArea import  GridFloodingAreas

class SimpleDambrk:
    """
    ============================================
    SimpleDambrk:
    ============================================
    """
    def __init__(self,myGDB_user):

        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)

        root_dir = self.plugin_dir + os.sep + '..'

        self.root_dir = os.path.realpath(root_dir)

        self.myGDB_user=os.path.realpath(myGDB_user)

        self.DamID=-1
        self.DTM=''
        self.PathFiles=''

        conn = sqlite3.connect(myGDB_user, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)

        # conn.cursor will return a cursor object, you can use this cursor to perform queries
        cur= conn.cursor()

        # import extention
        conn.enable_load_extension(True)
        conn.execute('SELECT load_extension("mod_spatialite")')

        self.conn=conn
        self.cur=cur


    def set_dam(self,DamID):

        """
        Set the current dam to be calculated
        """
        NotErr=bool('True')
        errMsg='Ok'

        NameTabella='DAMS'

        sql='SELECT Name FROM %s WHERE DamID=%d' % (NameTabella,DamID)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        try:
            self.Name=ChkDiga[0]
            self.DamID=DamID
            self.PathFiles=self.root_dir+os.sep+ str(self.DamID)
        except:
            errMsg = "Error %s dam data does not exists" % DamID
            NotErr= bool()

        return NotErr, errMsg


    def add_dam(self,DamID,Name,DamType,ResVolMcm,Height_m,BreachWidth,x,y):

        """
        Add the data of a new dam
        """

        NameTabella='DAMS'

        # Reference system check
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NameTabella.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        # check if it already exists
        sql='SELECT Name FROM %s WHERE DamID=%d' % (NameTabella,DamID)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        pt=ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(x, y)
        pt.FlattenTo2D()
        WKT_Point=pt.ExportToWkt()

        GeomWKT="GeomFromText('%s',%d)" % (WKT_Point,TargetEPSG)

        if ChkDiga==None:
            sql='INSERT INTO %s (DamID,Name,DamType,ResVolMcm,Height_m,BreachWidth,geom) VALUES (%d,"%s","%s",%s,%s,%s'  %  (NameTabella,DamID,Name,DamType,ResVolMcm,Height_m,BreachWidth)
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET Name="%s", DamType="%s", ResVolMcm=%s, Height_m=%s, BreachWidth=%s'  %  (NameTabella,Name,DamType,ResVolMcm,Height_m,BreachWidth)
            sql+=', geom=%s' % GeomWKT
            sql+=" WHERE DamID=%d;" % DamID
        self.cur.execute(sql)

        self.conn.commit()

        self.DamID=DamID
        self.Name=Name
        PathFiles=self.root_dir+os.sep+ str(self.DamID)
        self.PathFiles=os.path.abspath(PathFiles)
        if not os.path.exists(self.PathFiles):
            os.mkdir(self.PathFiles)

    def add_dam_from_shp(self,DamID,shpfile):

        """
        Add dam data from a shapefile
        """
        NotErr=bool('True')
        errMsg='OK'

        NameTabella='DAMS'

        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NameTabella.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        FieldsList=[]
        FieldsListType={}


        sql='PRAGMA table_info(%s);' % NameTabella
        self.cur.execute(sql)
        records=self.cur.fetchall()
        for rec in records:
            fieldname=rec[1]
            if  fieldname !='DamID' and fieldname !='PKUID' and fieldname !='geom':
                FieldsList.append(fieldname)
                FieldsListType[fieldname]=rec[2]



        # Reading shapefile path
        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        # load shapefile driver
        driver = ogr.GetDriverByName('ESRI Shapefile')

        # open the shapefile
        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg


        # reading source ayer
        layer = ds.GetLayer()
        layerDefinition = layer.GetLayerDefn()
        shp_field_list=[]
        for i in range(layerDefinition.GetFieldCount()):
            shpfieldName =  layerDefinition.GetFieldDefn(i).GetName()
            shp_field_list.append(shpfieldName)

        try:
            filtro="DamID = %d" %(DamID)
            layer.SetAttributeFilter(filtro)
            numDam=layer.GetFeatureCount()
            # reading first feature
            feat = layer.GetNextFeature()
            point_geom=feat.GetGeometryRef()
            Spatialref = point_geom.GetSpatialReference()
            Spatialref.AutoIdentifyEPSG()
            OriginEPSG=int(Spatialref.GetAuthorityCode(None))
            point_geom.FlattenTo2D()
            WKT=point_geom.ExportToWkt()

            if TargetEPSG!=OriginEPSG:
                targetSR = osr.SpatialReference()
                targetSR.ImportFromEPSG(TargetEPSG)
                sourceSR = osr.SpatialReference()
                sourceSR.ImportFromEPSG(OriginEPSG)
                coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
                trasformare=bool('True')
            else:
                trasformare=bool()

            if trasformare:
                geom2 = ogr.CreateGeometryFromWkt(WKT)
                geom2.Transform(coordTrans)
                wkt=geom2.ExportToWkt()
                GeomWKT="GeomFromText('%s',%d)" % (wkt,TargetEPSG)
            else:
                GeomWKT="GeomFromText('%s',%d)" % (WKT,TargetEPSG)

            # reading data fields
            FieldValue={}
            for field in FieldsList:
                found=bool()
                for  shpfieldName in  shp_field_list:
                    if  field[:5].lower() in shpfieldName.lower():
                        Value= feat.GetField(shpfieldName)
                        FieldValue[field]=Value
                        found='True'
                        break
                if not found:
                    errMsg='Shape file %s error %s field not found' % (shpfile,field)
                    NotErr= bool()
                    return NotErr, errMsg
        except:
            errMsg='Shape file %s error in DamID field' % shpfile
            NotErr= bool()
            return NotErr, errMsg

        # close the connection
        ds.Destroy()

        # controllo se esiste gia
        sql='SELECT Name FROM %s WHERE DamID=%d' % (NameTabella,DamID)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()


        if ChkDiga==None:
            insert_values=''
            sql='INSERT INTO %s (DamID'  %  (NameTabella)
            insert_values=',geom) VALUES (%d' % DamID
            for field in FieldsList:
                sql+=',%s' % field
                if 'VARCHAR' in FieldsListType[field]:
                    insert_values+=',"%s"' % FieldValue[field]
                else:
                    insert_values+=',%s' % FieldValue[field]
            sql+=insert_values
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET '  %  (NameTabella)
            for field in FieldsList:
                if 'VARCHAR' in FieldsListType[field]:
                    sql+='%s="%s",' % (field,FieldValue[field])
                else:
                    sql+='%s=%s,' % (field,FieldValue[field])
            sql+=' geom=%s' % GeomWKT
            sql+=" WHERE DamID=%d;" % DamID
        self.cur.execute(sql)

        self.conn.commit()

        self.DamID=DamID
        self.Name=FieldValue['Name']
        PathFiles=self.root_dir+os.sep+ str(self.DamID)
        self.PathFiles=os.path.abspath(PathFiles)
        if not os.path.exists(self.PathFiles):
            os.mkdir(self.PathFiles)

        return NotErr, errMsg

    def add_StudyArea(self,shpfile):

        """
        Read from a shapefile and loads the polygon of the boundary of the study area into the gdb
        """

        NotErr=bool('True')
        errMsg='OK'

        DamID=self.DamID
        NameDiga=self.Name

        # Reading shapefile path
        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        # load shapefile driver
        driver = ogr.GetDriverByName('ESRI Shapefile')

        # open the shapefile
        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg


        # reading source ayer
        layer = ds.GetLayer()

        n = layer.GetFeatureCount()

        feat = layer.GetNextFeature()

        geom=feat.GetGeometryRef()
        Spatialref = geom.GetSpatialReference()
        Spatialref.AutoIdentifyEPSG()
        OriginEPSG=int(Spatialref.GetAuthorityCode(None))

        geom.FlattenTo2D()
        TotalArea=geom.GetArea()

        # check if MULTIPOLYGON
        TipoGeom=geom.GetGeometryName()
        if TipoGeom=='POLYGON':
            multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
            multipolygon.AddGeometry(geom)
            WKT=multipolygon.ExportToWkt()
        elif geom.GetGeometryName()=='MULTIPOLYGON':
            WKT=geom.ExportToWkt()

        ds.Destroy()

        NameTabella='StudyArea'

        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NameTabella.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        if TargetEPSG!=OriginEPSG:
            targetSR = osr.SpatialReference()
            targetSR.ImportFromEPSG(TargetEPSG)
            sourceSR = osr.SpatialReference()
            sourceSR.ImportFromEPSG(OriginEPSG)
            coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
            trasformare=bool('True')
        else:
            trasformare=bool()

        if trasformare:
            geom2 = ogr.CreateGeometryFromWkt(WKT)
            geom2.Transform(coordTrans)
            wkt=geom2.ExportToWkt()
            GeomWKT="GeomFromText('%s',%d)" % (wkt,TargetEPSG)
        else:
            GeomWKT="GeomFromText('%s',%d)" % (WKT,TargetEPSG)

        sql='SELECT Name FROM %s WHERE DamID=%d' % (NameTabella,DamID)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        if ChkDiga==None:
            sql='INSERT INTO %s (DamID,Name,Area,geom) VALUES (%d,"%s",%s'  %  (NameTabella,DamID,NameDiga,TotalArea)
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET Name="%s", Area=%s'  %  (NameTabella,NameDiga,TotalArea)
            sql+=', geom=%s' % GeomWKT
            sql+=" WHERE DamID=%d;" % DamID
        self.cur.execute(sql)

        self.conn.commit()

        return NotErr, errMsg

    def add_RiverPath(self,shpfile):

        """
        Read from a shapefile the line of the river path downstream of the dam and load it into the gdb
        """

        NotErr=bool('True')
        errMsg='OK'

        DamID=self.DamID
        NameDiga=self.Name

        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        driver = ogr.GetDriverByName('ESRI Shapefile')

        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg


        layer = ds.GetLayer()

        n = layer.GetFeatureCount()

        feat = layer.GetNextFeature()

        line_sez=feat.GetGeometryRef()
        Spatialref = line_sez.GetSpatialReference()
        Spatialref.AutoIdentifyEPSG()
        OriginEPSG=int(Spatialref.GetAuthorityCode(None))

        line_sez.FlattenTo2D()
        TotalLength=line_sez.Length()

        WKT=line_sez.ExportToWkt()

        ds.Destroy()


        NameTabellaLinee='Downstreampath'

        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NameTabellaLinee.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        if TargetEPSG!=OriginEPSG:
            targetSR = osr.SpatialReference()
            targetSR.ImportFromEPSG(TargetEPSG)
            sourceSR = osr.SpatialReference()
            sourceSR.ImportFromEPSG(OriginEPSG)
            coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
            trasformare=bool('True')
        else:
            trasformare=bool()

        if trasformare:
            line = ogr.CreateGeometryFromWkt(WKT)
            line.Transform(coordTrans)
            wkt=line.ExportToWkt()
            GeomWKT="GeomFromText('%s',%d)" % (wkt,TargetEPSG)
        else:
            GeomWKT="GeomFromText('%s',%d)" % (WKT,TargetEPSG)



        sql='SELECT Name FROM %s WHERE DamID=%d' % (NameTabellaLinee,DamID)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        if ChkDiga==None:
            sql='INSERT INTO %s (DamID,Name,TotalLength,geom) VALUES (%d,"%s",%s'  %  (NameTabellaLinee,DamID,NameDiga,TotalLength)
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET Name="%s", TotalLength=%s'  %  (NameTabellaLinee,NameDiga,TotalLength)
            sql+=', geom=%s' % GeomWKT
            sql+=" WHERE DamID=%d;" % DamID
        self.cur.execute(sql)

        self.conn.commit()

        return NotErr, errMsg

    def add_MainCrossSec(self,shpfile):

        """
        Reads from a shapefile the lines of the main sections and loads them into the gdb
        """

        NotErr=bool('True')
        errMsg='OK'

        DamID=self.DamID
        NameDiga=self.Name

        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        driver = ogr.GetDriverByName('ESRI Shapefile')

        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg

        NameTabella='MainCrossSec'

        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NameTabella.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        sql='SELECT id FROM %s WHERE DamID=%d' % (NameTabella,DamID)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()
        if ChkDiga==None:
            ReplaceCrossSec=bool()
        else:
            ReplaceCrossSec=bool('True')


        layer = ds.GetLayer()
        Spatialref = layer.GetSpatialRef()
        Spatialref.AutoIdentifyEPSG()
        OriginEPSG=int(Spatialref.GetAuthorityCode(None))

        if TargetEPSG!=OriginEPSG:
            targetSR = osr.SpatialReference()
            targetSR.ImportFromEPSG(TargetEPSG)
            sourceSR = osr.SpatialReference()
            sourceSR.ImportFromEPSG(OriginEPSG)
            coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
            trasformare=bool('True')
        else:
            trasformare=bool()


        featureCount = layer.GetFeatureCount()
        if featureCount>0:
            CrossSecId_wkt={}
            count=0
            for feature in layer:
                count+=1
                try:
                    SecId=feature.GetField("id")
                    if SecId==None:
                        SecId= count
                except:
                    SecId=count
                geom = feature.GetGeometryRef()

                if trasformare:
                    WKT  = geom.ExportToWkt()
                    line = ogr.CreateGeometryFromWkt(WKT)
                    line.Transform(coordTrans)
                    wkt=line.ExportToWkt()
                    CrossSecId_wkt[SecId]=wkt
                else:
                    CrossSecId_wkt[SecId]=geom.ExportToWkt()

        ds.Destroy()

        if featureCount>0:

            if ReplaceCrossSec:
                sql='DELETE FROM %d WHERE DamID=%d;' % (NameTabella,DamID)
                self.cur.execute(sql)
                self.conn.commit()

            for ii in CrossSecId_wkt:
                wkt=CrossSecId_wkt[ii]
                GeomWKT="GeomFromText('%s',%d)" % (wkt,TargetEPSG)

                sql='INSERT INTO %s (DamID,id,edit,geom) VALUES (%d'  %  (NameTabella,DamID)
                sql+=',%d' %  ii
                sql+=',%d' %  0
                sql+=',%s' % GeomWKT
                sql+=');'
                self.cur.execute(sql)
            self.conn.commit()


        return NotErr, errMsg

    def ChkDTM(self):

        NotErr=bool('True')
        errMsg='Ok'

        if not os.path.exists(self.DTM):
            errMsg = "File DTM %s does not exists"  % self.DTM
            NotErr= bool()
            return NotErr, errMsg

        return NotErr, errMsg

    def set_DTM(self,DTMfile):


        """
        Set the path of the DTM file
        """

        NotErr=bool('True')
        errMsg='Ok'

        if not os.path.exists(DTMfile):
            errMsg = "File DTM %s does not exists"  % DTMfile
            NotErr= bool()
            return NotErr, errMsg
        else:
            self.DTM=DTMfile

        return NotErr, errMsg

    def Calc_IntermediatePoints(self):
        """
        Calculate the location, on the river path, of intermediate points between the main cross sections
        """
        NotErr=bool('True')
        errMsg='Ok'

        DamID=self.DamID
        mydb_path_user=self.myGDB_user
        fileDEM=self.DTM
        PathFiles=self.PathFiles

##        NotErr, errMsg= SetIntermediatePoints(mydb_path_user,fileDEM,DamID,DistanzaSezInterp,DeltaSezPrincipale)
        NotErr, errMsg= SetIntermediatePoints(mydb_path_user,PathFiles,fileDEM,DamID)

        return NotErr, errMsg

    def Calc_IntermediateCrossSections(self):
        """
        Calculate the location, on the river path, of intermediate cross sections between the main cross sections
        """
        NotErr=bool('True')
        errMsg='Ok'

        DamID=self.DamID
        mydb_path_user=self.myGDB_user
        fileDEM=self.DTM
        PathFiles=self.PathFiles

        NotErr, errMsg= SetCrossSec_2(mydb_path_user,PathFiles,fileDEM,DamID)

        return NotErr, errMsg

    def Calc_ValleyGeometry(self):

        """
        Divide the study area into sections between the intermediate sections and
        evaluates the geometrical and hydraulic characteristics of the various sections starting from the DTM
        """
        NotErr=bool('True')
        errMsg='Ok'

        DamID=self.DamID
        mydb_path_user=self.myGDB_user
        PathFiles=self.PathFiles
        fileDEM=self.DTM

        # grid creation Distances.tif: lateral distance from the river path
        NotErr, errMsg= GridDistanzaFiume(mydb_path_user,DamID,PathFiles,fileDEM)

        if NotErr:

            NotErr, errMsg= ModDTM(mydb_path_user,DamID,PathFiles,fileDEM)

            if NotErr:

                NotErr, errMsg= CurveAreaAltezza(mydb_path_user,DamID,PathFiles)


                if NotErr:

                    NotErr, errMsg= ParamGeomHydro(mydb_path_user,DamID,PathFiles)

                    return NotErr, errMsg

                else:
                    return NotErr, errMsg

            else:
                return NotErr, errMsg
        else:
            return NotErr, errMsg

        return NotErr, errMsg


    def DamBreakPropagation(self):

        """
        Calculate the propagation of the hypothetical dam break wave
        taking into account the geometry of the valley
        """

        NotErr=bool('True')
        errMsg='Ok'

        DamID=self.DamID
        mydb_path_user=self.myGDB_user
        PathFiles=self.PathFiles

        if os.path.exists(mydb_path_user):

            try:

                NameFileSezioni=PathFiles+os.sep+'CrossSecMean.shp'

                if os.path.exists(NameFileSezioni):

                    grafico2=0

                    NotErr, errMsg= RunCinemat(mydb_path_user,DamID,PathFiles,grafico2)

                    if not NotErr:

                        errMsg = "Error cal Dam Break"
                        NotErr= bool()
                        return NotErr, errMsg

                else:
                    errMsg = "File %s does not exists"  % NameFileSezioni
                    NotErr= bool()
                    return NotErr, errMsg

            except:

                errMsg = "Error cal Dam Break"
                NotErr= bool()
                return NotErr, errMsg
        else:

            errMsg = "File %s does not exists"  % mydb_path_user
            NotErr= bool()
            return NotErr, errMsg

        return NotErr, errMsg

    def Chk_Q_H_max(self):

        """
        Check if propagation results exists
        """
        Ok=bool()

        if os.path.exists(self.myGDB_user):

            NameTabella='Q_H_max'

            sql='SELECT PixDist FROM %s WHERE DamID=%d' % (NameTabella,self.DamID)
            self.cur.execute(sql)
            ChkDiga=cur.fetchone()

            if ChkDiga==None:
                Ok=bool()
            else:
                Ok=bool('True')

        return Ok

    def CalcFloodingArea(self,UseEnergyHead):

        """
        Calculate the map of the floodable area
        """

        NotErr=bool('True')
        errMsg='Ok'

        DamID=self.DamID
        mydb_path_user=self.myGDB_user
        PathFiles=self.PathFiles

        if os.path.exists(mydb_path_user):

            try:

                NameFileSezioni=PathFiles+os.sep+'CrossSecMean.shp'

                if os.path.exists(NameFileSezioni):

                    NotErr, errMsg = GridFloodingAreas(mydb_path_user,PathFiles,DamID,UseEnergyHead)

                    if not NotErr:

                        errMsg = "Error calc floading area "
                        NotErr= bool()
                        return NotErr, errMsg

                else:
                    errMsg = "File %s does not exists"  % NameFileSezioni
                    NotErr= bool()
                    return NotErr, errMsg

            except:

                errMsg = "Error calc floading area"
                NotErr= bool()
                return NotErr, errMsg
        else:

            errMsg = "File %s does not exists"  % mydb_path_user
            NotErr= bool()
            return NotErr, errMsg

        return NotErr, errMsg


def main():

    # -------------------------------------
    # Example of use of SimpleDambrk class
    # -------------------------------------

    # SimpleDambrk class can alternatively also be imported and used by other python modules

    mydb_path_template='..'+ os.sep + 'template'+os.sep+ 'GeoDB_template.sqlite'

    myGDB_user='..' + os.sep + 'db'+os.sep+ 'USER_GeoDB.sqlite'

    myGDB_user_dir= os.path.dirname(myGDB_user)

    if not os.path.exists(myGDB_user_dir):
        os.mkdir(myGDB_user_dir)

    if not os.path.exists(myGDB_user):
        shutil.copy (mydb_path_template, myGDB_user)

    start_time = time.time()

    MySimpleDambrk=SimpleDambrk(myGDB_user)

    # ==================
    # San Giuliano  dam
    # ==================

    fromshp=1

    if fromshp>0:

        DamID=449
        shpfile='..' + os.sep+'shp'+os.sep+'dams.shp'

        NotErr, errMsg= MySimpleDambrk.add_dam_from_shp(DamID,shpfile)

        if NotErr:
            CurrentID=MySimpleDambrk.DamID

            txt='Current DamID=%d : name=%s' %  (CurrentID,MySimpleDambrk.Name)
        else:
            txt='Upload dam data from shpafile err: %s' % errMsg
        print (txt)

    else:

        DamID=449
        Name='SAN GIULIANO'
        DamType='G'
        ResVolMcm=94.7
        Height_m=38.3
        BreachWidth=105.0
        x=1138399.4
        y=4522061.8

        # ==================
        # add dam
        # ==================
        MySimpleDambrk.add_dam(DamID,Name,DamType,ResVolMcm,Height_m,BreachWidth,x,y)

        CurrentID=MySimpleDambrk.DamID

        txt='Current DamID=%d : name=%s' %  (CurrentID,MySimpleDambrk.Name)
        print (txt)

    # ==================
    # set current dam
    # ==================
    NotErr, errMsg =  MySimpleDambrk.set_dam(DamID)
    if NotErr:
        txt='Current DamID=%d' %  (MySimpleDambrk.DamID)
    else:
        txt='Setting DamID err: %s' % errMsg
        sys.exit(txt)
    print (txt)

    # ==================
    # Set DTM
    # ==================
    DTMfile='..' + os.sep+'raster'+os.sep+'DTM_clip.tif'
    NotErr, errMsg= MySimpleDambrk.set_DTM(DTMfile)

    if NotErr:
        txt='Current DamID=%d : set DTM' %  (MySimpleDambrk.DamID)
    else:
        txt='Setting DTM err: %s' % errMsg
    print (txt)

    # ==================
    # Upload StudyArea
    # ==================
    shpfile='..' + os.sep+'shp'+os.sep+'StudyArea.shp'
    NotErr, errMsg= MySimpleDambrk.add_StudyArea(shpfile)

    if NotErr:
        txt='Current DamID=%d : added study area' %  (MySimpleDambrk.DamID)
    else:
        txt='Adding study err: %s' % errMsg
    print (txt)


    # ==================
    # Add river path
    # ==================
    shpfile='..' + os.sep+'shp'+os.sep+'River_path.shp'
    NotErr, errMsg= MySimpleDambrk.add_RiverPath(shpfile)

    if NotErr:
        txt='Current DamID=%d : added river path' %  (MySimpleDambrk.DamID)
    else:
        txt='Adding river path err: %s' % errMsg
    print (txt)


    # ======================
    # Add river MainCrossSec
    # ======================
    shpfile='..' + os.sep+'shp'+os.sep+'MainCrossSec.shp'
    NotErr, errMsg= MySimpleDambrk.add_MainCrossSec(shpfile)

    if NotErr:
        txt='Current DamID=%d : added main cross sec' %  (MySimpleDambrk.DamID)
    else:
        txt='Adding main cross sec err: %s' % errMsg
    print (txt)

    # ========================
    # Calc intermediate points
    # ========================

    NotErr, errMsg= MySimpleDambrk.Calc_IntermediatePoints()
    if NotErr:
        txt='Current DamID=%d : made intermediate points' %  (MySimpleDambrk.DamID)
    else:
        txt='Making intermediate points err: %s' % errMsg
    print (txt)

    # =================================
    # Calc intermediate cross sections
    # =================================

    NotErr, errMsg= MySimpleDambrk.Calc_IntermediateCrossSections()
    if NotErr:
        txt='Current DamID=%d : made intermediate points' %  (MySimpleDambrk.DamID)
    else:
        txt='Making intermediate points err: %s' % errMsg
    print (txt)

    # =================================
    # Calc Geometry of Valley
    # =================================

    NotErr, errMsg= MySimpleDambrk.Calc_ValleyGeometry()
    if NotErr:
        txt='Current DamID=%d : made Geometry of Valley' %  (MySimpleDambrk.DamID)
    else:
        txt='Making Geometry of Valley err: %s' % errMsg
    print (txt)

    # =================================
    # Calc DamBreak propagation
    # =================================

    NotErr, errMsg= MySimpleDambrk.DamBreakPropagation()
    if NotErr:
        txt='Current DamID=%d : calc Dam Break propagation' %  (MySimpleDambrk.DamID)
    else:
        txt='Making cal Dam Break propagation err: %s' % errMsg
    print (txt)

    # =================================
    # Calc Floading area
    # =================================

    NotErr, errMsg= MySimpleDambrk.CalcFloodingArea()
    if NotErr:
        txt='Current DamID=%d : calc floading area' %  (MySimpleDambrk.DamID)
    else:
        txt='Making cal calc floading area err: %s' % errMsg
    print (txt)


    elapsed_time = time.time() - start_time
    print ('elapsed time=  %s sec' % elapsed_time)

if __name__ == '__main__':
    main()
