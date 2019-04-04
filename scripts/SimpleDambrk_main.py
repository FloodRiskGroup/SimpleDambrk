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

        self.ID_Diga=ID_Diga=-1
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


    def set_dam(self,ID_Diga):

        """
        Definisce la diga corrente di cui effettuare i calcoli
        """

        self.ID_Diga=ID_Diga

        NomeTabella='DAMS'

        sql='SELECT Nome FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()
        self.Name=ChkDiga[0]

        self.PathFiles=self.root_dir+os.sep+ str(self.ID_Diga)


    def add_dam(self,ID_Diga,Nome,TipoDiga,Volume_mlnm3,Altezza_m,Breccia_m,x,y):

        """
        Aggiunge i dati di una nuova diga
        """

        NomeTabella='DAMS'

        # Controllo sistema riferimento della linea rispetto ai Catchments
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        # controllo se esiste gia
        sql='SELECT Nome FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        pt=ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(x, y)
        pt.FlattenTo2D()
        WKT_Point=pt.ExportToWkt()

        GeomWKT="GeomFromText('%s',%d)" % (WKT_Point,TargetEPSG)

        if ChkDiga==None:
            sql='INSERT INTO %s (ID_Diga,Nome,TipoDiga,Volume_mlnm3,Altezza_m,Breccia_m,geom) VALUES (%d,"%s","%s",%s,%s,%s'  %  (NomeTabella,ID_Diga,Nome,TipoDiga,Volume_mlnm3,Altezza_m,Breccia_m)
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET Nome="%s", TipoDiga="%s", Volume_mlnm3=%s, Altezza_m=%s, Breccia_m=%s'  %  (NomeTabella,Nome,TipoDiga,Volume_mlnm3,Altezza_m,Breccia_m)
            sql+=', geom=%s' % GeomWKT
            sql+=" WHERE ID_Diga=%d;" % ID_Diga
        self.cur.execute(sql)

        self.conn.commit()

        self.ID_Diga=ID_Diga
        self.Name=Nome
        self.PathFiles=self.root_dir+os.sep+ str(self.ID_Diga)


    def add_StudyArea(self,shpfile):

        """
        Legge da uno shapefile e carica nel gdb il poligono del contorno dell'area di studio
        """

        NotErr=bool('True')
        errMsg='OK'

        ID_Diga=self.ID_Diga
        NomeDiga=self.Name

        # lettura shapefile
        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        # carico il driver dello shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')

        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg


        # leggo il layer dalla sorgente dei dati
        layer = ds.GetLayer()

        # conto il numero di Feature
        n = layer.GetFeatureCount()

        # leggo la prima feature
        feat = layer.GetNextFeature()

        geom=feat.GetGeometryRef()
        Spatialref = geom.GetSpatialReference()
        Spatialref.AutoIdentifyEPSG()
        OriginEPSG=int(Spatialref.GetAuthorityCode(None))

        geom.FlattenTo2D()
        TotalArea=geom.GetArea()

        # controllo di avere un MULTIPOLIGON
        TipoGeom=geom.GetGeometryName()
        if TipoGeom=='POLYGON':
            multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
            multipolygon.AddGeometry(geom)
            WKT=multipolygon.ExportToWkt()
        elif geom.GetGeometryName()=='MULTIPOLYGON':
            WKT=geom.ExportToWkt()

        # chiudo la connessione
        ds.Destroy()

        NomeTabella='PoligonoValleDiga'

        # codice del sistema di riferimento della tabella
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
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

        # controllo se esiste gia
        sql='SELECT Nome FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        if ChkDiga==None:
            sql='INSERT INTO %s (ID_Diga,Nome,Area,geom) VALUES (%d,"%s",%s'  %  (NomeTabella,ID_Diga,NomeDiga,TotalArea)
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET Nome="%s", Area=%s'  %  (NomeTabella,NomeDiga,TotalArea)
            sql+=', geom=%s' % GeomWKT
            sql+=" WHERE ID_Diga=%d;" % ID_Diga
        self.cur.execute(sql)

        self.conn.commit()

        return NotErr, errMsg

    def add_RiverPath(self,shpfile):

        """
        Legge da uno shapefile e carica nel gdb la linea dell'asse del fiume a valle della diga
        """

        NotErr=bool('True')
        errMsg='OK'

        ID_Diga=self.ID_Diga
        NomeDiga=self.Name

        # lettura shapefile
        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        # carico il driver dello shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')

        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg


        # leggo il layer dalla sorgente dei dati
        layer = ds.GetLayer()

        # conto il numero di Feature
        n = layer.GetFeatureCount()

        # leggo la prima feature
        feat = layer.GetNextFeature()

        line_sez=feat.GetGeometryRef()
        Spatialref = line_sez.GetSpatialReference()
        Spatialref.AutoIdentifyEPSG()
        OriginEPSG=int(Spatialref.GetAuthorityCode(None))

        line_sez.FlattenTo2D()
        TotalLength=line_sez.Length()

        WKT=line_sez.ExportToWkt()

        # chiudo la connessione
        ds.Destroy()


        NomeTabellaLinee='LineaValleDiga'

        # codice del sistema di riferimento della tabella
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabellaLinee.lower())
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



        # controllo se esiste gia la linea a valle della diga in oggetto
        sql='SELECT Nome FROM %s WHERE ID_Diga=%d' % (NomeTabellaLinee,ID_Diga)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()

        if ChkDiga==None:
            sql='INSERT INTO %s (ID_Diga,Nome,TotalLength,geom) VALUES (%d,"%s",%s'  %  (NomeTabellaLinee,ID_Diga,NomeDiga,TotalLength)
            sql+=', %s' % GeomWKT
            sql+=");"
        else:
            sql='UPDATE %s SET Nome="%s", TotalLength=%s'  %  (NomeTabellaLinee,NomeDiga,TotalLength)
            sql+=', geom=%s' % GeomWKT
            sql+=" WHERE ID_Diga=%d;" % ID_Diga
        self.cur.execute(sql)

        self.conn.commit()

        return NotErr, errMsg

    def add_MainCrossSec(self,shpfile):

        """
        Legge da uno shapefile e carica nel gdb le linee della traccia delle sezioni principali
        """

        NotErr=bool('True')
        errMsg='OK'

        ID_Diga=self.ID_Diga
        NomeDiga=self.Name

        # lettura shapefile
        if not os.path.exists(shpfile):
            errMsg = "File %s does not exists"  % shpfile
            NotErr= bool()
            return NotErr, errMsg

        # carico il driver dello shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')

        ds = driver.Open(shpfile, 0)
        if ds is None:
            errMsg='Could not open file %s' % shpfile
            NotErr= bool()
            return NotErr, errMsg

        NomeTabella='MainCrossSec'

        # codice del sistema di riferimento della tabella
        sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (NomeTabella.lower())
        self.cur.execute(sql)
        record=self.cur.fetchone()
        if record!=None:
            TargetEPSG=record[0]
        else:
            TargetEPSG=32632

        # controllo se esistono gia la linea a valle della diga in oggetto
        sql='SELECT id FROM %s WHERE ID_Diga=%d' % (NomeTabella,ID_Diga)
        self.cur.execute(sql)
        ChkDiga=self.cur.fetchone()
        if ChkDiga==None:
            ReplaceCrossSec=bool()
        else:
            ReplaceCrossSec=bool('True')


        # leggo il layer dalla sorgente dei dati
        # --------------------------------------
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
                sql='DELETE FROM %d WHERE ID_Diga=%d;' % (NomeTabella,ID_Diga)
                self.cur.execute(sql)
                self.conn.commit()

            for ii in CrossSecId_wkt:
                wkt=CrossSecId_wkt[ii]
                GeomWKT="GeomFromText('%s',%d)" % (wkt,TargetEPSG)

                sql='INSERT INTO %s (ID_Diga,id,edit,geom) VALUES (%d'  %  (NomeTabella,ID_Diga)
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
        Assegna il path del file del DTM
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
        Effettua il tracciamento, sull'asse del fiume, dei punti intermedi alle sezioni principali
        """
        NotErr=bool('True')
        errMsg='Ok'

        ID_Diga=self.ID_Diga
        mydb_path_user=self.myGDB_user
        fileDEM=self.DTM
        PathFiles=self.PathFiles

##        NotErr, errMsg= SetIntermediatePoints(mydb_path_user,ID_Diga,fileDEM,DistanzaSezInterp,DeltaSezPrincipale)
        NotErr, errMsg= SetIntermediatePoints(mydb_path_user,ID_Diga,PathFiles,fileDEM)

        return NotErr, errMsg

    def Calc_IntermediateCrossSections(self):
        """
        Effettua il tracciamento delle sezioni intermedie mediante interpolazione fra quelle pricipali
        """
        NotErr=bool('True')
        errMsg='Ok'

        ID_Diga=self.ID_Diga
        mydb_path_user=self.myGDB_user
        fileDEM=self.DTM
        PathFiles=self.PathFiles

        NotErr, errMsg= SetCrossSec_2(mydb_path_user,ID_Diga,PathFiles,fileDEM)

        return NotErr, errMsg

    def Calc_ValleyGeometry(self):

        """
        Divide la zona di studio in tratti compresi fra le sezioni intermedie e
        valuta le caratteristiche geometriche idrauliche dei vari tratti a partire dal DTM
        """
        NotErr=bool('True')
        errMsg='Ok'

        ID_Diga=self.ID_Diga
        mydb_path_user=self.myGDB_user
        PathFiles=self.PathFiles
        fileDEM=self.DTM

        # creazione grid Distances.tif: distanza laterale dall'asse del fiume
        NotErr, errMsg= GridDistanzaFiume(mydb_path_user,ID_Diga,PathFiles,fileDEM)

        if NotErr:
            # modifica delle quote del DTM affinche' ci sia sempre una pendenza
            # minima verso l'asse del fiume. Risultato: StreamDHFilled.tif

            NotErr, errMsg= ModDTM(mydb_path_user,ID_Diga,PathFiles,fileDEM)

            if NotErr:
                # CreaCurveAreaAltezza: effettua il conteggio, per ogni tratto,
                # e per ogni dh con passo 1 metro, del numero di celle
                # soggiacenti la differenza di altezza dh dal fiume
                # risultato: MatricePixel.csv

                NotErr, errMsg= CurveAreaAltezza(mydb_path_user,ID_Diga,PathFiles)


                if NotErr:
                    # ParametriGeomIdraulici: legge  MatricePixel.csv e crea
                    # MatriceAexp.csv: file con le curve geometriche ed idrauliche

                    NotErr, errMsg= ParamGeomHydro(mydb_path_user,ID_Diga,PathFiles)

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
        Calcola la propagazione dell'onda di ipotetico dam break
        tenendo conto della geometria della valle
        """

        NotErr=bool('True')
        errMsg='Ok'

        ID_Diga=self.ID_Diga
        mydb_path_user=self.myGDB_user
        PathFiles=self.PathFiles

        if os.path.exists(mydb_path_user):

            try:

                NomeFileSezioni=PathFiles+os.sep+'CrossSecMean.shp'

                if os.path.exists(NomeFileSezioni):

                    grafico2=0
                    # calcolo propagazione
                    NotErr, errMsg= RunCinemat(mydb_path_user,ID_Diga,PathFiles,grafico2)

                    if not NotErr:

                        errMsg = "Error cal Dam Break"
                        NotErr= bool()
                        return NotErr, errMsg

                else:
                    errMsg = "File %s does not exists"  % NomeFileSezioni
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
        Controlla l'esistenza risultati della propagazione
        """
        Ok=bool()

        if os.path.exists(self.myGDB_user):

            NomeTabella='Q_H_max'

            # controllo se esistono i risultati della propagazione
            sql='SELECT PixDist FROM %s WHERE ID_Diga=%d' % (NomeTabella,self.ID_Diga)
            self.cur.execute(sql)
            ChkDiga=cur.fetchone()

            if ChkDiga==None:
                Ok=bool()
            else:
                Ok=bool('True')

        return Ok

    def CalcFloodArea(self):

        """
        Calcola la mappa dell'area inondabile
        """
        # Calcola la propagazione dell'onda per la seconda volta
        # tenendo conto della geometria delle sezioni tratta dal DTM
        ID_Diga=self.ID_Diga
        mydb_path_user=self.myGDB_user
        PathFiles=self.PathFiles

        if os.path.exists(mydb_path_user):

            try:

                NomeFileSezioni=PathFiles+os.sep+'CrossSecMean.shp'

                if os.path.exists(NomeFileSezioni):


                    TargetTabella='AreaInondabileValleDiga'

                    # codice del sistema di riferimento della tabella
                    sql="SELECT srid FROM geometry_columns WHERE f_table_name='%s'" % (TargetTabella.lower())
                    cur.execute(sql)
                    record=cur.fetchone()
                    if record!=None:
                        OriginEPSG=record[0]
                    else:
                        OriginEPSG=32632

                    # calcolo area inondabile
                    NotErr, errMsg = GridAreeInondabili(mydb_path_user,ID_Diga)


                    if NotErr:

                        QMessageBox.information(None, "SimpleDambrk", self.tr("Eseguito calcolo area inondabile"))

                        self.ControlloTab3()
                        self.ControlloTab4()

##                        # effettuo il salvataggio
##                        sql='SELECT PKUID,id FROM %s WHERE ID_Diga=%d' % (TargetTabella,ID_Diga)
##                        cur.execute(sql)
##                        ListaTratti=cur.fetchall()
##                        if len(ListaTratti)>0:
##                            # cancello ed aggiungo
##                            sql='DELETE FROM %s WHERE ID_Diga=%d' % (TargetTabella,ID_Diga)
##                            cur.execute(sql)
##                            conn.commit()
##                        # aggiungo
##                        for rec in MatriceRisultati:
##                            sql='INSERT INTO %s (ID_Diga,id,DV,FloodSeverity,WarningTimeMin,FatalityRate,geom) VALUES (%d'  %  (TargetTabella,ID_Diga)
##                            sql+=',%d,%.2f,"%s",%d,%.3f' % (rec[0],rec[1],rec[2],rec[3],rec[4])
##                            poly=ogr.CreateGeometryFromWkt(rec[5])
##                            area=poly.GetArea()
##
##                            # controllo di avere un MULTIPOLIGON
##                            TipoGeom=poly.GetGeometryName()
##                            if TipoGeom=='POLYGON':
##                                multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
##                                multipolygon.AddGeometry(poly)
##                                wkt2=multipolygon.ExportToWkt()
##                            elif poly.GetGeometryName()=='MULTIPOLYGON':
##                                wkt2=poly.ExportToWkt()
##                            poly2=ogr.CreateGeometryFromWkt(wkt2)
##                            area2=poly2.GetArea()
##
##                            GeomWKT="GeomFromText('%s',%d)" % (wkt2,OriginEPSG)
##
##                            sql+=',%s' % GeomWKT
##                            sql+=');'
##                            cur.execute(sql)
##                            # prova
##                            if rec[0]==1:
##                                sql='INSERT INTO AreaInondabileValleDigaTmp (ID_Diga,geom) VALUES (%d,%s);'  %  (ID_Diga,GeomWKT)
##                                cur.execute(sql)
##                                conn.commit()
##
##
##                        conn.commit()

##                        QMessageBox.information(None, "SimpleDambrk", self.tr("Eseguito Salvataggio Area Inondabile"))
                        self.pushButton_ZoomAreaInondabile.setEnabled(True)

                    else:
                        QMessageBox.information(None, "SimpleDambrk", errMsg)
                        self.pushButton_ZoomAreaInondabile.setEnabled(False)


                    # Close communication with the database
                    cur.close()
                    conn.close()
                else:

                    errMsg =self.tr("Attenzione manca CrossSecMean.shp : effettuare prima il calcolo delle sezioni a valle !")


            except:
                pass


def main():

    mydb_path_template='..'+ os.sep + 'template'+os.sep+ 'GeoDB_template.sqlite'

    myGDB_user='..' + os.sep + 'db'+os.sep+ 'USER_GeoDB.sqlite'

    if not os.path.exists(myGDB_user):
        shutil.copy (mydb_path_template, myGDB_user)


    MySimpleDambrk=SimpleDambrk(myGDB_user)

    # ==================
    # San Giuliano  dam
    # ==================
    ID_Diga=449
    Nome='SAN GIULIANO'
    TipoDiga='G'
    Volume_mlnm3=94.7
    Altezza_m=38.3
    Breccia_m=105.0
    x=1138399.4
    y=4522061.8

    # ==================
    # add dam
    # ==================
    MySimpleDambrk.add_dam(ID_Diga,Nome,TipoDiga,Volume_mlnm3,Altezza_m,Breccia_m,x,y)

    CurrentID=MySimpleDambrk.ID_Diga

    txt='Current ID_Diga=%d : name=%s' %  (CurrentID,MySimpleDambrk.Name)
    print (txt)

    # ==================
    # set dam
    # ==================
    MySimpleDambrk.set_dam(ID_Diga)

    # ==================
    # Set DTM
    # ==================
    DTMfile='..' + os.sep+'data'+os.sep+'DTM_clip.tif'
    NotErr, errMsg= MySimpleDambrk.set_DTM(DTMfile)

    if NotErr:
        txt='Current ID_Diga=%d : set DTM' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Setting DTM err: %s' % errMsg
    print (txt)

    # ==================
    # Upload StudyArea
    # ==================
    shpfile='..' + os.sep+'shp'+os.sep+'StudyArea.shp'
    NotErr, errMsg= MySimpleDambrk.add_StudyArea(shpfile)

    if NotErr:
        txt='Current ID_Diga=%d : added study area' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Adding study err: %s' % errMsg
    print (txt)


    # ==================
    # Add river path
    # ==================
    shpfile='..' + os.sep+'shp'+os.sep+'River_path.shp'
    NotErr, errMsg= MySimpleDambrk.add_RiverPath(shpfile)

    if NotErr:
        txt='Current ID_Diga=%d : added river path' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Adding river path err: %s' % errMsg
    print (txt)


    # ======================
    # Add river MainCrossSec
    # ======================
    shpfile='..' + os.sep+'shp'+os.sep+'MainCrossSec.shp'
    NotErr, errMsg= MySimpleDambrk.add_MainCrossSec(shpfile)

    if NotErr:
        txt='Current ID_Diga=%d : added main cross sec' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Adding main cross sec err: %s' % errMsg
    print (txt)

    # ========================
    # Calc intermediate points
    # ========================

    NotErr, errMsg= MySimpleDambrk.Calc_IntermediatePoints()
    if NotErr:
        txt='Current ID_Diga=%d : made intermediate points' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Making intermediate points err: %s' % errMsg
    print (txt)

    # =================================
    # Calc intermediate cross sections
    # =================================

    NotErr, errMsg= MySimpleDambrk.Calc_IntermediateCrossSections()
    if NotErr:
        txt='Current ID_Diga=%d : made intermediate points' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Making intermediate points err: %s' % errMsg
    print (txt)

    # =================================
    # Calc Geometry of Valley
    # =================================

    NotErr, errMsg= MySimpleDambrk.Calc_ValleyGeometry()
    if NotErr:
        txt='Current ID_Diga=%d : made Geometry of Valley' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Making Geometry of Valley err: %s' % errMsg
    print (txt)

    # =================================
    # Calc DamBreak propagation
    # =================================

    NotErr, errMsg= MySimpleDambrk.DamBreakPropagation()
    if NotErr:
        txt='Current ID_Diga=%d : calc Dam Break propagation' %  (MySimpleDambrk.ID_Diga)
    else:
        txt='Making cal Dam Break propagation err: %s' % errMsg
    print (txt)


if __name__ == '__main__':
    main()
