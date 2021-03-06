# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SimpleDamBrk - upload data

        begin                : 2019-05-02
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

 The script loads, into the user geodatabase, the data necessary for
 the simplified study of the dambreak of a dam

"""
from SimpleDambrk_main import SimpleDambrk
import os
import shutil


def run_script(iface):

    Script_dir=os.path.dirname(os.path.abspath(__file__))
    os.chdir(Script_dir)

    # path of the geodatabase template, which is used the first time to
    # create a new user geodatabase that has the schema required by the tool
    mydb_path_template='..'+ os.sep + 'template'+os.sep+ 'GeoDB_template.sqlite'

    # name of the user geodatabse
    myGDB_user='..' + os.sep + 'db'+os.sep+ 'USER_GeoDB.sqlite'

    myGDB_user_dir= os.path.dirname(myGDB_user)

    if not os.path.exists(myGDB_user_dir):
        os.mkdir(myGDB_user_dir)

    if not os.path.exists(myGDB_user):
        shutil.copy (mydb_path_template, myGDB_user)

    # creating a new, unique instance of the class SimpleDambrk
    MySimpleDambrk=SimpleDambrk(myGDB_user)


    # ========================
    # Gleno  dam  data
    # =======================

    DamID=0
    shpfile='..' + os.sep+'Gleno_data'+os.sep+'Gleno_dam.shp'

    NotErr, errMsg= MySimpleDambrk.add_dam_from_shp(DamID,shpfile)

    if NotErr:
        CurrentID=MySimpleDambrk.DamID

        txt='Current DamID=%d : name=%s' %  (CurrentID,MySimpleDambrk.Name)
    else:
        txt='Upload dam data from shpafile err: %s' % errMsg
    print (txt)



    # ===============================================
    # Upload StudyArea polygon downstream of the dam
    # ===============================================
    # polygon shapefile of the polygon that defines the outline of the study area
    shpfile='..' + os.sep+'Gleno_data'+os.sep+'Gleno_StudyArea.shp'
    NotErr, errMsg= MySimpleDambrk.add_StudyArea(shpfile)

    if NotErr:
        txt='Current DamID=%d : added study area' %  (MySimpleDambrk.DamID)
    else:
        txt='Adding study err: %s' % errMsg
    print (txt)


    # =====================================
    # Add river path downstream of the dam
    # =====================================
    # linear shape file of the river downstream the dam
    # the line must be digitized from upstream to downstream
    shpfile='..' + os.sep+'Gleno_data'+os.sep+'Gleno_RiverPath.shp'
    NotErr, errMsg= MySimpleDambrk.add_RiverPath(shpfile)

    if NotErr:
        txt='Current DamID=%d : added river path' %  (MySimpleDambrk.DamID)
    else:
        txt='Adding river path err: %s' % errMsg
    print (txt)


    # ======================
    # Add river MainCrossSec
    # ======================
    # linear shapefile of the main river cross section downstream the dam
    # the lines of the sections in the shapefile must be inserted from upstream to downstream
    # each line must join two points, the first on the left bank and the second on the right bank
    shpfile='..' + os.sep+'Gleno_data'+os.sep+'Gleno_MainCrossSec.shp'
    NotErr, errMsg= MySimpleDambrk.add_MainCrossSec(shpfile)

    if NotErr:
        txt='Current DamID=%d : added main cross sec' %  (MySimpleDambrk.DamID)
    else:
        txt='Adding main cross sec err: %s' % errMsg
    print (txt)

if __name__ == '__main__':

    iface=''
    run_script(iface)
