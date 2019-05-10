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

 The script made a simplified study of the dambreak of a dam

"""
from SimpleDambrk_main import SimpleDambrk
import os, sys
import time


def run_script(iface):

    Script_dir=os.path.dirname(os.path.abspath(__file__))
    os.chdir(Script_dir)

    # name of the user geodatabse
    myGDB_user='..' + os.sep + 'db'+os.sep+ 'USER_GeoDB.sqlite'

    if not os.path.exists(myGDB_user):
        txt='error %s not found'  %  myGDB_user
        sys.exit(txt)

    start_time = time.time()

    # creating a new, unique instance of the class SimpleDambrk
    MySimpleDambrk=SimpleDambrk(myGDB_user)


    # --------------------------
    # =========================
    # Start setting current dam
    # =========================
    # --------------------------


    # =========================
    # set current dam
    # =========================
    # unique ID of the dam into the geodatabase
    DamID=0

    NotErr, errMsg =  MySimpleDambrk.set_dam(DamID)
    if NotErr:
        txt='Current DamID=%d' %  (MySimpleDambrk.DamID)
    else:
        txt='Setting DamID err: %s' % errMsg
        sys.exit(txt)
    print (txt)

    # ==================
    # Set current DTM
    # ==================

    # pathname of the DTM
    DTMfile='..' + os.sep+'Gleno_data'+os.sep+'Gleno_ClipDTM.tif'
    NotErr, errMsg= MySimpleDambrk.set_DTM(DTMfile)

    if NotErr:
        txt='Current DamID=%d : set DTM' %  (MySimpleDambrk.DamID)
    else:
        txt='Setting DTM err: %s' % errMsg
    print (txt)

    # --------------------------
    # =========================
    # Start of calculation
    # =========================
    # --------------------------

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
        txt='Current DamID=%d : made intermediate CrossSections' %  (MySimpleDambrk.DamID)
    else:
        txt='Making intermediate points err: %s' % errMsg
    print (txt)

    # =================================
    # Calc Geometry of Valley
    # =================================
    txt='Start calculating Geometry of Valley'
    print (txt)

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
        txt='Current DamID=%d : Dam Break propagation - done' %  (MySimpleDambrk.DamID)
    else:
        txt='Making cal Dam Break propagation err: %s' % errMsg
    print (txt)

    # =================================
    # Calc Floading area
    # =================================
##    UseEnergyHead= bool('True')
    UseEnergyHead= bool()
    NotErr, errMsg= MySimpleDambrk.CalcFloodingArea(UseEnergyHead)
    if NotErr:
        txt='Current DamID=%d : Flood Area Extent cal - done' %  (MySimpleDambrk.DamID)
    else:
        txt='Making calc floading area err: %s' % errMsg
    print (txt)


    elapsed_time = time.time() - start_time
    print ('elapsed time=  %s sec' % elapsed_time)

if __name__ == '__main__':

    iface=''
    run_script(iface)
