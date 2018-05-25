"""

/***************************************************************************
 VectorFieldCalc
                                 A QGIS plugin
 Processing of vector fields
                             -------------------
        begin                : 2011-11-27
        copyright            : (C) 2011-2018 by Mauro Alberti
        email                : alberti.m65@gmail.com
 ***************************************************************************/


/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

"""
from __future__ import absolute_import


def classFactory(iface):
    
    from .vfc_main import VectorFieldCal
    return VectorFieldCal(iface)


