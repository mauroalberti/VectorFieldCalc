"""

/***************************************************************************
 VectorFieldCalc
                                 A QGIS plugin
 Processing of vector fields
                             -------------------
        begin                : 2011-11-27
        current version      : 2012-09-16
        copyright            : (C) 2011-2012 by Mauro Alberti - www.malg.eu
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


def name():
    
    return "VectorFieldCalc"


def description():
    
    return "Calculations of vector field parameters (divergence, magnitude, flowlines, etc.)"


def version():
    
    return "1.1"


def icon():
    
    return "icon.png"


def qgisMinimumVersion():
    
    return "1.0"


def classFactory(iface):
    
    from vfc_main import VectorFieldCal
    return VectorFieldCal(iface)


