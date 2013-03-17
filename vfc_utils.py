

from PyQt4.QtCore import *

from qgis.core import *




# from module RASTERCALC by Barry Rowlingson

def lastUsedDir():
  settings = QSettings()
  return settings.value( "/VectorFieldCalc/lastDir", QVariant( "" ) ).toString()

def setLastUsedDir( lastDir ):
  path = QFileInfo( lastDir ).absolutePath()
  settings = QSettings()
  settings.setValue( "/VectorFieldCalc/lastDir", QVariant( path ) )
  

  
  