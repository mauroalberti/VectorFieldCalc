

from builtins import str
from qgis.PyQt.QtCore import *

from qgis.core import *


# from module RASTERCALC by Barry Rowlingson

def lastUsedDir():
  
  settings = QSettings()
  return settings.value(
    "/VectorFieldCalc/lastDir",
    str(""),
    type=str)

def setLastUsedDir(lastDir):
  
  path = QFileInfo(lastDir).absolutePath()
  settings = QSettings()
  settings.setValue(
    "/VectorFieldCalc/lastDir",
    str(path))

  

  

