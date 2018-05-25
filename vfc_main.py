from __future__ import absolute_import

from builtins import object
import numpy as np

from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import *
from qgis.PyQt.QtWidgets import *

from qgis.core import *

import ogr

#from . import resources
from .vfc_dialog import vfc_dialog
from .vfc_classes import *



class VectorFieldCal(object):

    def __init__(self, iface):        
        # Save reference to the QGIS interface
        
        self.iface = iface

    def initGui(self):
        # Create action that will start plugin configuration
        
        self.action = QAction(QIcon(":/plugins/VectorFieldCalc/icon.png"), \
            "VectorFieldCalc", self.iface.mainWindow())
        
        # connect the action to the run method
        self.action.triggered.connect( self.run )

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&VectorFieldCalc", self.action)

    def unload(self):
        # Remove the plugin menu item and icon
        
        self.iface.removePluginMenu("&VectorFieldCalc", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):

        # create the dialog        
        dlg = vfc_dialog() 
 
        # show the dialog
        dlg.show()        
        dlg.exec_() 
