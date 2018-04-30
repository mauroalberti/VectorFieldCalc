# -*- coding: utf-8 -*-


from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from qgis.core import *

import resources

from vfc_dialog import VfcDialog


plugin_name = "VectorFieldCalc"


class VectorFieldCal:

    def __init__(self, iface):        
        """
        Save reference to the QGIS interface
        
        :param iface: 
        """
        
        self.iface = iface

    def initGui(self):
        """
        Create action that will start plugin configuration
        
        :return: 
        """
        
        self.action = QAction(
            QIcon(":/plugins/VectorFieldCalc/icon.png"),
            "VectorFieldCalc", 
            self.iface.mainWindow())
        
        # connect the action to the run method
        self.action.triggered.connect(self.run)

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&VectorFieldCalc", self.action)

    def unload(self):
        """
        Remove the plugin menu item and icon

        :return:
        """

        self.iface.removePluginMenu("&VectorFieldCalc",self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        """
        Run the plugin

        :return:
        """

        # create the dialog        
        dlg = VfcDialog(plugin_name)
 
        # show the dialog
        dlg.show()        
        dlg.exec_() 
        
        
        
               
 