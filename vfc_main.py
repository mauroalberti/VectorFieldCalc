from __future__ import absolute_import

from builtins import object

import os

import numpy as np


import webbrowser

from qgis.PyQt.QtCore import *


from qgis.PyQt.QtGui import *
from qgis.PyQt.QtWidgets import *

from qgis.core import *


from . import resources

from .qgis_utils.gui import *


from .main_dialog import MainDialog, HelpDialog


_plugin_name_ = "VectorFieldCalc"


class VectorFieldCal(object):

    def __init__(self, iface):        
        # Save reference to the QGIS interface
        
        self.iface = iface
        self.plugin_name = _plugin_name_
        self.actions = []

    def initGui(self):
        # Create action that will start plugin configuration

        self.qactMain = create_action(
            ':/plugins/{}/icons/icon.png'.format(self.plugin_name),
            'Open VectorFieldCalc',
            self.run,
            whats_this="VectorFieldCalc",
            parent=self.iface.mainWindow())
        self.iface.addPluginToMenu(self.plugin_name,
                                   self.qactMain)

        self.qactOpenHelp = create_action(
            ':/plugins/{}/icons/help.ico'.format(self.plugin_name),
            'Help',
            self.open_html_help,
            whats_this="Topographic and geological profiles Help",
            parent=self.iface.mainWindow())
        self.iface.addPluginToMenu(self.plugin_name,
                                   self.qactOpenHelp)

    def unload(self):
        # Remove the plugin menu item and icon

        self.iface.removePluginMenu(self.plugin_name, self.qactMain)
        self.iface.removePluginMenu(self.plugin_name, self.qactOpenHelp)

    def run(self):

        # create the dialog        
        dlg = MainDialog()
 
        # show the dialog
        dlg.show()        
        dlg.exec_()

    def open_html_help(self):

        file_path = '{}/help/help.html'.format(os.path.dirname(__file__))

        # create the dialog
        dlg = HelpDialog(file_path)

        # show the dialog
        dlg.show()
        dlg.exec_()