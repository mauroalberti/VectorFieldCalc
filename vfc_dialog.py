from __future__ import absolute_import


from builtins import str
from builtins import zip
import os, webbrowser
from osgeo import ogr

from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import *
from qgis.PyQt.QtWidgets import *

from qgis.core import *

from .vfc_utils import *
from .vfc_classes import *


class vfc_dialog(QDialog):
    
    def __init__(self):
        
        QDialog.__init__(self)
        self.setup_ui()        
 
    def setup_ui(self):

        self.setWindowTitle("VectorFieldCalc")
        
        self.dialog_layout = QVBoxLayout()
        self.setLayout(self.dialog_layout)
        
        self.input_raster_widget = self.setup_input_raster_widget()        
        self.dialog_layout.addWidget(self.input_raster_widget)
        
        self.calculations_widget = self.setup_calculations_widget()        
        self.dialog_layout.addWidget(self.calculations_widget)
        self.calculations_widget.setCurrentIndex(0)

        self.manage_widget = self.setup_manage_widget()
        self.dialog_layout.addWidget(self.manage_widget)         

        self.adjustSize()
 
    def setup_input_raster_widget(self):
        
        inraster_groupbox = QGroupBox("Input rasters")        
        inraster_layout = QFormLayout()

        self.inraster_x_comboBox = QComboBox()
        inraster_layout.addRow(QLabel("X-axis components"), self.inraster_x_comboBox)
        self.inraster_y_comboBox = QComboBox()
        inraster_layout.addRow(QLabel("Y-axis components"), self.inraster_y_comboBox)
                
        # append loaded raster layers to combo boxes
        #self.layermap = QgsMapLayerRegistry.instance().mapLayers()
        self.layermap = QgsProject.instance().mapLayers()
        for (name, layer) in list(self.layermap.items()):
            if layer.type() == QgsMapLayer.RasterLayer: 
                self.inraster_x_comboBox.addItem(layer.name())
                self.inraster_y_comboBox.addItem(layer.name()) 
 
        inraster_groupbox.setLayout(inraster_layout)
               
        return inraster_groupbox

    def setup_vector_operators_widget(self):

        vectorfieldoperators_tab = QWidget()
        vectorfieldoperators_layout = QGridLayout()
        
        self.magnitude_calc_choice_checkBox = QCheckBox("Magnitude")
        vectorfieldoperators_layout.addWidget(self.magnitude_calc_choice_checkBox, 0, 0, 1, 1)
        
        self.magnitude_outraster_lineEdit = QLineEdit()
        self.magnitude_outraster_lineEdit.setPlaceholderText("magn.asc")
        vectorfieldoperators_layout.addWidget(self.magnitude_outraster_lineEdit, 0, 1, 1, 1)
        
        self.magnitude_set_outraster = QPushButton("...")
        self.magnitude_set_outraster.clicked.connect(self.select_output_rasterFile)
        vectorfieldoperators_layout.addWidget(self.magnitude_set_outraster, 0, 2, 1, 1)
        
        self.orientations_calc_choice_checkBox = QCheckBox("Orientations")
        vectorfieldoperators_layout.addWidget(self.orientations_calc_choice_checkBox, 1, 0, 1, 1)
        
        self.orientations_outraster_lineEdit = QLineEdit()
        self.orientations_outraster_lineEdit.setPlaceholderText("orient.asc")
        vectorfieldoperators_layout.addWidget(self.orientations_outraster_lineEdit, 1, 1, 1, 1)
        
        self.orientations_set_outraster = QPushButton("...")
        self.orientations_set_outraster.clicked.connect(self.select_output_rasterFile)
        vectorfieldoperators_layout.addWidget(self.orientations_set_outraster, 1, 2, 1, 1)
        
        self.divergence_calc_choice_checkBox = QCheckBox("Divergence")
        vectorfieldoperators_layout.addWidget(self.divergence_calc_choice_checkBox, 2, 0, 1, 1)
        
        self.divergence_outraster_lineEdit = QLineEdit()
        self.divergence_outraster_lineEdit.setPlaceholderText("div.asc")
        vectorfieldoperators_layout.addWidget(self.divergence_outraster_lineEdit, 2, 1, 1, 1)
        
        self.divergence_set_outraster = QPushButton("...")
        self.divergence_set_outraster.clicked.connect(self.select_output_rasterFile)
        vectorfieldoperators_layout.addWidget(self.divergence_set_outraster, 2, 2, 1, 1)
        
        self.curlmodule_calc_choice_checkBox = QCheckBox("Curl module")
        vectorfieldoperators_layout.addWidget(self.curlmodule_calc_choice_checkBox, 3, 0, 1, 1)
        
        self.curlmodule_outraster_lineEdit = QLineEdit()
        self.curlmodule_outraster_lineEdit.setPlaceholderText("curmod.asc")
        vectorfieldoperators_layout.addWidget(self.curlmodule_outraster_lineEdit, 3, 1, 1, 1)
        
        self.curlmodule_set_outraster = QPushButton("...")
        self.curlmodule_set_outraster.clicked.connect(self.select_output_rasterFile)
        vectorfieldoperators_layout.addWidget(self.curlmodule_set_outraster, 3, 2, 1, 1)
        
        self.calculate_results_vfop_pb = QPushButton("Calculate")
        self.calculate_results_vfop_pb.clicked.connect(self.calculate_vectorfieldops)
        vectorfieldoperators_layout.addWidget(self.calculate_results_vfop_pb, 4, 0, 1, 2)
        
        self.output_vfops_load_choice_checkBox = QCheckBox("Load output in project")
        vectorfieldoperators_layout.addWidget(self.output_vfops_load_choice_checkBox, 4, 2, 1, 1)
        
        vectorfieldoperators_tab.setLayout(vectorfieldoperators_layout)
        
        return vectorfieldoperators_tab
 
    def setup_magnitude_gradients_widget(self):
        
        magnitudegradients_tab = QWidget()
        magnitudegradients_layout = QGridLayout() 
        
        self.gradient_x_calc_choice_checkBox = QCheckBox("X-axis")
        magnitudegradients_layout.addWidget(self.gradient_x_calc_choice_checkBox, 0, 0, 1, 1)
        
        self.gradient_x_outraster_lineEdit = QLineEdit()
        self.gradient_x_outraster_lineEdit.setPlaceholderText("grad_x.asc")
        magnitudegradients_layout.addWidget(self.gradient_x_outraster_lineEdit, 0, 1, 1, 1)
        
        self.gradient_x_set_outraster = QPushButton("...")
        self.gradient_x_set_outraster.clicked.connect(self.select_output_rasterFile)
        magnitudegradients_layout.addWidget(self.gradient_x_set_outraster, 0, 2, 1, 1)
        
        self.gradient_y_calc_choice_checkBox = QCheckBox("Y-axis")
        magnitudegradients_layout.addWidget(self.gradient_y_calc_choice_checkBox, 1, 0, 1, 1)
        
        self.gradient_y_outraster_lineEdit = QLineEdit()
        self.gradient_y_outraster_lineEdit.setPlaceholderText("grad_y.asc")
        magnitudegradients_layout.addWidget(self.gradient_y_outraster_lineEdit, 1, 1, 1, 1)
        
        self.gradient_y_set_outraster = QPushButton("...")
        self.gradient_y_set_outraster.clicked.connect(self.select_output_rasterFile)
        magnitudegradients_layout.addWidget(self.gradient_y_set_outraster, 1, 2, 1, 1)
         
        self.gradient_flowlines_calc_choice_checkBox = QCheckBox("Flowlines")
        magnitudegradients_layout.addWidget(self.gradient_flowlines_calc_choice_checkBox, 2, 0, 1, 1)
        
        self.gradient_flowlines_outraster_lineEdit = QLineEdit()
        self.gradient_flowlines_outraster_lineEdit.setPlaceholderText("grad_flw.asc")
        magnitudegradients_layout.addWidget(self.gradient_flowlines_outraster_lineEdit, 2, 1, 1, 1)
        
        self.gradient_flowlines_set_outraster = QPushButton("...")
        self.gradient_flowlines_set_outraster.clicked.connect(self.select_output_rasterFile)
        magnitudegradients_layout.addWidget(self.gradient_flowlines_set_outraster, 2, 2, 1, 1)
        
        self.calculate_results_grad_pb = QPushButton("Calculate")
        self.calculate_results_grad_pb.clicked.connect(self.calculate_vectorfieldgrads)
        magnitudegradients_layout.addWidget(self.calculate_results_grad_pb, 4, 0, 1, 2)
        
        self.output_grads_load_choice_checkBox = QCheckBox("Load output in project")
        magnitudegradients_layout.addWidget(self.output_grads_load_choice_checkBox, 4, 2, 1, 1)
     
        magnitudegradients_tab.setLayout(magnitudegradients_layout)
        
        return magnitudegradients_tab              
 
    def setup_pathline_calculation_widget(self):

        pathline_tab = QWidget()
        pathline_layout = QGridLayout()

        pathline_layout.addWidget(QLabel("Input shapefile"), 0, 0, 1, 2)
        
        self.pathlines_inshapefile_input = QLineEdit()
        self.pathlines_inshapefile_input.setPlaceholderText("in_points.shp")
        pathline_layout.addWidget(self.pathlines_inshapefile_input, 0, 2, 1, 4)
        
        self.calc_pathlines_set_inshape = QPushButton("...")
        self.calc_pathlines_set_inshape.clicked.connect(self.select_input_vectorFile)
        pathline_layout.addWidget(self.calc_pathlines_set_inshape, 0, 6, 1, 1) 

        pathline_layout.addWidget(QLabel("Time step"), 1, 0, 1, 1)
                
        self.pathlines_timestep_input = QLineEdit()
        self.pathlines_timestep_input.setPlaceholderText("0.1")
        pathline_layout.addWidget(self.pathlines_timestep_input, 1, 1, 1, 2)

        pathline_layout.addWidget(QLabel("Total time"), 1, 3, 1, 1)
        
        self.pathlines_totaltime_input = QLineEdit()
        self.pathlines_totaltime_input.setPlaceholderText("50")
        pathline_layout.addWidget(self.pathlines_totaltime_input, 1, 4, 1, 1)
        
        pathline_layout.addWidget(QLabel("Max error"), 1, 5, 1, 1)
        
        self.pathlines_stepmaxerror_input = QLineEdit()
        self.pathlines_stepmaxerror_input.setPlaceholderText("1e-6")
        pathline_layout.addWidget(self.pathlines_stepmaxerror_input, 1, 6, 1, 1)        
                
        pathline_layout.addWidget(QLabel("Output shapefile"), 2, 0, 1, 2)
        
        self.pathlines_outshapefile_input = QLineEdit()
        self.pathlines_outshapefile_input.setPlaceholderText("pathline.shp")
        pathline_layout.addWidget(self.pathlines_outshapefile_input, 2, 2, 1, 4)
        
        self.calc_pathlines_set_outshape = QPushButton("...")
        self.calc_pathlines_set_outshape.clicked.connect(self.select_output_vectorFile)
        pathline_layout.addWidget(self.calc_pathlines_set_outshape, 2, 6, 1, 1)

        self.calculate_results_pathl_pb = QPushButton("Calculate")
        self.calculate_results_pathl_pb.clicked.connect(self.calculate_vectorfieldpathlines)
        pathline_layout.addWidget(self.calculate_results_pathl_pb, 3, 0, 1, 6)
        
        self.output_pathl_load_choice_checkBox = QCheckBox("Load output in project")
        pathline_layout.addWidget(self.output_pathl_load_choice_checkBox, 3, 6, 1, 1)
        
 
        pathline_tab.setLayout(pathline_layout)
        
        return pathline_tab
             
    def setup_calculations_widget(self):
        
        calculations_widget = QTabWidget() 

        vectorfieldoperators_tab = self.setup_vector_operators_widget()
        calculations_widget.addTab(vectorfieldoperators_tab, "Vector field operators")
         
        magnitudegradients_tab = self.setup_magnitude_gradients_widget()
        calculations_widget.addTab(magnitudegradients_tab, "Magnitude gradients")
        
        pathline_tab = self.setup_pathline_calculation_widget()
        calculations_widget.addTab(pathline_tab, "Pathlines")         

        return calculations_widget 
         
    def setup_manage_widget(self):
        
        manage_window = QWidget()
        manage_layout = QHBoxLayout()
        manage_window.setLayout(manage_layout)
        
        self.help_pb = QPushButton("Help")
        self.help_pb.clicked.connect(self.open_help)
        manage_layout.addWidget(self.help_pb)

        self.about_pb = QPushButton("About")
        self.about_pb.clicked.connect(self.about)
        manage_layout.addWidget(self.about_pb)       
                
        return manage_window
         
    def select_input_vectorFile(self):
        
        sender = self.sender()
            
        fileName, _ = QFileDialog.getOpenFileName(self, "Open shapefile", lastUsedDir(), "shp (*.shp *.SHP)")
        if not fileName:
            return
        setLastUsedDir(fileName)
    
        if sender == self.calc_pathlines_set_inshape: 
            self.pathlines_inshapefile_input.setText(fileName)
    
    def select_output_rasterFile(self):
        # modified after RASTERCALC module by Barry Rowlingson    
            
        sender = self.sender()
            
        fileName, _ = QFileDialog.getSaveFileName(self, "Save ESRI grid ascii file", lastUsedDir(), "asc (*.asc *.ASC)")
        if not fileName:
            return
          
        setLastUsedDir(fileName)
    
        if sender == self.magnitude_set_outraster: 
            self.magnitude_outraster_lineEdit.setText(fileName)
        elif sender == self.orientations_set_outraster: 
            self.orientations_outraster_lineEdit.setText(fileName)            
        elif sender == self.divergence_set_outraster: 
            self.divergence_outraster_lineEdit.setText(fileName)            
        elif sender == self.curlmodule_set_outraster: 
            self.curlmodule_outraster_lineEdit.setText(fileName) 
        elif sender == self.gradient_x_set_outraster: 
            self.gradient_x_outraster_lineEdit.setText(fileName)
        elif sender == self.gradient_y_set_outraster: 
            self.gradient_y_outraster_lineEdit.setText(fileName)
        elif sender == self.gradient_flowlines_set_outraster: 
            self.gradient_flowlines_outraster_lineEdit.setText(fileName)             
             
    def select_output_vectorFile(self):
        
        sender = self.sender()
            
        fileName, _ = QFileDialog.getSaveFileName(self,
                                               "Save shapefile", 
                                               lastUsedDir(), 
                                               "shp (*.shp *.SHP)")
        if not fileName:
            return
        setLastUsedDir(fileName)
    
        if sender == self.calc_pathlines_set_outshape: 
            self.pathlines_outshapefile_input.setText(fileName)                

    def get_raster_parameters(self):
        
        return self.inraster_x_comboBox.currentText(), self.inraster_y_comboBox.currentText()

    def get_vector_operator_parameters(self):
    
        # vector operator parameters            
        return (self.magnitude_calc_choice_checkBox.isChecked(), self.magnitude_outraster_lineEdit.text(),
                self.orientations_calc_choice_checkBox.isChecked(), self.orientations_outraster_lineEdit.text(),
                self.divergence_calc_choice_checkBox.isChecked(), self.divergence_outraster_lineEdit.text(),
                self.curlmodule_calc_choice_checkBox.isChecked(), self.curlmodule_outraster_lineEdit.text())

    def get_gradient_parameters(self):
            
        # vector field parameters            
        return (self.gradient_x_calc_choice_checkBox.isChecked(), self.gradient_x_outraster_lineEdit.text(),
                self.gradient_y_calc_choice_checkBox.isChecked(), self.gradient_y_outraster_lineEdit.text(),
                self.gradient_flowlines_calc_choice_checkBox.isChecked(), self.gradient_flowlines_outraster_lineEdit.text())

    def check_inrasters_def(self, inraster_x, inraster_y):
        
        if inraster_x is None or inraster_x == '':  
            return False, "No layer defined for x components" 
         
        if inraster_y is None or inraster_y == '':  
            return False, "No layer defined for y components"
        
        if inraster_x == inraster_y:                
            msgBox = QMessageBox()
            msgBox.setText("The two input rasters are the same.")
            msgBox.setInformativeText("Is this correct?")
            msgBox.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)
            msgBox.setDefaultButton(QMessageBox.Cancel)
            ret = msgBox.exec_()
            
            if ret == QMessageBox.Yes:
                return True, "Equal rasters"
            else:
                return False, "Redefine x- and y- input rasters"         

        return True, "Ok"  

    def get_raster_data(self, inraster_comp):        

        comp_params, comp_array = read_raster_layer(
            inraster_comp,
            iter(list(self.layermap.items())))
        params_correct, msg = comp_params.check_params()        
        if params_correct:
            return True, (comp_params, comp_array)
        else:
            return False, msg

    def calculate_vectorfieldops(self):

        # input rasters            
        inraster_x, inraster_y = self.get_raster_parameters()

        # vector field parameters            
        magnitude_calc_choice, magnitude_outraster_path, \
        orientations_calc_choice, orientations_outraster_path, \
        divergence_calc_choice, divergence_outraster_path, \
        curlmodule_calc_choice, curlmodule_outraster_path  = self.get_vector_operator_parameters() 

        # load results                       
        load_output = self.output_vfops_load_choice_checkBox.isChecked()  

        ### PRE-PROCESSINGS
        
        # check if any calculation            
        if not (magnitude_calc_choice or orientations_calc_choice or divergence_calc_choice or curlmodule_calc_choice):
            QMessageBox.critical(self, "Vector field processing", "No parameter to calculate")   
            return

        # check input values for vector field rasters
        success, msg = self.check_inrasters_def(inraster_x, inraster_y)
        if not success:
            QMessageBox.critical(self, "Input component rasters", msg)
            return
        
        # get x- and y-axis component data
        success, result = self.get_raster_data(inraster_x)
        if not success:
            QMessageBox.critical(self, "Input x-component raster", msg)
            return   
        else:
            comp_x_params, comp_x_array = result
                     
        success, result = self.get_raster_data(inraster_y)
        if not success:
            QMessageBox.critical(self, "Input y-component raster", msg)
            return   
        else:
            comp_y_params, comp_y_array = result        
            
        # check geometric and geographic equivalence of the two components rasters
        if not comp_x_params.geo_equiv(comp_y_params):
            QMessageBox.critical(self, "Raster component grids", "The two rasters have different geographic extent and/or cell sizes")   
            return 
        
        # pre-processes vector field parameters            
        params_choices = (magnitude_calc_choice, orientations_calc_choice, divergence_calc_choice, curlmodule_calc_choice)     
        params_names = ('magnitude','orientation','divergence','curl module')                      
        params_savefiles = (magnitude_outraster_path, orientations_outraster_path, divergence_outraster_path, curlmodule_outraster_path)
        params_functions = ('magnitude', 'orientations', 'divergence', 'curl_module')
    
        # verify input choices for vector field parameters            
        for vfp_name, vfp_choice, vfp_savefile in zip(params_names, params_choices, params_savefiles):
            if vfp_choice and (vfp_savefile is None or vfp_savefile == ''):
                QMessageBox.critical(self, vfp_name + " values", "No output layer defined")   
                return 
        
        ### PROCESSINGS
        
        # create velocity vector field            
        vector_array = np.zeros((comp_x_params.get_rows(), comp_x_params.get_cols(), 2))
        vector_array[:,:,0], vector_array[:,:,1] = comp_x_array, comp_y_array

        # calculates vector field parameters 
        for vfp_name, vfp_choice, vfp_savefile, vfp_function in zip(params_names, params_choices, params_savefiles, params_functions):  
            if vfp_choice:
                exec("curr_fld = vector_field.%s()" % vfp_function)
                exec("curr_fld.write_esrigrid('%s')" % vfp_savefile)
        
        # add required layer to the map canvas - modified after RasterCalc module                
        if load_output:
            # vector field parameters
            for vfc_choice, vfc_savefile in zip(params_choices, params_savefiles):  
                if vfc_choice:
                    newLayer = QgsRasterLayer(vfc_savefile, QFileInfo(vfc_savefile).baseName())
                    QgsProject.instance().addMapLayer(newLayer)

        # all done
        QMessageBox.information(self, "Vector operators output", "Processings completed.")

    def calculate_vectorfieldgrads(self):
        
        # input rasters            
        inraster_x, inraster_y = self.get_raster_parameters()

        # gradients parameters                                           
        gradient_x_calc_choice, gradient_x_outraster_path, \
        gradient_y_calc_choice, gradient_y_outraster_path, \
        gradient_flowlines_calc_choice, gradient_flowlines_outraster_path  = self.get_gradient_parameters() 

        # load results                       
        load_output = self.output_grads_load_choice_checkBox.isChecked()

        ### PRE-PROCESSINGS
        
        # check if any calculation                   
        if not (gradient_x_calc_choice or gradient_y_calc_choice or gradient_flowlines_calc_choice):
            QMessageBox.critical(self, "Vector field gradient processing", "No parameter to calculate")   
            return

        # check input values for vector field rasters
        success, msg = self.check_inrasters_def(inraster_x, inraster_y)
        if not success:
            QMessageBox.critical(self, "Input component rasters", msg)
            return
        
        # get x- and y-axis component data
        success, result = self.get_raster_data(inraster_x)
        if not success:
            QMessageBox.critical(self, "Input x-component raster", msg)
            return   
        else:
            comp_x_params, comp_x_array = result
                     
        success, result = self.get_raster_data(inraster_y)
        if not success:
            QMessageBox.critical(self, "Input y-component raster", msg)
            return   
        else:
            comp_y_params, comp_y_array = result        
            
        # check geometric and geographic equivalence of the two components rasters
        if not comp_x_params.geo_equiv(comp_y_params):
            QMessageBox.critical(self, "Raster component grids", "The two rasters have different geographic extent and/or cell sizes")   
            return 
        
        # pre-processes gradients parameters            
        vfgrads_choices = (gradient_x_calc_choice, gradient_y_calc_choice, gradient_flowlines_calc_choice)     
        vfgrads_names = ('gradient along x axis','gradient along y axis','gradient along flow lines')                      
        vfgrads_savefiles = (gradient_x_outraster_path, gradient_y_outraster_path, gradient_flowlines_outraster_path)
        vfgrads_functions = ('grad_xaxis', 'grad_yaxis', 'grad_flowlines')                     
                 
        # verify input choices for gradients parameters            
        for vfg_name, vfg_choice, vfg_savefile in zip(vfgrads_names, vfgrads_choices, vfgrads_savefiles):
            if vfg_choice and (vfg_savefile is None or vfg_savefile == ''):
                QMessageBox.critical(self, vfg_name + " values", "No output layer defined")   
                return 
        
        ### PROCESSINGS
        
        # create velocity vector field            
        vector_array = np.zeros((comp_x_params.get_rows(), comp_x_params.get_cols(), 2))
        vector_array[:, :, 0], vector_array[:, :, 1] = comp_x_array, comp_y_array

        # calculates gradients parameters 
        for vfg_name, vfg_choice, vfg_savefile, vfg_function in zip(vfgrads_names, vfgrads_choices, vfgrads_savefiles, vfgrads_functions):  
            if vfg_choice:
                exec("curr_fld = vector_field.%s()" % vfg_function)
                exec("curr_fld.write_esrigrid('%s')" % vfg_savefile)
    
        # add required layer to the map canvas - modified after RasterCalc module                
        if load_output:
            # gradient parameters                                              
            for vfg_choice, vfg_savefile in zip(vfgrads_choices, vfgrads_savefiles):  
                if vfg_choice:
                    newLayer = QgsRasterLayer(vfg_savefile, QFileInfo(vfg_savefile).baseName())
                    QgsMapLayerRegistry.instance().addMapLayer(newLayer)    
   
        # all done
        QMessageBox.information(self, "Gradients output", "Processings completed.")

    def calculate_vectorfieldpathlines(self):

        # input rasters            
        inraster_x, inraster_y = self.get_raster_parameters()

        # pathline calculation parameters                    
        pathlines_input_shapefile = self.pathlines_inshapefile_input.text() 
        time_step = self.pathlines_timestep_input.text() 
        total_time = self.pathlines_totaltime_input.text() 
        error_max_tolerance = self.pathlines_stepmaxerror_input.text()
        pathlines_output_shapefile = self.pathlines_outshapefile_input.text()
        
        # load results                       
        load_output = self.output_pathl_load_choice_checkBox.isChecked()

        ### PRE-PROCESSINGS
        
        # check input values for vector field rasters
        success, msg = self.check_inrasters_def(inraster_x, inraster_y)
        if not success:
            QMessageBox.critical(self, "Input component rasters", msg)
            return
        
        # get x- and y-axis component data
        success, result = self.get_raster_data(inraster_x)
        if not success:
            QMessageBox.critical(self, "Input x-component raster", msg)
            return   
        else:
            comp_x_params, comp_x_array = result
                     
        success, result = self.get_raster_data(inraster_y)
        if not success:
            QMessageBox.critical(self, "Input y-component raster", msg)
            return   
        else:
            comp_y_params, comp_y_array = result        
            
        # check geometric and geographic equivalence of the two components rasters
        if not comp_x_params.geo_equiv(comp_y_params):
            QMessageBox.critical(self, "Raster component grids", "The two rasters have different geographic extent and/or cell sizes")   
            return 

        # verify input parameters                
        if pathlines_input_shapefile is None or pathlines_input_shapefile == '':
            QMessageBox.critical(self, "Pathline calculation", "No input point layer defined")   
            return 
        try:
            time_step = float(time_step)
            total_time = float(total_time)
            assert time_step * total_time > 0.0
        except:
            QMessageBox.critical(self, "Pathline calculation", "Time step and total time should be both positive or both negative")   
            return 
        
        try:
            error_max_tolerance = float(error_max_tolerance)
            assert error_max_tolerance > 0.0
        except:
            QMessageBox.critical(self, "Pathline calculation", "Max error should be a positive number")   
            return                        
                       
        if pathlines_output_shapefile is None or pathlines_output_shapefile == '':
            QMessageBox.critical(self, "Pathline calculation", "No output point layer defined")   
            return 
            
        ### PROCESSINGS
        
        # create velocity vector field            
        vector_array = np.zeros((comp_x_params.get_rows(), comp_x_params.get_cols(), 2))
        vector_array[:,:,0], vector_array[:,:,1] = comp_x_array, comp_y_array
        
        vector_field = Grid(comp_x_params, vector_array)            
        
        pathlines_input_shapefile = str(pathlines_input_shapefile)        
        pathlines_output_shapefile = str(pathlines_output_shapefile)

        # get input and output shapefiles        
        driver = ogr.GetDriverByName('ESRI Shapefile') 
      
        in_shapefile = driver.Open(pathlines_input_shapefile, 0)                        
        if in_shapefile is None:
            QMessageBox.critical(
                self,
                "Error in pathline calculation",
                "Unable to get input shapefile")
            return
   
        if os.path.exists(pathlines_output_shapefile):
            driver.DeleteDataSource(pathlines_output_shapefile)
    
        out_shape = driver.CreateDataSource(pathlines_output_shapefile)
        if out_shape is None:
            QMessageBox.critical(
                self,
                "Pathline calculation",
                "Unable to create output shapefile: %s" % pathlines_output_shapefile)
            return
        
        out_layer = out_shape.CreateLayer('out_pathlines', geom_type=ogr.wkbPoint)
        if out_layer is None:
            QMessageBox.critical(
                self,
                "Pathline calculation",
                "Unable to create output shapefile: %s" % pathlines_output_shapefile)
            return        
    
        # add fields to the output shapefile    
        path_id_fieldDef = ogr.FieldDefn('path_id', ogr.OFTInteger)
        out_layer.CreateField(path_id_fieldDef)
    
        point_id_fieldDef = ogr.FieldDefn('point_id', ogr.OFTInteger)
        out_layer.CreateField(point_id_fieldDef)
        
        x_fieldDef = ogr.FieldDefn('x', ogr.OFTReal)
        out_layer.CreateField(x_fieldDef)
    
        y_fieldDef = ogr.FieldDefn('y', ogr.OFTReal)
        out_layer.CreateField(y_fieldDef)
    
        delta_s_fieldDef = ogr.FieldDefn('ds', ogr.OFTReal)
        out_layer.CreateField(delta_s_fieldDef)
        
        s_fieldDef = ogr.FieldDefn('s', ogr.OFTReal)
        out_layer.CreateField(s_fieldDef)
    
        vx_fieldDef = ogr.FieldDefn('vx', ogr.OFTReal)
        out_layer.CreateField(vx_fieldDef)
    
        vy_fieldDef = ogr.FieldDefn('vy', ogr.OFTReal)
        out_layer.CreateField(vy_fieldDef)
    
        vmagn_fieldDef = ogr.FieldDefn('vmagn', ogr.OFTReal)
        out_layer.CreateField(vmagn_fieldDef)
    
        dtime_fieldDef = ogr.FieldDefn('d_time', ogr.OFTReal)
        out_layer.CreateField(dtime_fieldDef)
    
        t_fieldDef = ogr.FieldDefn('t', ogr.OFTReal)
        out_layer.CreateField(t_fieldDef)
    
        error_fieldDef = ogr.FieldDefn('error', ogr.OFTReal)
        out_layer.CreateField(error_fieldDef)    
    
        # get the layer definition of the output shapefile
        outshape_featdef = out_layer.GetLayerDefn()
    
        # get info from input layer
        ptLayer = in_shapefile.GetLayer(0) 
    
        # inizialization of pathline id
        pathline_id = 0
        
        # get next feature
        pt_feature = ptLayer.GetNextFeature()
    
        while pt_feature:
            
            # update of current pathline id
            pathline_id += 1
            
            # inizializations for new point pathline
            pathline_cumulated_time = 0.0
            pathline_cumulated_length = 0.0
            pathline_point_id = 0
            interp_Pt_error_estimate = 0.0
            delta_time = time_step
            delta_s = 0.0
            
            # new point with coords and path_id from input layer
            curr_Pt = Point(pt_feature.GetGeometryRef().GetX(), pt_feature.GetGeometryRef().GetY())
            
            # case where the start point is outside the velocity field
            if not vector_field.include_point_location(curr_Pt):
                # get next feature
                pt_feature = ptLayer.GetNextFeature()
                continue                        
                
            # pathline cycle
            while abs(pathline_cumulated_time) < abs(total_time):
                
                if vector_field.include_point_location(curr_Pt):            
                    
                    # projects current point in grid coords
                    curr_Pt_gridcoord = vector_field.geog2gridcoord(curr_Pt) 
    
                    # interpolate velocity components for current point 
                    curr_Pt_vx = vector_field.interpolate_level_bilinear(0, curr_Pt_gridcoord)
                    curr_Pt_vy = vector_field.interpolate_level_bilinear(1, curr_Pt_gridcoord)
                               
                    # velocity magnitude
                    curr_v_magnitude = sqrt(curr_Pt_vx**2 + curr_Pt_vy**2)            
    
                    # update total length
                    if abs(pathline_cumulated_time) > 0:
                        delta_s = curr_Pt.distance(self.prev_Pt)
                        pathline_cumulated_length += delta_s
    
                    # pre-processing for new feature in output layer
                    curr_Pt_geom = ogr.Geometry(ogr.wkbPoint)
                    curr_Pt_geom.AddPoint(float(curr_Pt.x), float(curr_Pt.y))
                        
                    # create a new feature
                    curr_Pt_shape = ogr.Feature(outshape_featdef)
                    curr_Pt_shape.SetGeometry(curr_Pt_geom)
                    curr_Pt_shape.SetField('path_id', pathline_id) 
                    curr_Pt_shape.SetField('point_id', pathline_point_id)                            
                    curr_Pt_shape.SetField('x', curr_Pt.x)
                    curr_Pt_shape.SetField('y', curr_Pt.y) 
                    curr_Pt_shape.SetField('ds', delta_s) 
                    curr_Pt_shape.SetField('s', pathline_cumulated_length)                            
                    curr_Pt_shape.SetField('t', pathline_cumulated_time)
                    curr_Pt_shape.SetField('vx', curr_Pt_vx)
                    curr_Pt_shape.SetField('vy', curr_Pt_vy)
                    curr_Pt_shape.SetField('vmagn', curr_v_magnitude)        
                    curr_Pt_shape.SetField('d_time', delta_time)            
                    curr_Pt_shape.SetField('error', interp_Pt_error_estimate)         
    
                    # add the feature to the output layer
                    out_layer.CreateFeature(curr_Pt_shape)            
                    
                    # destroy no longer used objects
                    curr_Pt_geom.Destroy()
                    curr_Pt_shape.Destroy()
                    
                    # store coordinates for pathline length calculation
                    self.prev_Pt = Point(curr_Pt.x, curr_Pt.y)
                    
                    # when possible, doubles the delta time value 
                    if interp_Pt_error_estimate < error_max_tolerance/100.0:
                        delta_time *= 2.0                
                    
                    # interpolate new location
                    interp_Pt, interp_Pt_error_estimate = vector_field.interpolate_RKF(delta_time, curr_Pt)                            
                    if interp_Pt is None or interp_Pt_error_estimate is None:
                        break                             
                    while interp_Pt_error_estimate > error_max_tolerance:
                        delta_time /= 2.0
                        interp_Pt, interp_Pt_error_estimate = vector_field.interpolate_RKF(delta_time, curr_Pt)                
                       
                    # update current Point coords with those of the interpolated point
                    curr_Pt.x = interp_Pt.x
                    curr_Pt.y = interp_Pt.y
                      
                    # update total time and point counter
                    pathline_cumulated_time += delta_time
                    pathline_point_id += 1
                    
            # get next feature
            pt_feature = ptLayer.GetNextFeature()
        
        # destroy output geometry
        out_shape.Destroy()

        # add required layer to the map canvas - modified after RasterCalc module                
        if load_output:
            pathlines_outshape = QgsVectorLayer(pathlines_output_shapefile, 
                                                QFileInfo(pathlines_output_shapefile).baseName(), 
                                                "ogr")                    
            QgsMapLayerRegistry.instance().addMapLayer(pathlines_outshape)

        # all done
        QMessageBox.information(self, "Pathline output", "Processings completed.")

    def open_help(self):
        # modified after CADTOOLS module   
             
        help_path = os.path.join(os.path.dirname(__file__), 'help', 'help.html')
        webbrowser.open(help_path) 

    def about(self):
        """
        Visualize an About window.
        """
        
        QMessageBox.about(self, "About VectorFieldCalc", 
        """
            <p>VectorFieldCalc version 1.4<br />License: GPL v. 3</p>
            <p>Mauro Alberti</p> 
            <p>This application calculates vector field parameters (e.g., divergence, curl module, gradients)
            and pathlines.            
            </p>
             <p>Please report any bug to <a href="mailto:alberti.m65@gmail.com">alberti.m65@gmail.com</a></p>
        """)              

