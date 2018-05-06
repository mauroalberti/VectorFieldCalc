# -*- coding: utf-8 -*-


import os, webbrowser
from osgeo import ogr
    
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

from qgis.core import *

from pygsf.spatial.rasters.gdal_utils import *
from pygsf.spatial.vectorial.vectorial import Point

from libs_utils.qt.filesystem import *


class VfcDialog(QDialog):
    
    def __init__(self, plugin_name):
        
        super().__init__(self)
        self.plugin_nm = plugin_name
        
        self.qlytDialogLayout = QVBoxLayout()
        self.setLayout(self.qlytDialogLayout)
        
        self.qgrpbxInputRasterWidget = self.setupInputRasterWidget()
        self.qlytDialogLayout.addWidget(self.qgrpbxInputRasterWidget)
        
        self.qtabwCalculationsWidget = self.setupCalculationsWidget()
        self.qlytDialogLayout.addWidget(self.qtabwCalculationsWidget)
        self.qtabwCalculationsWidget.setCurrentIndex(0)

        self.qwdgtManageWdgt = self.setupManageWidget()
        self.qlytDialogLayout.addWidget(self.qwdgtManageWdgt)

        self.adjustSize()

        self.setWindowTitle(self.plugin_nm)

    def setupInputRasterWidget(self):
        
        qgrpbxInRaster = QGroupBox("Input rasters")
        qfrmlytInRaster = QFormLayout()

        self.qcmbbxInRasterX = QComboBox()
        qfrmlytInRaster.addRow(QLabel("X-axis components"), self.qcmbbxInRasterX)
        self.qcmbbxInRasterY = QComboBox()
        qfrmlytInRaster.addRow(QLabel("Y-axis components"), self.qcmbbxInRasterY)
                
        # append loaded rasters layers to combo boxes
        self.qmapLayersMap = QgsMapLayerRegistry.instance().mapLayers()
        for name, layer in self.qmapLayersMap.iteritems():
            if layer.type() == QgsMapLayer.RasterLayer: 
                self.qcmbbxInRasterX.addItem(layer.name())
                self.qcmbbxInRasterY.addItem(layer.name())
 
        qgrpbxInRaster.setLayout(qfrmlytInRaster)
               
        return qgrpbxInRaster

    def setupVectorOperatorsWidget(self):

        qwdgVectorFieldOperators = QWidget()
        qlytVectorFieldOperators = QGridLayout()
        
        self.qchkbxMagnitudeCalcChoice = QCheckBox("Magnitude")
        qlytVectorFieldOperators.addWidget(self.qchkbxMagnitudeCalcChoice, 0, 0, 1, 1)
        
        self.qlndtMagnitudeOutRaster = QLineEdit()
        self.qlndtMagnitudeOutRaster.setPlaceholderText("magn.asc")
        qlytVectorFieldOperators.addWidget(self.qlndtMagnitudeOutRaster, 0, 1, 1, 1)
        
        self.qpshbtMagnitudeOutRaster = QPushButton("...")
        self.qpshbtMagnitudeOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytVectorFieldOperators.addWidget(self.qpshbtMagnitudeOutRaster, 0, 2, 1, 1)

        self.qchkbxOrientationsCalcChoice = QCheckBox("Orientations")
        qlytVectorFieldOperators.addWidget(self.qchkbxOrientationsCalcChoice, 1, 0, 1, 1)
        
        self.qlnedtOrientationsOutRaster = QLineEdit()
        self.qlnedtOrientationsOutRaster.setPlaceholderText("orient.asc")
        qlytVectorFieldOperators.addWidget(self.qlnedtOrientationsOutRaster, 1, 1, 1, 1)
        
        self.qpshbtOrientationsOutRaster = QPushButton("...")
        self.qpshbtOrientationsOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytVectorFieldOperators.addWidget(self.qpshbtOrientationsOutRaster, 1, 2, 1, 1)
        
        self.qchkbxDivergenceCalcChoice = QCheckBox("Divergence")
        qlytVectorFieldOperators.addWidget(self.qchkbxDivergenceCalcChoice, 2, 0, 1, 1)
        
        self.qlnedtDivergenceOutRaster = QLineEdit()
        self.qlnedtDivergenceOutRaster.setPlaceholderText("div.asc")
        qlytVectorFieldOperators.addWidget(self.qlnedtDivergenceOutRaster, 2, 1, 1, 1)
        
        self.qpshbtDivergenceOutRaster = QPushButton("...")
        self.qpshbtDivergenceOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytVectorFieldOperators.addWidget(self.qpshbtDivergenceOutRaster, 2, 2, 1, 1)
        
        self.qchkbxCurlModuleCalcChoice = QCheckBox("Curl module")
        qlytVectorFieldOperators.addWidget(self.qchkbxCurlModuleCalcChoice, 3, 0, 1, 1)
        
        self.qlnedtCurlModuleOutRaster = QLineEdit()
        self.qlnedtCurlModuleOutRaster.setPlaceholderText("curmod.asc")
        qlytVectorFieldOperators.addWidget(self.qlnedtCurlModuleOutRaster, 3, 1, 1, 1)
        
        self.qpshbtCurlModuleOutRaster = QPushButton("...")
        self.qpshbtCurlModuleOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytVectorFieldOperators.addWidget(self.qpshbtCurlModuleOutRaster, 3, 2, 1, 1)
        
        self.qpshbtCalculateResultsVfop = QPushButton("Calculate")
        self.qpshbtCalculateResultsVfop.clicked.connect(self.calculateVectorFieldOps)
        qlytVectorFieldOperators.addWidget(self.qpshbtCalculateResultsVfop, 4, 0, 1, 2)
        
        self.qchkbxOutputVfopsLoadChoice = QCheckBox("Load output in project")
        qlytVectorFieldOperators.addWidget(self.qchkbxOutputVfopsLoadChoice, 4, 2, 1, 1)
        
        qwdgVectorFieldOperators.setLayout(qlytVectorFieldOperators)
        
        return qwdgVectorFieldOperators
 
    def setupMagnitudeGradientsWidget(self):
        
        qwdgtMagnitudeGradients = QWidget()
        qlytMagnitudeGradients = QGridLayout()
        
        self.qchkbxGradientXCalcChoice = QCheckBox("X-axis")
        qlytMagnitudeGradients.addWidget(self.qchkbxGradientXCalcChoice, 0, 0, 1, 1)
        
        self.qlnedtGradientXOutRaster = QLineEdit()
        self.qlnedtGradientXOutRaster.setPlaceholderText("grad_x.asc")
        qlytMagnitudeGradients.addWidget(self.qlnedtGradientXOutRaster, 0, 1, 1, 1)
        
        self.qpshbtGradientXOutRaster = QPushButton("...")
        self.qpshbtGradientXOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytMagnitudeGradients.addWidget(self.qpshbtGradientXOutRaster, 0, 2, 1, 1)
        
        self.qchkbxGradientYCalcChoice = QCheckBox("Y-axis")
        qlytMagnitudeGradients.addWidget(self.qchkbxGradientYCalcChoice, 1, 0, 1, 1)
        
        self.qlnedtGradientYOutRaster = QLineEdit()
        self.qlnedtGradientYOutRaster.setPlaceholderText("grad_y.asc")
        qlytMagnitudeGradients.addWidget(self.qlnedtGradientYOutRaster, 1, 1, 1, 1)
        
        self.qpshbtGradientYOutRaster = QPushButton("...")
        self.qpshbtGradientYOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytMagnitudeGradients.addWidget(self.qpshbtGradientYOutRaster, 1, 2, 1, 1)
         
        self.qchkbxGradientFlowlinesCalcChoice = QCheckBox("Flowlines")
        qlytMagnitudeGradients.addWidget(self.qchkbxGradientFlowlinesCalcChoice, 2, 0, 1, 1)
        
        self.qlnedtGradientFlowlinesOutRaster = QLineEdit()
        self.qlnedtGradientFlowlinesOutRaster.setPlaceholderText("grad_flw.asc")
        qlytMagnitudeGradients.addWidget(self.qlnedtGradientFlowlinesOutRaster, 2, 1, 1, 1)
        
        self.qpshbtGradientFlowlinesOutRaster = QPushButton("...")
        self.qpshbtGradientFlowlinesOutRaster.clicked.connect(self.selectOutputRasterFile)
        qlytMagnitudeGradients.addWidget(self.qpshbtGradientFlowlinesOutRaster, 2, 2, 1, 1)
        
        self.qpshbtCalculateResultsGrad = QPushButton("Calculate")
        self.qpshbtCalculateResultsGrad.clicked.connect(self.calculate_vectorfieldgrads)
        qlytMagnitudeGradients.addWidget(self.qpshbtCalculateResultsGrad, 4, 0, 1, 2)
        
        self.qchkbxOutputGradsLoadChoice = QCheckBox("Load output in project")
        qlytMagnitudeGradients.addWidget(self.qchkbxOutputGradsLoadChoice, 4, 2, 1, 1)
     
        qwdgtMagnitudeGradients.setLayout(qlytMagnitudeGradients)
        
        return qwdgtMagnitudeGradients

    def setupPathlineCalculationWidget(self):

        qwdgtPathline = QWidget()
        qlytPathline = QGridLayout()

        qlytPathline.addWidget(QLabel("Input shapefile"), 0, 0, 1, 2)
        
        self.qlnedtPathlinesInShapefile = QLineEdit()
        self.qlnedtPathlinesInShapefile.setPlaceholderText("in_points.shp")
        qlytPathline.addWidget(self.qlnedtPathlinesInShapefile, 0, 2, 1, 4)
        
        self.qpshbtCalcPathlinesInShape = QPushButton("...")
        self.qpshbtCalcPathlinesInShape.clicked.connect(self.selectInputVectorFile)
        qlytPathline.addWidget(self.qpshbtCalcPathlinesInShape, 0, 6, 1, 1)

        qlytPathline.addWidget(QLabel("Time step"), 1, 0, 1, 1)
                
        self.qlnedtPathlinesTimeStepInput = QLineEdit()
        self.qlnedtPathlinesTimeStepInput.setPlaceholderText("0.1")
        qlytPathline.addWidget(self.qlnedtPathlinesTimeStepInput, 1, 1, 1, 2)

        qlytPathline.addWidget(QLabel("Total time"), 1, 3, 1, 1)
        
        self.qlnedtPathlinesTotalTimeInput = QLineEdit()
        self.qlnedtPathlinesTotalTimeInput.setPlaceholderText("50")
        qlytPathline.addWidget(self.qlnedtPathlinesTotalTimeInput, 1, 4, 1, 1)
        
        qlytPathline.addWidget(QLabel("Max error"), 1, 5, 1, 1)
        
        self.qlnedtPathlinesStepMaxErrorInput = QLineEdit()
        self.qlnedtPathlinesStepMaxErrorInput.setPlaceholderText("1e-6")
        qlytPathline.addWidget(self.qlnedtPathlinesStepMaxErrorInput, 1, 6, 1, 1)
                
        qlytPathline.addWidget(QLabel("Output shapefile"), 2, 0, 1, 2)
        
        self.qlnedtPathlinesOutShapefileInput = QLineEdit()
        self.qlnedtPathlinesOutShapefileInput.setPlaceholderText("pathline.shp")
        qlytPathline.addWidget(self.qlnedtPathlinesOutShapefileInput, 2, 2, 1, 4)
        
        self.qpshbtCalcPathlinesOutShape = QPushButton("...")
        self.qpshbtCalcPathlinesOutShape.clicked.connect(self.selectOutputVectorFile)
        qlytPathline.addWidget(self.qpshbtCalcPathlinesOutShape, 2, 6, 1, 1)

        self.qpshbtCalculateResultsPathline = QPushButton("Calculate")
        self.qpshbtCalculateResultsPathline.clicked.connect(self.calculateVectorFieldPathlines)
        qlytPathline.addWidget(self.qpshbtCalculateResultsPathline, 3, 0, 1, 6)
        
        self.qchkbxOutputPathlineLoadChoice = QCheckBox("Load output in project")
        qlytPathline.addWidget(self.qchkbxOutputPathlineLoadChoice, 3, 6, 1, 1)

        qwdgtPathline.setLayout(qlytPathline)
        
        return qwdgtPathline

    def setupCalculationsWidget(self):
        
        qtbwgtCalculations = QTabWidget()

        qwgtVectorFieldOperators = self.setupVectorOperatorsWidget()
        qtbwgtCalculations.addTab(qwgtVectorFieldOperators, "Vector field operators")
         
        qwgtMagnitudeGradients = self.setupMagnitudeGradientsWidget()
        qtbwgtCalculations.addTab(qwgtMagnitudeGradients, "Magnitude gradients")
        
        qwgtPathline = self.setupPathlineCalculationWidget()
        qtbwgtCalculations.addTab(qwgtPathline, "Pathlines")

        return qtbwgtCalculations

    def setupManageWidget(self):
        
        qwgtManage = QWidget()
        qlytManage = QHBoxLayout()
        qwgtManage.setLayout(qlytManage)
        
        self.qpshbtHelp = QPushButton("Help")
        self.qpshbtHelp.clicked.connect(self.open_help)
        qlytManage.addWidget(self.qpshbtHelp)

        self.qpshbtAbout = QPushButton("About")
        self.qpshbtAbout.clicked.connect(self.about)
        qlytManage.addWidget(self.qpshbtAbout)
                
        return qwgtManage

    def selectInputVectorFile(self):
        
        sender = self.sender()
            
        fileName = QFileDialog.getOpenFileName(
            self,
            "Open shapefile",
            lastUsedDir(self.plugin_nm),
            "shp (*.shp *.SHP)")

        if not fileName:
            return
        setLastUsedDir(
            self.plugin_nm,
            fileName)
    
        if sender == self.qpshbtCalcPathlinesInShape:
            self.qlnedtPathlinesInShapefile.setText(fileName)
    
    def selectOutputRasterFile(self):
        """
        Modified after RASTERCALC module by Barry Rowlingson

        :return:
        """

        sender = self.sender()
            
        fileName = QFileDialog.getSaveFileName(
            self,
            "Save ESRI grid ascii file",
            lastUsedDir(self.plugin_nm),
            "asc (*.asc *.ASC)")

        if not fileName:
            return
          
        setLastUsedDir(
            self.plugin_nm,
            fileName)
    
        if sender == self.qpshbtMagnitudeOutRaster:
            self.qlndtMagnitudeOutRaster.setText(fileName)
        elif sender == self.qpshbtOrientationsOutRaster: 
            self.qlnedtOrientationsOutRaster.setText(fileName)
        elif sender == self.qpshbtDivergenceOutRaster: 
            self.qlnedtDivergenceOutRaster.setText(fileName)            
        elif sender == self.qpshbtCurlModuleOutRaster: 
            self.qlnedtCurlModuleOutRaster.setText(fileName) 
        elif sender == self.qpshbtGradientXOutRaster:
            self.qlnedtGradientXOutRaster.setText(fileName)
        elif sender == self.qpshbtGradientYOutRaster:
            self.qlnedtGradientYOutRaster.setText(fileName)
        elif sender == self.qpshbtGradientFlowlinesOutRaster:
            self.qlnedtGradientFlowlinesOutRaster.setText(fileName)

    def selectOutputVectorFile(self):
        
        sender = self.sender()
            
        fileName = QFileDialog.getSaveFileName(
            self,
            "Save shapefile",
            lastUsedDir(self.plugin_nm),
            "shp (*.shp *.SHP)")
        if not fileName:
            return
        setLastUsedDir(
            self.plugin_nm,
            fileName)
    
        if sender == self.qpshbtCalcPathlinesOutShape:
            self.qlnedtPathlinesOutShapefileInput.setText(fileName)


    def get_raster_parameters(self):
        
        return self.qcmbbxInRasterX.currentText(), self.qcmbbxInRasterY.currentText()
   

    def get_vector_operator_parameters(self):
    
        # vector operator parameters            
        return (self.qchkbxMagnitudeCalcChoice.isChecked(), self.qlndtMagnitudeOutRaster.text(),
                self.qchkbxOrientationsCalcChoice.isChecked(), self.qlnedtOrientationsOutRaster.text(),
                self.qchkbxDivergenceCalcChoice.isChecked(), self.qlnedtDivergenceOutRaster.text(),
                self.qchkbxCurlModuleCalcChoice.isChecked(), self.qlnedtCurlModuleOutRaster.text())
    

    def get_gradient_parameters(self):
            
        # vector field parameters            
        return (self.qchkbxGradientXCalcChoice.isChecked(), self.qlnedtGradientXOutRaster.text(),
                self.qchkbxGradientYCalcChoice.isChecked(), self.qlnedtGradientYOutRaster.text(),
                self.qchkbxGradientFlowlinesCalcChoice.isChecked(), self.qlnedtGradientFlowlinesOutRaster.text())
    
    
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

        comp_params, comp_array = read_raster_layer(inraster_comp,
                                                    self.qmapLayersMap.iteritems())
        params_correct, msg = comp_params.check_params()        
        if params_correct:
            return True, (comp_params, comp_array)
        else:
            return False, msg


    def calculateVectorFieldOps(self):

        # input rasters            
        inraster_x, inraster_y = self.get_raster_parameters()

        # vector field parameters            
        magnitude_calc_choice, magnitude_outraster_path, \
        orientations_calc_choice, orientations_outraster_path, \
        divergence_calc_choice, divergence_outraster_path, \
        curlmodule_calc_choice, curlmodule_outraster_path  = self.get_vector_operator_parameters() 

        # load results                       
        load_output = self.qchkbxOutputVfopsLoadChoice.isChecked()  
 

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
            QMessageBox.critical(self, "Input x-component rasters", msg)
            return   
        else:
            comp_x_params, comp_x_array = result
                     
        success, result = self.get_raster_data(inraster_y)
        if not success:
            QMessageBox.critical(self, "Input y-component rasters", msg)
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
        
        vector_field = Grid(comp_x_params, vector_array)            
        
        # calculates vector field parameters 
        for vfp_name, vfp_choice, vfp_savefile, vfp_function in zip(params_names, params_choices, params_savefiles, params_functions):  
            if vfp_choice:
                exec "curr_fld = vector_field.%s()" % vfp_function
                exec "curr_fld.write_esrigrid('%s')" % vfp_savefile
        
        # add required layer to the map canvas - modified after RasterCalc module                
        if load_output:
            #vector field parameters
            for vfc_choice, vfc_savefile in zip(params_choices, params_savefiles):  
                if vfc_choice:
                    newLayer = QgsRasterLayer(vfc_savefile, QFileInfo(vfc_savefile).baseName())
                    QgsMapLayerRegistry.instance().addMapLayer(newLayer)
    
   
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
        load_output = self.qchkbxOutputGradsLoadChoice.isChecked()

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
            QMessageBox.critical(self, "Input x-component rasters", msg)
            return   
        else:
            comp_x_params, comp_x_array = result
                     
        success, result = self.get_raster_data(inraster_y)
        if not success:
            QMessageBox.critical(self, "Input y-component rasters", msg)
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
        vector_array[:,:,0], vector_array[:,:,1] = comp_x_array, comp_y_array
        
        vector_field = Grid(comp_x_params, vector_array)            
        
        # calculates gradients parameters 
        for vfg_name, vfg_choice, vfg_savefile, vfg_function in zip(vfgrads_names, vfgrads_choices, vfgrads_savefiles, vfgrads_functions):  
            if vfg_choice:
                exec "curr_fld = vector_field.%s()" % vfg_function
                exec "curr_fld.write_esrigrid('%s')" % vfg_savefile
    
        # add required layer to the map canvas - modified after RasterCalc module                
        if load_output:
            # gradient parameters                                              
            for vfg_choice, vfg_savefile in zip(vfgrads_choices, vfgrads_savefiles):  
                if vfg_choice:
                    newLayer = QgsRasterLayer(vfg_savefile, QFileInfo(vfg_savefile).baseName())
                    QgsMapLayerRegistry.instance().addMapLayer(newLayer)    
   
        # all done
        QMessageBox.information(self, "Gradients output", "Processings completed.")
    

    def calculateVectorFieldPathlines(self):

        # input rasters            
        inraster_x, inraster_y = self.get_raster_parameters()

        # pathline calculation parameters                    
        pathlines_input_shapefile = self.qlnedtPathlinesInShapefile.text()
        time_step = self.qlnedtPathlinesTimeStepInput.text()
        total_time = self.qlnedtPathlinesTotalTimeInput.text()
        error_max_tolerance = self.qlnedtPathlinesStepMaxErrorInput.text()
        pathlines_output_shapefile = self.qlnedtPathlinesOutShapefileInput.text()
        
        # load results                       
        load_output = self.qchkbxOutputPathlineLoadChoice.isChecked()

        ### PRE-PROCESSINGS
        
        # check input values for vector field rasters
        success, msg = self.check_inrasters_def(inraster_x, inraster_y)
        if not success:
            QMessageBox.critical(self, "Input component rasters", msg)
            return
        
        # get x- and y-axis component data
        success, result = self.get_raster_data(inraster_x)
        if not success:
            QMessageBox.critical(self, "Input x-component rasters", msg)
            return   
        else:
            comp_x_params, comp_x_array = result
                     
        success, result = self.get_raster_data(inraster_y)
        if not success:
            QMessageBox.critical(self, "Input y-component rasters", msg)
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
            QMessageBox.critical(self, "Error in pathline calculation", "Unable to get input shapefile")   
            return
   
        if os.path.exists(pathlines_output_shapefile):
            driver.DeleteDataSource(pathlines_output_shapefile)
    
        out_shape = driver.CreateDataSource(pathlines_output_shapefile)
        if out_shape is None:
            QMessageBox.critical(self, 
                                 "Pathline calculation", 
                                 "Unable to create output shapefile: %s" % pathlines_output_shapefile)   
            return
        
        out_layer = out_shape.CreateLayer('out_pathlines', geom_type=ogr.wkbPoint)
        if out_layer is None:
            QMessageBox.critical(self, 
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
            <p>VectorFieldCalc version 1.3<br />License: GPL v. 3</p>
            <p>M. Alberti, <a href="http://www.malg.eu">www.malg.eu</a></p> 
            <p>This application calculates vector field parameters (e.g., divergence, curl module, gradients)
            and pathlines.            
            </p>
             <p>Please report any bug to <a href="mailto:alberti.m65@gmail.com">alberti.m65@gmail.com</a></p>
        """)              
            
        
        
        
                    

     
            
            
            
