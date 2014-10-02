# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 19:44:52 2014

@author: joseph
"""
from Soil3 import * #import tool classes
from TrainingSet3 import TrainingSet

import numpy as np

#from guiqwt.plot import CurveDialog


import guidata
guidata.qapplication() # not required if a QApplication has already been created

import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di

Set=TrainingSet(50)

class Processing(dt.DataSet):
    """Example"""
    
    
    chosensoil1 = di.ChoiceItem("Soil Name (blue line)",Set.NameList)
    chosensoil2 = di.ChoiceItem("Soil Name (red line)",Set.NameList)
    

    def SpectrumVisu(self,item,value,parent):
        #app=guidata.qapplication()
        from guiqwt.curve import CurvePlot
        plot=CurvePlot(title="spectrum",xlabel="x",ylabel="y")

        number1=self.chosensoil1
        spectrum1=Set.SpectrumList[number1]
        x1=range(len(spectrum1))
        number2=self.chosensoil2
        spectrum2=Set.SpectrumList[number2]
        x2=range(len(spectrum2))

        from guiqwt.builder import make
        
        curve1=make.curve(x1,spectrum1,title="curve1",color='b')
        plot.add_item(curve1)
        curve2=make.curve(x2,spectrum2,title="curve2",color='r')
        plot.add_item(curve2)
        #app.exec_()
        plot.show()
        print(self.chosensoil1)
        
   
        
    button=di.ButtonItem("launch method",SpectrumVisu)
        
   
param = Processing()
param.edit()

