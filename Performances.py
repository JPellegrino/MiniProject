# -*- coding: utf-8 -*-
"""
Created on Thu Oct 02 13:15:42 2014

@author: joseph
"""
from TrainingSet3 import*
from LearningAlgo import*
import numpy as np


Set50=TrainingSet(50)
Result=Learning(Set50)


Set100=TrainingSet(100)

EstimateList=[]
for soil in Set100.SoilList:
    EstimateList.append(Result.Estimate(Set50,soil))

pyplot.figure()
pyplot.plot(EstimateList,'ro')
plt.plot(Set100.pHList,'go')

Result.Launch()

EstimateList2=[]
for soil in Set100.SoilList:
    EstimateList2.append(Result.Estimate(Set50,soil))
    
pyplot.figure()
pyplot.plot(EstimateList2,'bo')
plt.plot(Set100.pHList,'go')