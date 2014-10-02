# -*- coding: utf-8 -*-
"""
Created on Thu Oct 02 13:15:42 2014

@author: joseph
"""
from TrainingSet4 import*
from LearningAlgo2 import*
import numpy as np


Set50=TrainingSet(20)
Result=Learning(Set50)


Set100=TrainingSet(40)

EstimateList=[]
for soil in Set100.SoilList:
    EstimateList.append(Result.Estimate(Set50,soil))

plt.figure()
plt.plot(EstimateList,'ro')
plt.plot(Set100.pHList,'go')
plt.xlabel('Soil number')
plt.ylabel('pH')
plt.legend(['Prediction before learning','Real values'])

Result.Launch()

EstimateList2=[]
for soil in Set100.SoilList:
    EstimateList2.append(Result.Estimate(Set50,soil))
    
plt.figure()
plt.plot(EstimateList2,'bo')
plt.plot(Set100.pHList,'go')
plt.xlabel('Soil number')
plt.ylabel('pH')
plt.legend(['Prediction after learning','Real values'])