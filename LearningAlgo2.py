# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 17:42:33 2014

@author: joseph
"""

from Soil3 import * #import tool classes
from TrainingSet3 import *
import numpy as np
from numpy import linalg as la
from scipy.optimize import minimize

class Learning:
    """Learning algorithm based on the chosen training set.
    The idea is to find a topology that seperates the graph of soils (where each 
    edge is weighted by the formal distance between the vertices)
    into two distinct components corresponding to soils having pH>0 and pH<0.
    
    The formal distance "d(m,n)" between soils m and n is given by:
    d(m,n)=param1*d1(m,n)+...+paramN*dN(m,n)
    where d1(), ..., d2() are distances between soils features.    
    """
    
    def __init__(self,TrainingSet):
        """ The attributes of the "Learning" class are:
        - Param : Parameters to train
        - M0 : Matrix of inter-soils formal distances we idealy want to obtain from 
        the data thanks to the trained parameters
        - M1, ...,M5 : features distances matrices
        """
        self.Param=[0.2,0.2,0.2,0.2,0.2] #set of parameters to train
        
        m=TrainingSet.NumberOfElements
        M0=np.zeros((m,m))
        M1=np.zeros((m,m))
        M2=np.zeros((m,m))
        M3=np.zeros((m,m))
        M4=np.zeros((m,m))
        M5=np.zeros((m,m))
        
             
        for i in range(0,m):
            for j in range(i,m):
                if TrainingSet.pHList[i]*TrainingSet.pHList[j]<0: #if pH values of the two selcted soils have opposite signs then the distance is set to 1, and if not to 0.
                    M0[i,j]=10
                M1[i,j]=TrainingSet.SpectrumDistance2(i,j,'band1')
                M2[i,j]=TrainingSet.SpectrumDistance2(i,j,'band2')
                M3[i,j]=TrainingSet.SpectrumDistance2(i,j,'band3')
                M4[i,j]=TrainingSet.SpectrumDistance2(i,j,'band4')
                M5[i,j]=TrainingSet.SpectrumDistance2(i,j,'band5')

        
        self.MatrixSize=m
        self.MatrixGoal=np.matrix(M0)
        self.MatrixDist1=np.matrix(M1)
        self.MatrixDist2=np.matrix(M2)
        self.MatrixDist3=np.matrix(M3)
        self.MatrixDist4=np.matrix(M4)
        self.MatrixDist5=np.matrix(M5)
        #self.AllMatrices=[np.matrix(M0),np.matrix(M1),np.matrix(M2),np.matrix(M3),np.matrix(M4),np.matrix(M5)]
        self.AllMatrices=[self.MatrixGoal/la.norm(self.MatrixGoal),self.MatrixDist1/la.norm(self.MatrixDist1),self.MatrixDist2/la.norm(self.MatrixDist2),self.MatrixDist3/la.norm(self.MatrixDist3),self.MatrixDist4/la.norm(self.MatrixDist4),self.MatrixDist5/la.norm(self.MatrixDist5)]
        
    def Launch(self):
        """ The idea is to use a constrained optimization procedure to obtain the best
        formal distance that will separates soils with pH>0 from soils with pH<0.
        For this we will "train" the parameters of the formal distance.
        """
        
        #Initialisation stage
        
        M=self.AllMatrices #Inter-feature distances matrices
        Param=self.Param
        

        # A and B will be used to compute de gradient of the function to minimize
        a=np.ones((5,6))
        B=np.matrix(np.ones((5,1)))
        for p in range(1,6):
            for k in range(0,6):
                a[p-1,k]=2*np.trace(M[k]*M[p].T) 
            B[p-1]=-a[p-1,0]
        A=np.matrix(a[:,1:])
        
        
        #Training stage
        
        def func(x): #function to minimize      
            return la.norm(M[0]-(x[0]*M[1]+x[1]*M[2]+x[2]*M[3]+x[3]*M[4]+x[4]*M[5]),'fro')**2 
            
        def func_grad(x): #gradient of the function to minimize    
            xm=np.matrix(x)
            G=A*xm.transpose()+B
            return np.array(G.transpose()) 
            
        #constraints on the minimization (sum of the parameters equal 1 et all positive)
        cons=(
        {'type':'eq',
         'fun': lambda x: np.array([x[0]+x[1]+x[2]+x[3]+x[4]-1]),
        'jac': lambda x: np.array([1.,1.,1.,1.,1.])},
        {'type':'ineq',
         'fun': lambda x: np.array([x[0]]),
        'jac': lambda x: np.array([1.,0.,0.,0.,0.])},
        {'type':'ineq',
         'fun': lambda x: np.array([x[1]]),
        'jac': lambda x: np.array([0.,1.,0.,0.,0.])},
        {'type':'ineq',
         'fun': lambda x: np.array([x[2]]),
        'jac': lambda x: np.array([0.,0.,1.,0.,0.])},
        {'type':'ineq',
         'fun': lambda x: np.array([x[3]]),
        'jac': lambda x: np.array([0.,0.,0.,1.,0.])},
        {'type':'ineq',
         'fun': lambda x: np.array([x[4]]),
        'jac': lambda x: np.array([0.,0.,0.,0.,1.])}
        )
        
        print(func(Param))
        res=minimize(func,Param,jac=func_grad,constraints=cons,tol=1e-6,method='SLSQP',options={'disp':True})   
        self.Param=res.x        
        print(res.x)
        #print(func(res.x))
        Mfinal=res.x[0]*M[1]+res.x[1]*M[2]+res.x[2]*M[3]+res.x[3]*M[4]+res.x[4]*M[5]
        
        return Mfinal, max(np.array(abs(M[0]-Mfinal)).flat)

    def Visualize(self,TrainingSet,SoilNumber):
        """Plot the formal distance between a given Soil (correponding to SoilNumber) 
        and all the other soils of the training set, versus their pH values
        Outputs:
            -tSame: list of soil having pH value with the same sign than the pH value of the chosen soil
            -tDiff: list of soil having pH value with the different sign than the pH value of the chosen soil
        """
        
        Param=self.Param
        m=self.MatrixSize
        M=self.AllMatrices
        D=Param[0]*M[1]+Param[1]*M[2]+Param[2]*M[3]+Param[3]*M[4]+Param[4]*M[5]
        Dtot=D+D.transpose() #total inter-soil distances matrix
        
        pH=TrainingSet.pHDistancematrix(2)
        pHtot=pH+pH.T
        
        #Visualize distance versus pH
        figure()#2D plot to visualize distance between soils
        plt.plot(np.array(pHtot[SoilNumber]).T,np.array(Dtot[SoilNumber].T),'bo')
        xlabel('pH class')
        ylabel('Formal distance')
        
        xSame=[]
        ySame=[]
        tSame=[]
        xDiff=[]
        yDiff=[]
        tDiff=[]
        
        #2D plot to visualize distance between soils
        for p in range(m):
            if pHtot[SoilNumber,p]==1:
                tDiff.append(Dtot[SoilNumber,p])
                r=random()
                xDiff.append(Dtot[SoilNumber,p]*cos(2*pi*r))
                yDiff.append(Dtot[SoilNumber,p]*sin(2*pi*r))
            r=random()
            tSame.append(Dtot[SoilNumber,p])
            xSame.append(Dtot[SoilNumber,p]*cos(2*pi*r))
            ySame.append(Dtot[SoilNumber,p]*sin(2*pi*r))
        
        figure()
        plt.hold(True)
        plt.plot(xSame,ySame,'bo',label='same pH')
        plt.plot(xDiff,yDiff,'ro',label='opposite pH')
        legend()
        
        return tSame, tDiff
        
    def VisualizeAll(self,TrainingSet):
        """Plot the formal distance between all the soils of the training set, 
         versus their pH values"""
        
        Param=self.Param
        M=self.AllMatrices
        D=Param[0]*M[1]+Param[1]*M[2]+Param[2]*M[3]+Param[3]*M[4]+Param[4]*M[5]
        Dtot=D+D.transpose()
        
        pH=TrainingSet.pHDistancematrix(1)
        pHtot=pH+pH.T
        
        figure()
        pYplot.plot(np.array(pHtot).T,np.array(Dtot),'bo')
        xlabel('pH ')
        ylabel('Formal distance')
        
        
        
    def Estimate(self,TrainingSet,Soil):
        """Uses trained parameters to estimate the pH value of the soil"""
        
        param=self.Param
        dist=[]

        for trained in TrainingSet.SoilList:
            d=param[0]*spectrumdistance(Soil,trained,'band1')+param[1]*spectrumdistance(Soil,trained,'band2')+param[2]*spectrumdistance(Soil,trained,'band3')+param[3]*spectrumdistance(Soil,trained,'band4')+param[4]*spectrumdistance(Soil,trained,'band5')
            dist.append(d)
        index=np.argmin(dist) 
        

        return TrainingSet.SoilList[index].pH
        
        
        

  