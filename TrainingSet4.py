# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 22:44:05 2014

@author: joseph
"""

from Soil3 import * #import tool classes
import numpy as np

class TrainingSet:
    """Class defining an African Soil training set """
            
    def __init__(self,Nlines):
        
        # filename='C:\\Users\\joseph\\Desktop\\Kaggle\\train\\training2.csv'
        filename='C:\\Users\\joseph\\Desktop\\Kaggle\\project\\training2.csv'
        

        #initializing list
        Soil_list=[] #list of the soils contained in the 'training set'
        name_list=[]
        spectrum_list=[]        
        pH_list=[]
        pH_list_binary=[]
        BSAN_list=[]
        BSAS_list=[]
        BSAV_list=[]
        CTI_list=[]
        ELEV_list=[]
        EVI_list=[]
        

       
        with open(filename,'r') as f: #open the csv file
            line0=f.readline().split(',') #the first line contains the name of the features and properties of each soil 
            for k in range(1,Nlines+1):
                
                line=f.readline().split(',') #features and properties of the current soil
                name=line[0]
                spectrum=np.asarray([float(x) for x in line[1:3578]])
                spectrum_list.append(spectrum)
                BSAN=float(line[-21])
                BSAS=float(line[-20])
                BSAV=float(line[-19])
                CTI=float(line[-18])
                ELEV=float(line[-17])
                EVI=float(line[-16])
                LSTD=float(line[-15])
                LSTN=float(line[-14])
                REF1=float(line[-13])
                REF2=float(line[-12])
                REF3=float(line[-11])
                REF7=float(line[-10])
                RELI=float(line[-9])
                TMAP=float(line[-8])
                TMFI=float(line[-7])
                Depth=line[-6]
                Ca=float(line[-5])
                P=float(line[-4])
                pH=float(line[-3])
                SOC=float(line[-2])
                Sand=float(line[-1])
                
                name_list.append(name)
                Soil_list.append(Soil(line[0],spectrum,BSAN,BSAS,BSAV,CTI,ELEV,EVI,LSTD,LSTN,REF1,REF2,REF3,REF7,RELI,TMAP,TMFI,Depth,Ca,P,pH,SOC,Sand)) #
                BSAN_list.append(BSAN)
                BSAS_list.append(BSAS)
                BSAV_list.append(BSAV)
                CTI_list.append(CTI)
                ELEV_list.append(ELEV)
                EVI_list.append(EVI)
                pH_list.append(pH)
                
                if pH>0:
                    pH_list_binary.append(1)
                else:
                    pH_list_binary.append(-1)
                    
        self.NumberOfElements=Nlines
        self.NameList=name_list
        self.SoilList=Soil_list #list of all Soil objects
        self.ColumnNames=line0 #Names of the features and properties of Soils
        self.BSANList=BSAN_list
        self.BSASList=BSAS_list
        self.BSAVList=BSAV_list
        self.CTIList=CTI_list
        self.ELEVList=ELEV_list
        self.EVIList=EVI_list        
        self.pHList=pH_list #List of the pH values
        self.pHListBinary=pH_list_binary #List of the binary pH values
        self.SpectrumList=spectrum_list #List of the spectra
       
    
    def SpectrumDistance(self,num1,num2,band='band0'):
        """Compute distance between the two spectra of SpectrumList[num1]
        and SpectrumList[num2] in the given "wavelength" band:
        'band0' = full specrum
        'band1' = 0-800
        'band2' = 801-1600
        'band3' = 1601-2600
        'band4' = 2601-3200
        'band5' = 3201-3577 """
        
        if band=='band0':
            s1=sum(abs(self.SpectrumList[num1]-self.SpectrumList[num2]))/float(len(self.SpectrumList[num1]))
        elif band=='band1':
            s1=sum(abs(self.SpectrumList[num1][0:800]-self.SpectrumList[num2][0:800]))/float(len(self.SpectrumList[num1][0:800]))
        elif band=='band2':
            s1=sum(abs(self.SpectrumList[num1][801:1600]-self.SpectrumList[num2][801:1600]))/float(len(self.SpectrumList[num1][801:1600]))
        elif band=='band3':
            s1=sum(abs(self.SpectrumList[num1][1601:2600]-self.SpectrumList[num2][1601:2600]))/float(len(self.SpectrumList[num1][1601:2600]))
        elif band=='band4':
            s1=sum(abs(self.SpectrumList[num1][2601:3200]-self.SpectrumList[num2][2601:3200]))/float(len(self.SpectrumList[num1][2601:3200]))
        elif band=='band5':
            s1=sum(abs(self.SpectrumList[num1][3201:3577]-self.SpectrumList[num2][3201:3577]))/float(len(self.SpectrumList[num1][3201:3577]))
        return s1
        
    def SpectrumDistance2(self,num1,num2,band='band0'):
        """Compute distance between the two spectra of SpectrumList[num1]
        and SpectrumList[num2] (normalized to unit vectors) in the given "wavelength" band:
        'band0' = full specrum
        'band1' = 0-800
        'band2' = 801-1600
        'band3' = 1601-2600
        'band4' = 2601-3200
        'band5' = 3201-3577 """
        
        if band=='band0':
            s1=sum(abs(self.SpectrumList[num1]-self.SpectrumList[num2]))/(sum(abs(self.SpectrumList[num1]))*sum(abs(self.SpectrumList[num2])))
        elif band=='band1':
            s1=sum(abs(self.SpectrumList[num1][0:800]-self.SpectrumList[num2][0:800]))/(sum(abs(self.SpectrumList[num1][0:800]))*sum(abs(self.SpectrumList[num2][0:800])))
        elif band=='band2':
            s1=sum(abs(self.SpectrumList[num1][801:1600]-self.SpectrumList[num2][801:1600]))/(sum(abs(self.SpectrumList[num1][801:1600]))*sum(abs(self.SpectrumList[num2][801:1600])))
        elif band=='band3':
            s1=sum(abs(self.SpectrumList[num1][1601:2600]-self.SpectrumList[num2][1601:2600]))/(sum(abs(self.SpectrumList[num1][1601:2600]))*sum(abs(self.SpectrumList[num2][1601:2600])))
        elif band=='band4':
            s1=sum(abs(self.SpectrumList[num1][2601:3200]-self.SpectrumList[num2][2601:3200]))/(sum(abs(self.SpectrumList[num1][2601:3200]))*sum(abs(self.SpectrumList[num2][2601:3200])))
        elif band=='band5':
            s1=sum(abs(self.SpectrumList[num1][3201:3577]-self.SpectrumList[num2][3201:3577]))/(sum(abs(self.SpectrumList[num1][3201:3577]))*sum(abs(self.SpectrumList[num2][3201:3577])))
        return s1
    
    def SpectrumDistanceList(self,band='band0'):
        """Compute all the inter-spectra distances with specified "wavelength" band"""
        SDL=[]
        for n in range(0,len(self.SpectrumList)):
            for m in range(n+1,len(self.SpectrumList)):
                SDL.append(self.SpectrumDistance(n,m,band))
        return SDL
    
    def pHDistanceList(self,option=1):
        """Compute all the inter-pH distances"""
        PDL=[]
        for n in range(0,len(self.pHList)):
            for m in range(n+1,len(self.pHList)):
                if option==1:
                    PDL.append(abs(self.pHList[n]-self.pHList[m]))
                else:
                    PDL.append(0.5*abs(self.pHListBinary[n]-self.pHListBinary[m]))
                
        return PDL
        
    def pHDistancematrix(self,option=1):
        """Compute all the inter-pH distances"""
        N=len(self.pHList)
        PDL=np.zeros((N,N))
        for n in range(0,N):
            for m in range(n+1,N):
                if option==1:
                    PDL[n,m]=abs(self.pHList[n]-self.pHList[m])
                else:
                    PDL[n,m]=0.5*abs(self.pHListBinary[n]-self.pHListBinary[m])
                
        return np.matrix(PDL)
        
                
            
        

        
        
        
