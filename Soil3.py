# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 09:19:32 2014

@author: joseph
"""
import matplotlib.pyplot as pyplot
import numpy as np



class Soil:
    """Class defining an African Soil characterized by:
    - name (unique identifier PIDN)
    - spectrum (3578 points from 7497.96 cm-1 to 599.76 cm-1)
    (BSAN,BSAS,BSAV,CTI,ELEV,EVI,LSTD,LSTN,REF1,REF2,REF3,REF7,RELI,TMAP,TMFI)
    - Depth
    - Ca
    - P
    - pH
    - SOC
    - Sand """
    
    def __init__(self,name,spectrum,BSAN,BSAS,BSAV,CTI,ELEV,EVI,LSTD,LSTN,REF1,REF2,REF3,REF7,RELI,TMAP,TMFI,depth,Ca,P,pH,SOC,Sand): #constructeur
        self.name=name
        self.spectrum=spectrum
        self.spectrum_norm=np.divide(spectrum,sum(spectrum)) #normalized spectrum
        self.depth=depth
        self.Ca=Ca
        self.P=P
        self.pH=pH
        self.SOC=SOC
        self.Sand=Sand
        self.TMFI=TMFI
        self.TMAP=TMAP
        self.RELI=RELI
        self.REF7=REF7
        self.REF3=REF3
        self.REF2=REF2
        self.REF1=REF1
        self.LSTN=LSTN
        self.LSTD=LSTD
        self.EVI=EVI
        self.ELEV=ELEV
        self.CTI=CTI
        self.BSAV=BSAV
        self.BSAS=BSAS
        self.BSAN=BSAN
    
    def spectrumdisplay(self,norm=False):#method displaying the spectrum
        pyplot.figure()
        if norm==False:
            pyplot.plot(self.spectrum,'bo')
        else:
            pyplot.plot(self.spectrum_norm,'bo')
            
def spectrumdistance(soil1,soil2,band='band0'):
    """Compute distance between the two spectra of SpectrumList[num1] and SpectrumList[num2] (normalized to unit vectors) in the given "wavelength" band:
    'band0' = full specrum
    'band1' = 0-800
    'band2' = 801-1600
    'band3' = 1601-2600
    'band4' = 2601-3200
    'band5' = 3201-3577 """
        
    if band=='band0':
        s1=sum(abs(soil1.spectrum-soil2.spectrum))/(sum(abs(soil1.spectrum))*sum(abs(soil2.spectrum)))
    elif band=='band1':
        s1=sum(abs(soil1.spectrum[0:800]-soil2.spectrum[0:800]))/(sum(abs(soil1.spectrum[0:800]))*sum(abs(soil2.spectrum[0:800])))
    elif band=='band2':
        s1=sum(abs(soil1.spectrum[801:1600]-soil2.spectrum[801:1600]))/(sum(abs(soil1.spectrum[801:1600]))*sum(abs(soil2.spectrum[801:1600])))
    elif band=='band3':
        s1=sum(abs(soil1.spectrum[1601:2600]-soil2.spectrum[1601:2600]))/(sum(abs(soil1.spectrum[1601:2600]))*sum(abs(soil2.spectrum[1601:2600])))
    elif band=='band4':
        s1=sum(abs(soil1.spectrum[2601:3200]-soil2.spectrum[2601:3200]))/(sum(abs(soil1.spectrum[2601:3200]))*sum(abs(soil2.spectrum[2601:3200])))
    elif band=='band5':
        s1=sum(abs(soil1.spectrum[3201:3577]-soil2.spectrum[3201:3577]))/(sum(abs(soil1.spectrum[3201:3577]))*sum(abs(soil2.spectrum[3201:3577])))
    return s1    
        
        
        
        
            
        