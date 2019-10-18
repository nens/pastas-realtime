# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 09:43:40 2019

@author: Karl.Schutt

Function for calculating reference evaporation. Cropfactors for each month
are required for the location given as and data frame with months as index, the column name does not matter. 
The reference evaporation should be given as pd.series. Returns eref as pd.series
"""


def calc_epot(cropfactors, eref):
    #Import depences
    import numpy as np
    import pandas as pd
    
    # If eref is pd.Series
    if type(eref) == pd.core.series.Series:
        eref = pd.DataFrame(eref)
    
    # set the column name of cropfactors to a decent name
    cropfactors.columns = ['cf']
    eref.columns = ['EV24']
    
    # Calculate month out of index eref
    eref['month']= eref.index.month
     
    # Create column for storing the epot
    eref['etpot']=0
    
    # Now perform the calculation
    for i in range(len(cropfactors)):
        month = str(cropfactors.index[i])
        eref['cf'+month]=np.where(eref['month']==cropfactors.index[i],cropfactors['cf'][cropfactors.index[i]],0)
        eref['etpot']= (eref['EV24']*eref['cf'+month])+eref['etpot']
    
    # Now get rit of all the excessive columns
    eref = eref['etpot']
    eref.columns = ['etpot']
    
    # Returns eref as 
    return(eref)
    