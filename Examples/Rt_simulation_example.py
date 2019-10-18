"""
Example rt_simulation:

In this example a Pastas model is set up and integrated in the real time 
simulation class to perform a forecast and an update. 
NOTE: To make the script work, update the working  
directory

"""


import os
import pylab
import pandas as pd
import pastas as ps
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys
from scipy.signal import fftconvolve
import math
from hydroeval import nse, kge, rmse, pbias
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from math import sqrt
import pickle



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set the location of the working directory below (location of the pastas-realtime  folder)
dr = str(r'C:\Users\Karl.schutt\Desktop\Pastas real time application')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Import tools and real time simulation class
sys.path.append(dr+'/Tools')
from Eref_to_epot import calc_epot # Function for calculating potential evaporation
sys.path.append(dr+'/real time functionalities')
from rt_simulation import Rt_simulation # Function for performing real-time forecasts


"""
Import data

"""

# Import groundwater data and edit
Pbcode = 'B58C0352'
filternumber = str(1)
gwl = ps.read_dino(dr+'/Data/Real time simulation example/Groundwater/'+Pbcode+'00'+filternumber+'_1.csv')
gwl = gwl.series
gwl = gwl.to_frame()['2001-01-01 00:00:00':]
gwl.columns = ['Stand (cm t.o.v. NAP)']
for i in range(len(gwl)):
    if gwl['Stand (cm t.o.v. NAP)'][i] < gwl['Stand (cm t.o.v. NAP)'].mean()-(3*gwl['Stand (cm t.o.v. NAP)'].std()) or gwl['Stand (cm t.o.v. NAP)'][i] > gwl['Stand (cm t.o.v. NAP)'].mean()+(3*gwl['Stand (cm t.o.v. NAP)'].std()):
        gwl['Stand (cm t.o.v. NAP)'][i] = 'nan'
gwl = gwl.interpolate()

# Distinguish historical data and year with most recent observations 
gwl1 = gwl[:len(gwl)-360]
gwl2 = gwl[len(gwl)-360:]


# Import meteorological input data and edit
#knmi_stn = ps.read.KnmiStation(stns=377, interval = ' daily'  ) 
#knmi_stn.download()
#rain = knmi_stn.data['RH'].resample('D').sum()['2001-01-01 00:00:00':]
#evap_r = knmi_stn.data['EV24'].resample('D').sum()['2001-01-01 00:00:00':]

# When the ps.read.KnmiStation function does not work, import the data from files:
rain = ps.read.knmi.read_knmi(str(dr+'/Data/Real time simulation example/Meteo/Ell.txt'),variables = 'RH').series.resample('D').sum()['2001-01-01 00:00:00':]
evap_r = ps.read.knmi.read_knmi(str(dr+'/Data/Real time simulation example/Meteo/Ell.txt'),variables = 'EV24').series.resample('D').sum()['2001-01-01 00:00:00':]
cf = pd.DataFrame(index = [1,2,3,4,5,6,7,8,9,10,11,12])
cf['cf'] = [0.9,0.9,0.9,1,1,1,1,1,0.9,0.9,0.9,0.9]
evap_p = calc_epot(cf,evap_r) # Calculate potential evaporation with crop factors
recharge = rain-evap_p

# Distinguish historical data and most recent year with observations 
rain1 = rain.resample('D').sum()
rain2 = rain.resample('D').sum()
evap1_r = evap_r.resample('D').sum()
evap2_r = evap_r.resample('D').sum()
evap1_p = evap_p.resample('D').sum()
evap2_p = evap_p.resample('D').sum().loc[gw_val.index,]

# Calculate precipitation excess
recharge1 = recharge[ : gw.index[len(gw.index)-1]]
recharge2 = recharge[gw_val.index[0] : gw_val.index[len(gw_val.index)-1]]


"""
Create initial Pastas model with historical dataset

"""

ml = ps.Model(gwl1, name = "GWL") # Create the Pastas model
ts1 = ps.StressModel(recharge1, ps.FourParam, name = 'recharge') # create stressmodel for precipitation excess *
ml.add_stressmodel(ts1) # Add stressmodel to Pastas model
ml.solve() # Solve the Pastas model

# Analyse results
ml.stats.summary()
ml.plots.results()



"""
Real time simulation

"""
# Create real-time simulation object
rtc = Rt_simulation(ml,"test")

# Forecasts of the groundwater level for a year ahead based on observed precipitation excess: 
f=recharge2.to_frame() # Create dataframe with forcing data
f.columns = ['recharge'] # Name of stress should correspond with the name of the stressmodel in the Pastas model *
forecast = rtc.forecast(f).to_frame() # Perform the forecast

# Plot forecast
forecast['Simulation'].plot()
gwl2['Stand (cm t.o.v. NAP)'].plot()
plt.legend(labels = ['forecasted','observed'])


# Update the Pastas model with a year of newly observed data
cwlen = 'total' # Length of the calibration period. (In this case the total set of historical data, but can also specified as a number of days)
rtc.update(f,gwl2,cwlen) # Update the model (give the observed forcing data, observed groundwater level data and calibration period length as arguments)
rtc.model.plot() # Plot the results
