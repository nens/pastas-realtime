"""
Pastas project real-time simulation example:

In this example a Pastas project is created and subsequently a real-time
simulation is run. For a year, forecasts are done every month and
the model is updated afterwards.
NOTE: the working directory should be editted to make the script work 

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
from sklearn.metrics import mean_absolute_error
import pickle


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set the location of the working directory below (location of the pastas-realtime  folder)
dr = str(r'C:\Users\Karl.schutt\Desktop\Pastas real time application')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


sys.path.append(dr+'/Real time functionalities')
from Project_rt import * # Pastas project aangevuld met voorspellingsmodus
sys.path.append(dr+'/Tools')
from Eref_to_epot import calc_epot # Functie om de potentieele verdamping uit te rekenen



###############################################################################
# Usefull functions for project setup / data analysis
###############################################################################

# Fill gaps in groundwater timeseries:
def fill_gaps(series): # input argument should be a pandas series object
    idx = pd.date_range(series.index[0],series.index[len(series)-1])
    filled = series.reindex(idx, fill_value = None)
    filled = filled.interpolate()
    return(filled)

# Pre edit pastas groundwatertimeseries before importation in a pastas project
def pre_edit_gw_series(series,tmin,tmax): # input argument should be a pastas timeseries
    series_filled = fill_gaps(series.series[tmin:tmax])
    series_pre_edit = ps.TimeSeries(series = series_filled, metadata = series.metadata)
    return(series_pre_edit)
    

# Function for handling the forecast data for plotting
def merge_forecasts(forc,idx,loc):
    """
    Explanation
    forc: dictionary holding dataframes with forecasts for all locations
    idx: index of the forecasts given as array
    location: list of locations for which forecasts were performed
    
    """
    forecasts = dict()
    for i in range(len(locations)):
        location = loc[i]
        loc_forc = pd.DataFrame(columns = ['Simulation'])
        for ID in range(len(idx)):
            loc_forc = loc_forc.append(forc[idx[ID]][location])
        forecasts[location]=loc_forc
    
    return(forecasts)

###############################################################################
# Import data
###############################################################################

"""
Oseries (Groundwater levels)
"""
# Create dataframe containing the names of locations, historic data (ps.timeseries) and recent data(pd.series)
# historic data needs to be ps.Timeseries because the models should infer the metadata of the historic series
Oseries = pd.DataFrame(index = ['B58C0352','B58A0093','B58B0260','B52C0005','B52E3234'])
Oseries['historic_series']=0 
Oseries['forecast_val_series']=0
for i in range(len(Oseries.index)):
    
    historic = pd.read_csv(dr+'\\Data\\Pastas project real time example\\Groundwater\\'+Oseries.index[i]+'_historic.csv',index_col=0)
    historic.index = pd.to_datetime(historic.index, format = '%Y-%m-%d')
    with open(dr+'\\Data\\Pastas project real time example\\Groundwater\\'+Oseries.index[i]+'metadata.pkl', 'rb') as handle:
            metadata = pickle.load(handle)
    historic_series = ps.TimeSeries(historic, name=Oseries.index[i], metadata=metadata)
    Oseries['historic_series'][Oseries.index[i]]=historic_series
    
    forecast_val = pd.read_csv(dr+'\\Data\\Pastas project real time example\\Groundwater\\'+Oseries.index[i]+'_recent.csv',index_col=0)
    forecast_val.index = pd.to_datetime(forecast_val.index, format = '%Y-%m-%d')
    forecast_val_series = ps.TimeSeries(forecast_val, name=Oseries.index[i], metadata=metadata)
    Oseries['forecast_val_series'][Oseries.index[i]]=forecast_val_series


"""
Stresses (Precipitation, Evaporation)
"""

# Import precipitation, and define metadata:

Projection = 'epsg:28992'
Pseries = pd.DataFrame(index = ['Ell_p','Arcen_p','Ysselsteyn_p','Heibloem_p'])
Pseries['X']=['181262.123899','211002.705853','189993.000000','190001.000000']
Pseries['Y']=['356393.981379','390542.713609','389000.000000','367002.000000']
Pseries['historic_series']=0 
Pseries['forecast_val_series']=0
for i in range(len(Pseries.index)):
    
    historic = pd.read_csv(dr+'\\Data\\Pastas project real time example\\Meteo\\'+Pseries.index[i]+'_historic.csv',index_col=0)
    historic.index = pd.to_datetime(historic.index, format = '%Y-%m-%d')

    # Metadata is added to the Timeseries
    historic_series = ps.TimeSeries(historic, name=Pseries.index[i])
    historic_series.metadata['x']=Pseries['X'][i]
    historic_series.metadata['y']=Pseries['Y'][i]
    historic_series.metadata['projection']=Projection
    Pseries['historic_series'][Pseries.index[i]]=historic_series
    
    forecast_val = pd.read_csv(dr+'\\Data\\Pastas project real time example\\Meteo\\'+Pseries.index[i]+'_recent.csv',index_col=0)
    forecast_val.index = pd.to_datetime(forecast_val.index, format = '%Y-%m-%d')
    forecast_val_series = ps.TimeSeries(forecast_val, name=Pseries.index[i])
    forecast_val_series.metadata['x']=Pseries['X'][i]
    forecast_val_series.metadata['y']=Pseries['Y'][i]
    forecast_val_series.metadata['projection']=Projection
    Pseries['forecast_val_series'][Pseries.index[i]]=forecast_val_series

# Import Evaporation (In this case reference evaporation is used)
Projection = 'epsg:28992'
ETseries = pd.DataFrame(index = ['Ell_et','Arcen_et'])
ETseries['X']=['181262.123899','211002.705853']
ETseries['Y']=['356393.981379','390542.713609']
ETseries['historic_series']=0 
ETseries['forecast_val_series']=0
for i in range(len(ETseries.index)):
    
    historic = pd.read_csv(dr+'\\Data\\Pastas project real time example\\Meteo\\'+ETseries.index[i]+'_historic.csv',index_col=0)
    historic.index = pd.to_datetime(historic.index, format = '%Y-%m-%d')
    
    # Metadata is added to the Timeseries
    historic_series = ps.TimeSeries(historic, name=ETseries.index[i])
    historic_series.metadata['x']=ETseries['X'][i]
    historic_series.metadata['y']=ETseries['Y'][i]
    historic_series.metadata['projection']=Projection
    ETseries['historic_series'][ETseries.index[i]]=historic_series
    
    forecast_val = pd.read_csv(dr+'\\Data\\Pastas project real time example\\Meteo\\'+ETseries.index[i]+'_recent.csv',index_col=0)
    forecast_val.index = pd.to_datetime(forecast_val.index, format = '%Y-%m-%d')
    forecast_val_series = ps.TimeSeries(forecast_val, name=ETseries.index[i])
    forecast_val_series.metadata['x']=ETseries['X'][i]
    forecast_val_series.metadata['y']=ETseries['Y'][i]
    forecast_val_series.metadata['projection']=Projection
    ETseries['forecast_val_series'][ETseries.index[i]]=forecast_val_series    

  
# Create dataframe with stresses
Stresses = Pseries.append(ETseries)
# Define the kind of the stresses
Stresses['kind']=['prec','prec','prec','prec','evap','evap']

# Create monthly forecast sets:
start = ['2015-01-01','2015-02-01','2015-03-01','2015-04-01','2015-05-01','2015-06-01',
         '2015-07-01','2015-08-01','2015-09-01','2015-10-01','2015-11-01','2015-12-01']
end = ['2015-01-31','2015-02-28','2015-03-31','2015-04-30','2015-05-31','2015-06-30',
         '2015-07-31','2015-08-31','2015-09-30','2015-10-31','2015-11-30','2015-12-31']
month = ['jan','feb','mar','apr','may','jun','jul',
         'aug','sep','okt','nov','dec']

for i in range(12):
    title = month[i]+'_calseries'
    Oseries[title]=0
    for name in Oseries.index:
        Selection=Oseries['forecast_val_series'][name].series[start[i]:end[i]]
        Metadata = Oseries['forecast_val_series'][name].metadata
        s_name = Oseries['forecast_val_series'][name].name 
        Oseries[title][name]=ps.TimeSeries(Selection, metadata = Metadata, name = s_name)

for i in range(12):
    title = month[i]+'_forcseries'
    Stresses[title]=0
    for name in Stresses.index:
        Selection=Stresses['forecast_val_series'][name].series[start[i]:end[i]]
        Metadata = Stresses['forecast_val_series'][name].metadata
        s_name = Stresses['forecast_val_series'][name].name 
        Stresses[title][name]=ps.TimeSeries(Selection, metadata = Metadata, name = s_name)

###############################################################################
# Set up the project
###############################################################################


# Create project:
pr = Project(name = 'Test')

# Import the observational groundwater timeseries in the Pastas project
for i in range(len(Oseries.index)):
    pr.add_series(Oseries['historic_series'][i],kind='oseries') 

# Import the precipitation timeseries in the Pastas project
for i in range(len(Stresses.index)):
    pr.add_series(Stresses['historic_series'][i],kind=Stresses['kind'][i], settings = Stresses['kind'][i])     
        
# Plot the locations of the monitoring locations
f,ax= plt.subplots()
ax.axis('equal')
pr.maps.series(kind='oseries')
pr.maps.series(kind='prec')
pr.maps.series(kind='evap')
pr.maps.series(kind='waterlevel')
pr.maps.series(kind='flux')

# Create models 
for name in pr.oseries.index:
    ml = pr.add_model(name)
    pr.add_recharge(ml,rfunc = ps.FourParam)
    ml.solve(report=False)
    ml.plots.results()




###############################################################################
# Perform forecasts with 'future' forcing data
###############################################################################

"""
Real-time forecasts: 
The forecasts are performed for the year 2015. Each month, a forecast is done
for the groundwater level of one month ahead and subsequently the model is 
updated with observed data. The continuous forecasts are plotted for each 
location
"""

Forecast_data = Stresses.drop(['X','Y','historic_series','forecast_val_series','kind'],axis = 1)
Cal_Oseries= Oseries.drop(['historic_series','forecast_val_series'], axis = 1)

fdata_true = dict()

# Real-time forecasts
for i in range(len(Forecast_data.columns)):
    forc = pd.DataFrame(Forecast_data[Forecast_data.columns[i]])
    forc['kind'] = 0
    forc.columns = ['true','kind']
    
    # Ad the predictions to the dictionaries for the different scenario's
    forc_true = pd.DataFrame(forc['true'])
    forc_true.columns = ['series']
    forc_true['kind'] = forc['kind']
    pred_true =  pr.forecast(forc_true)
    fdata_true[month[i]]=pred_true

    # Import the observational groundwater timeseries in the Pastas project
    for u in range(len(Oseries.index)):
        name = Oseries['historic_series'][u].name
        metadata =  Oseries['historic_series'][u].metadata
        pr.del_oseries(name)
        series_his = Oseries['historic_series'][u].series
        series_forc = Cal_Oseries[Cal_Oseries.columns[i]][u].series
        new_series = ps.TimeSeries(series_his.append(series_forc),name = name, metadata = metadata)
        pr.add_series(new_series,kind='oseries') 
        Oseries['historic_series'][u]=new_series
    
    # Import the precipitation timeseries in the Pastas project
    for u in range(len(Stresses.index)):
        kind = Stresses['kind'][u]
        name = Stresses['historic_series'][u].name
        metadata =  Stresses['historic_series'][u].metadata
        pr.del_stress(name)
        series_his = Stresses['historic_series'][u].series
        series_forc = Forecast_data[Forecast_data.columns[i]][u].series
        new_series = ps.TimeSeries(series_his.append(series_forc),name = name, metadata = metadata)
        pr.add_series(new_series,kind=kind, settings = kind) 
        Stresses['historic_series'][u]=new_series
    
    # Update de model series()
    pr.update_model_series()     
    
    # Recalibrate the model
    for name in pr.oseries.index:
        ml = pr.add_model(name)
        pr.add_recharge(ml,rfunc = ps.FourParam)
        ml.solve(report=False)
            
 
                                  

###############################################################################
# Analyze results
###############################################################################        

# Define locations and months of forecasts
locations = ['B58C0352','B58A0093','B58B0260','B52C0005','B52E3234']

# Merge the monthly predicitions for each location:
forc_true_data =  merge_forecasts(forc = fdata_true, idx=month,loc=locations)


    
# Plot the results:
fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)

col = 'dodgerblue'

ax.plot(Oseries['forecast_val_series'][locations[0]].series,color = 'black', marker = '.', ms = '2', ls = '', label = 'observed')
ax.plot(forc_true_data[locations[0]],color = col, label = 'pred')

ax.title.set_text(locations[0])

ax2.plot(Oseries['forecast_val_series'][locations[1]].series,color = 'black', marker = '.', ms = '2', ls = '', label = 'observed')
ax2.plot(forc_true_data[locations[1]],color = col, label = 'pred')

ax2.title.set_text(locations[1])

ax3.plot(Oseries['forecast_val_series'][locations[2]].series,color = 'black', marker = '.', ms = '2', ls = '', label = 'observed')
ax3.plot(forc_true_data[locations[2]],color = col, label = 'pred')

ax3.title.set_text(locations[2])

ax4.plot(Oseries['forecast_val_series'][locations[3]].series,color = 'black', marker = '.', ms = '2', ls = '', label = 'observed')
ax4.plot(forc_true_data[locations[3]],color = col, label = 'pred')
 
ax4.title.set_text(locations[3])

ax5.plot(Oseries['forecast_val_series'][locations[4]].series,color = 'black', marker = '.', ms = '2', ls = '', label = 'observed')
ax5.plot(forc_true_data[locations[4]],color = col, label = 'pred')


ax5.title.set_text(locations[4])
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
L = plt.legend( loc='center left', bbox_to_anchor=(1.5, 0.5))
L.get_texts()[0].set_text('Observed')
L.get_texts()[1].set_text('Forecasted')

