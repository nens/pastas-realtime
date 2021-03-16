"""
Description: 
    
This file contains a class that can be used for controlling a 
Pastas model. It can be useful for real-time simulations. The class is 
basically a wrapper around a Pastas model which can be used for updating the 
model with newly observed data (oseries and stresses) and performing forecasts 
based on timeseries containing expected values of the stresses in the future. 

It should be noted that this class can only handle the following response 
functions at the moment: Gamma, Hantush, Exponential, One, FourParam, 
DoubleExponential and Polder. The script does not work yet for a linear 
trend.

Author: Karl Schutt

Last time editted: 16-10-2019
"""

import sys
import pandas as pd
import pastas as ps
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

class rt_simulation:
    """
    Wrapper around Pastas model for updating and forecasting which can be 
    used for real-time simulations
    """   
    def __init__(self, model, name):
        """ 
        Initialize a real-time simulation. 
        
        Parameters
        ----------
        model: (Pastas model)
            A solved Pastas model 
        name: (str)
            Name of the simulation
        
        """
        self.name = name
        self.model = model
        self.stressmodelnames = self.model.get_stressmodel_names()
        self.pars = self.model.get_parameters()
        self.residuals = self.model.residuals()
        self.std_residuals = np.std(np.absolute(self.residuals))
        
        # Get the type of response functions of the model 
        self.stressmodeltypes = []
        for i in range(len(model.stressmodels)):
            name = self.model.get_stressmodel_names()[i]
            self.stressmodeltypes.append(self.model.stressmodels[name].rfunc._name)
    
    def update(self,forcing, oseries, selection):
        """ 
        Update the model with newly observed data
        
        Parameters
        ----------
        forcing: Pandas dataframe
            Dataframe that holds the newly observed values of the stresses. 
            NOTE: column names should correspond to the names of the stresses
            in the stressmodels of the Pastas model. The index column should consist
            of the dates of observations (Timestamps). The first date should succeed 
            the last date of the oseries & forcing dataset that are currently 
            included in the Pastas model.
        oseries: Pandas dataframe
            Dataframe that holds the newly observed values of the oseries. It
            should consist of one column containing the oseries.
            NOTE: the column name should equal the column name of the oseries 
            in the Pastas model. The index should consist of the dates of 
            observation (Timestamps). The first date should succeed the last date 
            of the oseries & forcing dataset that are currently included in the 
            Pastas model.
        """

        # Update forcing data
        forc = forcing.copy()
        initial_data = self.model.get_stress(self.stressmodelnames[0])[:].to_frame()
        initial_data.columns = [self.stressmodelnames[0]]
        for i in range(len(list(forcing.columns))-1):
            initial_data[self.stressmodelnames[i+1]] = self.model.get_stress(self.stressmodelnames[i+1])[:]
        forc = initial_data.append(forc[:])
        
        # Update oseries 
        new_oseries = self.model.oseries.series.to_frame().append(oseries[:])
        if selection != 'total': 
            new_oseries = new_oseries[new_oseries.index[len(new_oseries)-1]-timedelta(days = selection):]
        
        # Create new model
        new_model = ps.Model(new_oseries, name = "forecast")

        # Add stressmodels
        for i in range(len(self.stressmodelnames)):
            # The forecast function automaticly infers the response functions from the stressmodels
            if self.stressmodeltypes[i] == 'Gamma':
                responsefunc = ps.Gamma
            elif self.stressmodeltypes[i] == 'Hantush':
                responsefunc = ps.Hantush                      
            elif self.stressmodeltypes[i] == 'Exponential':
                responsefunc = ps.Exponential
            elif self.stressmodeltypes[i] == 'One':
                responsefunc = ps.One
            elif self.stressmodeltypes[i] == 'FourParam':
                responsefunc = ps.FourParam
            elif self.stressmodeltypes[i] == 'DoubleExponential':
                responsefunc = ps.DoubleExponential
            elif self.stressmodeltypes[i] == 'Polder':
                responsefunc = ps.Polder
            ts = ps.StressModel(forc[self.stressmodelnames[i]], responsefunc, name=self.stressmodelnames[i])
            
            # Add the stressmodels to the new model
            new_model.add_stressmodel(ts)
        
            if self.model.constant:
                c = ps.Constant()
                new_model.add_constant(c)
            if self.model.noisemodel:
                n = ps.NoiseModel()
                new_model.add_noisemodel(n)
            if self.model.transform:
                tt = ps.ThresholdTransform()
                new_model.add_transform(tt)
            
        # Solve the new model and initialise the real-time simulation overnew 
        new_model.solve()
        self.__init__(new_model,self.name)
                 
    def forecast(self,forcing):
        """ 
        Perform a forecast with timeseries containing the expected future values
        of the stresses. The timeseries of stresses in the Pastas model are 
        temporarly complemented with future values and a simulation is performed
        using the current parameter settings.
        
        Parameters:
        -----------
        forcing: Pandas Dataframe
            Dataframe that holds expected values of the stresses in the future. 
            NOTE: column names should correspond to the names of the stresses
            in the stressmodels of the Pastas model. The index column should consist
            of the dates (Timestamps). The first date should succeed the last date 
            of the oseries & forcing dataset that are currently included in the 
            Pastas model.
            
        Returns: Dataframe with forecasted values 
            
        """
        
        # Check wether stresses in forcing dataframe correspond with the stresses in the model
        if self.stressmodelnames != list(forcing.columns):
            return ("stresses in forcing data do not correspond with stresses in stressmodels")
        
        # Create forecasting dataset (stresses)
        initial_data = self.model.get_stress(self.stressmodelnames[0])[:].to_frame()
        initial_data.columns = [self.stressmodelnames[0]]
        for ii in range(len(list(forcing.columns))-1):
            initial_data[self.stressmodelnames[ii+1]] = self.model.get_stress(self.stressmodelnames[ii+1])[:]
        forc = initial_data.append(forcing[:])
         
        # Create a temporary model for forecasting
        forecast_ml = ps.Model(self.model.oseries, name = "forecast")    
        for i in range(len(self.stressmodelnames)):
            # The forecast function automaticly infers the response functions from the stressmodels
            if self.stressmodeltypes[i] == 'Gamma':
                responsefunc = ps.Gamma
            elif self.stressmodeltypes[i] == 'Hantush':
                responsefunc = ps.Hantush                      
            elif self.stressmodeltypes[i] == 'Exponential':
                responsefunc = ps.Exponential
            elif self.stressmodeltypes[i] == 'One':
                responsefunc = ps.One
            elif self.stressmodeltypes[i] == 'FourParam':
                responsefunc = ps.FourParam
            elif self.stressmodeltypes[i] == 'DoubleExponential':
                responsefunc = ps.DoubleExponential
            elif self.stressmodeltypes[i] == 'Polder':
                responsefunc = ps.Polder
            ts = ps.StressModel(forc[self.stressmodelnames[i]], responsefunc, name=self.stressmodelnames[i])
            # The created stressmodels are added to the forecast model
            forecast_ml.add_stressmodel(ts)   
        if self.model.constant:
            forecast_ml.constant = self.model.copy().constant
        if self.model.noisemodel:
            forecast_ml.noisemodel = self.model.copy().noisemodel
        if self.model.transform:
            forecast_ml.transform = self.model.copy().transform
                
        # Perform a simulation with the forecast model and return the forecasted values of the oseries
        forecast = forecast_ml.simulate(parameters = self.p)
        return(forecast[len(forc)-len(forcing):])
