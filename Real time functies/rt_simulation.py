# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 16:07:25 2019

@author: Karl.Schutt
"""

import sys
import pandas as pd
import pastas as ps
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

#Add path of RD_converter file directory
sys.path.append(r'C:\Users\Karl.schutt\Documents\Tijdreeksmodel gws WL\Python Scripts\Real-time calibration control scripts')

from RD_converter import Converter


class rt_simulation:
    """
    Important: this real time simulation model only knows the response functions that are at the moment included
    in Pastas (Gamma, Hantush, Exponential, One)
    """
    import sys
    import pandas as pd
    import pastas as ps
    import numpy as np
    import matplotlib.pyplot as plt
    from datetime import datetime
    from datetime import timedelta
    
    def __init__(self, model, name):
        """ The rt simulation requires the following input: a pastas model and 
        a name. Note: the model should be solved on forehand. Important: stressmodel names should correspond with
        the names of the stresses
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
        """ By calling this function the model is updated. Recent forcing data can be added to the 
        current forcing dataset, and the model is recalibrated. 'Selection' holds the length of the calibration period.
        When the historical data should be continuosly be accumulated, and not be cropped, then 'total' should be given 
        as input argument for selection.
        """
        from datetime import timedelta
        # Update forcing data
        forc = forcing.copy()
        initial_data = self.model.get_stress(self.stressmodelnames[0])[:].to_frame()
        initial_data.columns = [self.stressmodelnames[0]]
        for ii in range(len(list(forcing.columns))-1):
            initial_data[self.stressmodelnames[ii+1]] = self.model.get_stress(self.stressmodelnames[ii+1])[:]
        forc = initial_data.append(forc[:])
        
        # Update oseries (Note)
        new_oseries = self.model.oseries.series.to_frame().append(oseries[:])
        if selection != 'total': 
            new_oseries = new_oseries[new_oseries.index[len(new_oseries)-1]-timedelta(days = selection):]
        
        # Create model
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
            
            # The created stressmodels are added to the forecast model
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
            
        
        new_model.solve()
        self.__init__(new_model,self.name)
                 
    def forecast(self,forcing):
        """ By calling this function a forecast is done for a specified
        length. Besides, a forcing dataframe for the predictions should be added. The length of the 
        forecast is inferred from the length of the forcing data. If the 
        past value of the groundwaterlevel is included as input variable, the column in the
        dataframe should be named gw_previous. Optionally add column with observed values for plotting. 
        """
        
        # Check wether stresses in forcing dataframe correspond with the stresses in the model
        if self.stressmodelnames != list(forcing.columns):
            return ("stresses in forcing data do not correspond with stresses in stressmodels")
            
        initial_data = self.model.get_stress(self.stressmodelnames[0])[:].to_frame()
        initial_data.columns = [self.stressmodelnames[0]]
        for ii in range(len(list(forcing.columns))-1):
            initial_data[self.stressmodelnames[ii+1]] = self.model.get_stress(self.stressmodelnames[ii+1])[:]
            
        forc = initial_data.append(forcing[:])
            
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
                
        # Simulate with 
        forecast = forecast_ml.simulate(parameters = self.pars)
        return(forecast[len(forc)-len(forcing):])
