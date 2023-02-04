# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 11:43:06 2022
disaggrega the state-level CFs to bus level.
Assumption is that 
all buses in the same state have the same CF
@author: Rahman Khorramfar
"""

from IPython import get_ipython;
get_ipython().magic('reset -f') # to clear the namespace
get_ipython().magic('clear');
import numpy as np;
import pandas as pd;
import random;
from geopy.distance import distance;
import os;

dfs = pd.read_csv('profile_solar_hourly.csv');
dfw = pd.read_csv('profile_wind_hourly.csv');

dfb = pd.read_csv('Region_buses.csv');

# Clus = dfc['Cluster'].unique();

sh = np.zeros((len(dfs),len(dfb)));
wh = np.zeros((len(dfs),len(dfb)));

for i in range(len(dfb)):
    s1 = dfb['zone_id'].iloc[i]-1;
    sh[:,i] = dfs.iloc[:,s1];    
    wh[:,i] = dfw.iloc[:,s1];

dfs = pd.DataFrame(sh);
dfs.to_csv('solar_CF-188_node.csv');

dfw = pd.DataFrame(wh);
dfw.to_csv('wind_CF-188_node.csv');

























