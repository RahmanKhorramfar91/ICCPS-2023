# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 12:50:15 2022
re-configure the data for the new graph segementation
@author: Rahman Khorramfar
"""

#%% call packages
# from IPython import get_ipython;
# get_ipython().magic('reset -f') # to clear the namespace
# get_ipython().magic('clear');
import numpy as np;
import pandas as pd;
import random;
# from geopy.distance import distance;
import os;
import elec_funcs;

def Cluster_reconfig(nC,method):
#% Read base grid files
    if nC==88:
        name =os.getcwd()+ '/power-network-clusters' +'/'+ str(nC)+'-nodes.csv';
    else:
        name =os.getcwd()+ '/power-network-clusters' +'/'+ str(nC)+'-nodes-' + method+ '.csv';
    dfc = pd.read_csv(name);
    node2zone = dict(zip(dfc['Node'],dfc['Cluster']));
    path = os.getcwd()+'/all_cvs_files-full-network';
    
    dfb = pd.read_csv(path+'/Region_buses.csv');
    nB = len(dfc['Cluster'].unique());
    cls2 = np.array(dfc['Cluster'].unique());
    Nodes = np.zeros((nB,2));
    count = 0;
    for i in cls2:
        s1 = dfc[dfc['Cluster']==i];
        s1 = np.array(s1['Node']);
        Nodes[count,0] = np.mean(dfb['lat'].iloc[s1]);
        Nodes[count,1] = np.mean(dfb['lon'].iloc[s1]);
        
        count += 1;
    
    #% branch data
    
    dfb = pd.read_csv(path+'/Branches.csv');
    dfb2 = pd.DataFrame(columns=dfb.columns);
    count = 0;
    for i in range(len(dfb)):
        dfb.loc[i,'from bus'] = node2zone[dfb['from bus'].iloc[i]];
        dfb.loc[i,'to bus'] = node2zone[dfb['to bus'].iloc[i]];
        f1 =dfb['from bus'].iloc[i];
        t1 =dfb['to bus'].iloc[i];
        if t1!=f1:
            dfb2.loc[count] = np.zeros(len(dfb2.columns));
            dfb2.iloc[count,:] = dfb.iloc[i,:];
            count += 1;
    
    del f1,t1,count,dfb,i;
    
    #% plant data
    dfp = pd.read_csv(path+'/Region_plant.csv');
    dfp['zone_id'] = np.zeros(len(dfp));
    plt_cap = {'ng':173,'solar':6.3,'wind':42,'hydro':23,'nuclear':933};
    for i in range(len(dfp)):
        dfp.loc[i,'zone_id']= int(node2zone[dfp['bus_num'].iloc[i]]+1);
    
    s1 = {'Pmax':'sum','Pmin':'sum','GenFuelCost':'mean',
              'GenIOB':'mean','GenIOC':'mean','ramp_30':'sum','count':'sum'};   
    dfp1 = dfp.groupby(['zone_id','plant_type']).agg(s1);
    
    dfp1['adjusted_count'] = np.zeros(len(dfp1),dtype=int);
    dfp1 = dfp1.reset_index();
    dfp2 = pd.DataFrame(columns = dfp1.columns);
    count = 0;
    for i in range(len(dfp1)):   
        if dfp1['plant_type'].iloc[i]=='dfo' or dfp1['plant_type'].iloc[i]=='coal' or dfp1['plant_type'].iloc[i]=='wind_offshore':
            continue;
        
        dfp1.loc[i,'adjusted_count'] = max(1, np.round(dfp1['Pmax'].iloc[i]/plt_cap[dfp1['plant_type'].iloc[i]]));
        dfp2.loc[count] = np.zeros(len(dfp1.columns));
        dfp2.iloc[count,:] = dfp1.iloc[i,:];
        count += 1;
    
    del dfp, dfp1,i,count,s1;
    
    #% CFs and load 
    dfs = pd.read_csv(path+'/solar-CF-188-nodes-synthetic.csv');
    dfw = pd.read_csv(path+'/solar-CF-188-nodes-synthetic.csv');
    dfl = pd.read_csv(path+'/bus_load_RM_2050.csv');
    nB = len(dfc['Cluster'].unique())
    CFs = np.zeros((len(dfs),nB));
    CFw = np.zeros((len(dfw),nB));
    Load = np.zeros((len(dfl),nB));
    for i in range(nB):
        s1 = np.array(dfc['Cluster']);
        s1 = np.where(s1==i);
        s1 = s1[0]+1;
        CFs[:,i] = np.mean(np.array(dfs.iloc[:,s1]),1);  
        CFw[:,i] = np.mean(np.array(dfw.iloc[:,s1]),1);
        Load[:,i] = np.sum(np.array(dfl.iloc[:,s1]),1);
    
    dfs = pd.DataFrame(CFs);
    dfw = pd.DataFrame(CFw);
    dfl = pd.DataFrame(Load);
    
    
    #% get adjE of NG nodes
    path = os.getcwd()+'/Gas_System_Data';
    dfg=pd.read_csv(path+'/new_england_ng_nodes.csv');
    path = os.getcwd()+'/Power_System_Data';
    df_b = pd.read_csv(path+'/Region_Nodes.csv');
    adjE = [];
    adjE_dist=[];
    for i in range(len(dfg)):
        adjE.append([]);
        adjE_dist.append([]);
    
    
    nB = len(df_b);
    near_ng = 3;
    # for every bus, find near_ng nearest NG node
    for b in range(nB):
        lat_lon = np.array([df_b['lat'].iloc[b],df_b['lon'].iloc[b]]);    
        
        dist = np.zeros(len(dfg));
        for i in range(len(dfg)):
            #dist[i] = distance(lat_lon,[dfg['Lat'].iloc[i],dfg['Lon'].iloc[i]]).miles;
            dist[i] = np.linalg.norm(lat_lon-np.array([dfg['Lat'].iloc[i],dfg['Lon'].iloc[i]]));
            dist[i] = np.round(dist[i],1);
        args = np.argsort(dist);
        s1 = np.sort(dist);
        for j in range(near_ng):
            adjE[args[j]].append(b);
            adjE_dist[args[j]].append(s1[j]);
        
    df_adjE = pd.DataFrame(adjE);
    df_adjE_dist = pd.DataFrame(adjE_dist);
    
    
        
    #% Save files in the Power_System_Data folder
    path = os.getcwd()+'/Power_System_Data';
    dfb2.to_csv(path+'/Branches.csv');
    
    dfp2.to_csv(path+'/Region_plants.csv')
    
    dfs.to_csv(path+'/profile_solar_hourly.csv');
    dfw.to_csv(path+'/profile_wind_hourly.csv');
    dfl.to_csv(path+'/zonal_load-RM.csv');
    
    
    dfn = pd.DataFrame(Nodes,columns=['lat','lon']);
    dfn.to_csv(path+'/Region_Nodes.csv',index=False);
    
    path = os.getcwd()+'/Gas_System_Data';
    df_adjE.to_csv(path+'/ng_adjE.csv',index=False);
        
