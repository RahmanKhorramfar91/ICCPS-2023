# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:18:26 2021

@author: Rahman Khorramfar
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 18:36:05 2021
filter by high-voltage buses rather than branches with rateA>100
@author: Rahman Khorramfar  
"""
# -*- coding: utf-8 -*-
# from IPython import get_ipython;
# get_ipython().magic('reset -f') # to clear the namespace
# get_ipython().magic('clear');
import numpy as np;
import pandas as pd;
# from geopy.geocoders import Nominatim;
# from geopy.distance import distance;
import os;
# import matplotlib.pyplot as plt;
# import geopandas as gpd;
# import plotly.express as px;
# from geopy.distance import distance;
#%% Functions

def dict_creator(df_bus,df_b2s,df_zone):
    '''
    Returns dictionaries that mapt (buss to sub),
    (sub to zone_id),(bust to zone_id) and (zone_id to state name)   
    '''        
    # zone_id,zone_count = np.unique(df_bus['zone_id'].tolist(),return_counts=True);
    # zone_names = df_zone['zone_name'].unique();
    # zone_id = df_zone['zone_id'].unique();
    # res1 = zone_names.repeat(zone_count)
    # aa = df_b2s;
    # aa['zone_name'] = res1;
    # res1 = zone_id.repeat(zone_count);
    # aa['zone_id'] = res1;
    
    aa = df_bus[['bus_id','sub_id','zone_id']];
    
    s2 = aa.groupby(['sub_id','zone_id']).mean();
    s2 = s2.reset_index();
    bus2zone = dict(zip(aa['bus_id'],aa['zone_id']));
    sub2zone = dict(zip(s2['sub_id'],s2['zone_id']));
    bus2sub= dict(zip(aa['bus_id'],aa['sub_id']));
    return bus2sub, sub2zone,bus2zone;
   
def Region_buses(df_zone,df_bus,df_sub,States):
    
    Reg_zone_ids = df_zone[df_zone['state'].isin(States)];
    Reg_zone_ids=Reg_zone_ids['zone_id'];
    Reg_bus= df_bus[df_bus['zone_id'].isin(Reg_zone_ids)];
    Reg_hv_bus = Reg_bus[Reg_bus['baseKV']>=345];
    
    lat = np.zeros(len(Reg_hv_bus),dtype=float);
    lon = np.zeros(len(Reg_hv_bus),dtype=float);
    for i in range(len(lat)):
        lat[i] = df_sub['lat'][Reg_hv_bus['sub_id'].iloc[i]];
        lon[i] = -df_sub['lon'][Reg_hv_bus['sub_id'].iloc[i]];
        
    Reg_hv_bus['lat'] = lat;
    Reg_hv_bus['lon'] = lon;
    Reg_hv_bus=Reg_hv_bus.reset_index();
    Reg_hv_bus['bus_num'] = np.arange(len(Reg_hv_bus));
    cols = Reg_hv_bus.columns.tolist();
    cols = cols[-1:] + cols[:-1];
    Reg_hv_bus= Reg_hv_bus[cols];
    
    s1 = np.arange(len(Reg_hv_bus));
    s2 = Reg_hv_bus['bus_id'].tolist();
    aa = dict(zip(s2,s1));
    Reg_hv_bus.drop('index', axis=1, inplace=True);
    return Reg_hv_bus,aa;

 
def Region_all_plants(df_zone,df_bus,df_sub,df_plt,
                      df_gen_cost,bus2zone,df_b2s,sub2zone,zone_id2state,States):
    
    NE_zone_ids = df_zone[df_zone['state'].isin(States)];
    NE_zone_ids=NE_zone_ids['zone_id'];
    NE_bus= df_bus[df_bus['zone_id'].isin(NE_zone_ids)];
    NE_hv_bus_id = NE_bus[NE_bus['baseKV']>=0];
    NE_hv_bus_id = NE_hv_bus_id['bus_id'];

    # get plant data, add parameters, and filter for in-service and non-other plants
    df2 = df_plt;  
    df2 = df2[df2['bus_id'].isin(NE_hv_bus_id)];
    df2['startup_cost'] = df_gen_cost['startup'];
    df2['shutdown_cost'] =df_gen_cost['shutdown'];
    df2['fixed_cost'] = df_gen_cost['c0']; # constant component of cost function ($/h)
    df2['variable_cost'] = df_gen_cost['c1']; # linear component
    df2 = df2.rename(columns={'type':'plant_type'});
    # only keep those with status=1 (in-service plants)
    df2 = df2[df2['status']==1];
    df2 = df2.sort_values(by=['bus_id'])
    df3 = df_bus[['bus_id','type','sub_id']];
    df3 = df3.rename(columns = {'type':'bus_type'});
    df2['bus_type'] = np.zeros((len(df2)));
    df4 = df3[df3['bus_id'].isin(df2['bus_id'])]
    bus_ids, bus_id_count = np.unique(df2['bus_id'].tolist(), return_counts=True);
    res1 = df4['bus_type'];
    res1 = res1.repeat(bus_id_count);
    df2['bus_type'] = res1.values;
    
    #remove plants with negligible output
    a = df2;
    a1 = a[a['plant_type']=='ng'];
    a2 = a1['Pmax']<10;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    a = df2;
    a1 = a[a['plant_type']=='wind'];
    a2 = a1['Pmax']<2;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    a = df2;
    a1 = a[a['plant_type']=='nuclear'];
    a2 = a1['Pmax']<10;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    a = df2;
    a1 = a[a['plant_type']=='coal'];
    a2 = a1['Pmax']<10;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    a = df2;
    a1 = a[a['plant_type']=='dfo'];
    a2 = a1['Pmax']<10;
    a3 = a2[a2];
    df2 = a.drop(a3.index);

    a = df2;
    a1 = a[a['plant_type']=='solar'];
    a2 = a1['Pmax']<2;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    a = df2;
    a1 = a[a['plant_type']=='hydro'];
    a2 = a1['Pmax']<2;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    a = df2;
    a1 = a[a['plant_type']=='wind_offshore'];
    a2 = a1['Pmax']<2;
    a3 = a2[a2];
    df2 = a.drop(a3.index);
    
    s1 = (0.5*df2['GenIOB']*df2['Pmax']**2
          ) + (2/3)*df2['GenIOC']*df2['Pmax']**3-(
              (0.5*df2['GenIOB']*df2['Pmin']**2
               ) + (2/3)*df2['GenIOC']*df2['Pmin']**3);
    s2 = 0.5*df2['Pmax']**2-0.5*df2['Pmin']**2;
    df2['heat_rate'] = s1/s2;
    
    s1 = (0.5*df2['GenIOB']*df2['Pmax']**2
          ) + (2/3)*df2['GenIOC']*df2['Pmax']**3-(
              (0.5*df2['GenIOB']*(df2['Pmin']-0.1)**2
               ) + (2/3)*df2['GenIOC']*(df2['Pmin']-0.1)**3);
    s2 = 0.5*df2['Pmax']**2-0.5*(df2['Pmin']-0.1)**2;
    df2['heat_rate_near_Pmax'] = s1/s2;
    
    s1 = (df2['GenIOB']*df2['Pmax']+df2['GenIOC']*df2['Pmax']**2
          )-(df2['GenIOB']*df2['Pmin']-df2['GenIOC']*df2['Pmin']**2);
    s2 = df2['Pmax']-df2['Pmin'];
    df2['heat_rate_linear_slope'] = s1/s2;
    
    
    #df2 = df2[df2['plant_type']=='ng'];
    #df2.to_csv('Texas-all-plants.csv',encoding='utf-8',index = True);
    # aggregate based on bus_id and plant_type
    agg_cols = {'Qmax':['mean'],'Qmin':['mean'],
              'Vg':['mean'],'mBase':['mean'],
              'Pmax':['sum'],'Pmin':['sum'],
              'Pc1':['mean'],'Pc2':['mean'],
              'Qc1min':['mean'],'Qc1max':['mean'],
              'Qc2min':['mean'],'Qc2max':['mean'],
              'ramp_agc':['sum'],'ramp_10':['sum'],
              'ramp_30':['sum'],'ramp_q':'sum',
              'apf':['mean'],'mu_Pmax':'mean',
              'mu_Pmin':['mean'],'mu_Qmax':'mean',
              'mu_Qmin':['mean'],'GenFuelCost':['sum'],
              'GenIOB':['mean'],'GenIOC':['mean'],'GenIOD':['mean'],
              'startup_cost':['sum'],'shutdown_cost':['sum'],
              'fixed_cost':['sum'],'variable_cost':['sum']};
      
    res2 = df2.groupby(['bus_id','plant_type']).agg(agg_cols);
    
    res2.columns =  [ 'Qmax','Qmin','Vg','mBase',
        'Pmax','Pmin','Pc1','Pc2','Qc1min',
        'Qc1max','Qc2min','Qc2max',
        'ramp_agc','ramp_10','ramp_30','ramp_q',
        'apf','mu_Pmax','mu_Pmin','mu_Qmax',
        'mu_Qmin','GenFuelCost','GenIOB','GenIOC','GenIOD','startup_cost','shutdown_cost',
        'fixed_cost','variable_cost'];
    res2 = res2.reset_index();
    
    
    # remove type=other plants
    aa = res2['plant_type']=='other'
    aa = np.array( ~np.array(aa));
    df2=res2[aa];
    
      
    # <-
    df2 = df2.reset_index();
    df2['state'] = np.array(len(df2));
    df2['sub_id'] = np.zeros(len(df2),dtype=int);
    df2['lat'] = np.zeros(len(df2),dtype=float);
    df2['lon'] = np.zeros(len(df2),dtype=float);
    df2['zone_id'] =np.zeros(len(df2),dtype=int); 
    #df2['bus_num'] = np.zeros(len(df2),dtype=int);
    for i in range(len(df2)):   
        df2['zone_id'][i] = bus2zone[df2['bus_id'][i]];
        a2 = np.where(df_b2s['bus_id']==df2['bus_id'][i])[0];
        df2['sub_id'][i] = df_b2s['sub_id'][a2[0]];
        #df2['sub_id'][i] = (df_b2s['sub_id'].iloc[df2['bus_id'][i]-1]);
        df2['lat'][i] = df_sub['lat'].iloc[df2['sub_id'][i]];
        df2['lon'][i] = -df_sub['lon'].iloc[df2['sub_id'][i]];
        zone_id = sub2zone[df2['sub_id'][i]];
        df2['state'][i]= zone_id2state[zone_id];   
        #df2['bus_num'][i] = busID2num[df2['bus_id'].iloc[i]];
        
    return df2;   

def Region_plants_in_bus(Reg_plants,Reg_hv_bus):
    pltBus = np.zeros(len(Reg_plants),dtype=int);

    for i in range(len(Reg_plants)):
        plt_loc = np.array([Reg_plants['lat'].iloc[i],Reg_plants['lon'].iloc[i]]);
        dist = np.zeros(len(Reg_hv_bus));
        for j in range(len(Reg_hv_bus)):
            dist[j] = np.linalg.norm(plt_loc-np.array([Reg_hv_bus['lat'].iloc[j],Reg_hv_bus['lon'].iloc[j]]));
    
        arg_sort = np.argsort(dist);
        pltBus[i] = int(arg_sort[0]);
    
    
    Reg_plants['bus_num'] = pltBus;
    Reg_plants['count'] =np.ones(len(pltBus));
    
    
    
    # apply groupby to bus_num and type
    
    s1 = {'Pmax':'sum','Pmin':'sum','GenFuelCost':'mean',
          'GenIOB':'mean','GenIOC':'mean','ramp_30':'sum','count':'sum'};
    df_nep = Reg_plants.groupby(['bus_num','plant_type']).agg(s1);
    df_nep2 = Reg_plants.groupby(['bus_num','plant_type']).size().reset_index(name='counts');
    df_nep = df_nep.reset_index();
    
    return df_nep.reset_index();

   
def add_FIPS(df_nec,NE_hv_bus):
    dict1 = {'Connecticut':'CT','Maine':'ME','Massachusetts':'MA',
             'New Hampshire':'NH','Rhode Island':'RI','Vermont':'VT',
            'New York':'NY', 'New Brunswick':'NB','New Brunswick / Nouveau-Brunswick':'NB'};
    cts = np.zeros(len(df_nec));
    c2f = {};
    for i in range(len(cts)):
        s1 = df_nec['State'].iloc[i]+'-'+df_nec['County'].iloc[i];
        c2f[s1] = df_nec['FIPS'].iloc[i];
    geolocator = Nominatim(user_agent="geoapiExercises");
    fips = np.zeros(len(NE_hv_bus));
    cty = np.empty(len(NE_hv_bus),dtype=object);
    # buses from 4830 to 4843 are located in the ocean so are beware
    #NE_hv_bus.drop(index=range(4830,4843), inplace=True);
    for i in range(len(NE_hv_bus)):#
        if i==4830: # this is a location in ocean near Maine
            fips[i] = c2f['ME-Washington'];
            cty[i] = 'ME-Washington';
            continue;
        elif i==4831 :
            fips[i] = c2f['ME-Hancock'];
            cty[i] = 'ME-Hancock';
            continue;
        elif i==4832:
            fips[i] = c2f['ME-Lincoln'];
            cty[i] = 'ME-Lincoln';
            continue;
        elif i==4833:
            fips[i] = c2f['ME-York'];
            cty[i] = 'ME-York';
            continue;
        elif i==4834:
            fips[i] = c2f['MA-Barnstable'];
            cty[i] = 'MA-Barnstable';
            continue;
        elif i in range(4835,4840):
            fips[i] = c2f['MA-Dukes'];
            cty[i] = 'MA-Dukes';
            continue;
        elif i in range(4840,4844):
            fips[i] = c2f['RI-Washington'];
            cty[i] = 'RI-Washington';
            continue;    
        
        location = geolocator.reverse(str(NE_hv_bus['lat'].iloc[i])+","+str(-NE_hv_bus['lon'].iloc[i]))
    
        s1 = location.raw['address']['county']
        s1 = s1.replace(" County","");
        s2 = dict1[location.raw['address']['state']]; # get the symbol
        s2 = s2+'-'+s1;
        #print(i,"\t",s2)
        if s2=='NB-Victoria' or s2=='NB-Madawaska':
            fips[i] = c2f['ME-Aroostook'];
            cty[i] ='ME-Aroostook';
        elif s2=="RI-South":
            fips[i] = c2f['RI-Washington'];
        elif s2=='NY-Westchester' or s2=='NY-Suffolk' or s2=='NY-Putnam' or s2=='NY-Nassau' or s2=='NY-Dutchess':
            fips[i] = c2f['CT-Fairfield'];
            cty[i] = 'CT-Fairfield';
        elif s2=='NY-Washington':
            fips[i] = c2f['VT-Rutland'];
            cty[i] = 'VT-Rutland';
        else:
            fips[i] = c2f[s2];
            cty[i] = s2;
            
    NE_hv_bus['fips'] = fips;
    NE_hv_bus['st-cty'] = cty;
    
    return NE_hv_bus;
    


def Plants_param_and_count(Reg_all_plants):
    
    Reg_all_plants['count'] = np.ones(len(Reg_all_plants));
    pt =np.array(['ng','dfo','solar','wind','wind_offshore','hydro','coal','nuclear']);    
    s1 = {'Pmax':'sum','Pmin':'sum','GenFuelCost':'mean',
          'GenIOB':'mean','GenIOC':'mean','GenIOD':'mean','ramp_30':'sum','count':'sum'};   
    Reg_plants = Reg_all_plants.groupby(['zone_id','plant_type']).agg(s1);
     #Reg_plants['count'] =np.ones(len(Reg_plants));
    Reg_plants = Reg_plants.reset_index();

    Maxp = np.zeros(len(pt));
    max_total_cap = np.zeros(len(pt));
    min_total_cap = np.zeros(len(pt));
    AveP = np.zeros(len(pt));
    rampU = np.zeros(len(pt));
    heat_rate = np.zeros(len(pt));
    fuel_cost = np.zeros(len(pt));
    Minp = np.zeros(len(pt));
       
    for i in range(len(pt)):    
        a1 = Reg_all_plants[Reg_all_plants['plant_type']==pt[i]];
        b1 = a1[a1['ramp_30']!=0];
        rampU[i] = 2*np.mean(b1['ramp_30']/b1['Pmax']); #per hour
        Maxp[i] = a1['Pmax'].mean();
        max_total_cap[i] =a1['Pmax'].sum()*8760;
        min_total_cap[i] =a1['Pmin'].sum()*8760;
        Minp[i] = a1['Pmin'].mean();
        AveP[i] = a1['Pmax'].sum()/len(a1);
        heat_rate[i] = a1['GenIOB'].mean()+np.mean(((a1['Pmax'])**2*(a1['GenIOC'])))/Maxp[i];
        fuel_cost[i] = a1['GenFuelCost'].mean();
    
    #df_nep['count'] = np.array(df_nep2['counts']);
    Reg_plants['adjusted_count'] = np.zeros(len(Reg_plants));
    for i in range(len(Reg_plants)):
        # identify the plant type
        ind1 = np.where(pt==Reg_plants['plant_type'].iloc[i]);
        ind1 = ind1[0];
        Reg_plants['adjusted_count'][i] =max(1,np.round(Reg_plants['Pmax'].iloc[i]/AveP[ind1]));
    
    return Reg_plants,Maxp,Minp,rampU,heat_rate,fuel_cost,max_total_cap,min_total_cap;

def Region_branch_buses(NE_bus,df_branch,df_b2s,busID2num):   
    
    # #zone_id2state,bus2sub, sub2zone,bus2zone=Funcs_hvb.dict_creator(df_bus,df_b2s,df_zone);             
    # NorthEast_states = ['Connecticut', 'Maine', 'Massachusetts',
    #     'New Hampshire', 'Rhode Island', 'Vermont','New York'];    
    # NorthEast_states = ['Connecticut', 'Maine', 'Massachusetts',
    #                     'New Hampshire', 'Rhode Island', 'Vermont'];
    # #NorthEast_states = ['Vermont'];
    # NE_zone_ids = df_zone[df_zone['state'].isin(NorthEast_states)];
    # NE_zone_ids=NE_zone_ids['zone_id'];
    # NE_bus= df_bus[df_bus['zone_id'].isin(NE_zone_ids)];
    # NE_hv_bus_id = NE_bus[NE_bus['baseKV']>=115];
    NE_hv_bus_id = NE_bus['bus_id'];
    
    NE_highVol_br = df_branch[df_branch['from_bus_id'].isin(NE_hv_bus_id)];
    NE_highVol_br = NE_highVol_br[NE_highVol_br['to_bus_id'].isin(NE_hv_bus_id)];
    
    NE_highVol_br['suscept'] = 1/NE_highVol_br['x'];    
    oper = {'r':'sum','suscept':'sum','rateA':'sum','x':'min','b':'sum',
                'rateB':'mean','rateC':'mean','ratio':'mean',
                'angle':'mean','angmin':'mean','angmax':'mean'};
    
    df4 = NE_highVol_br.groupby(['from_bus_id','to_bus_id'],as_index = False).agg(oper);
    df4['from_bus_num'] = np.zeros(len(df4),dtype=int);
    df4['to_bus_num'] = np.zeros(len(df4),dtype=int);
    df4['length'] = np.zeros(len(df4),dtype=float);
    df4['from_kv'] = np.zeros(len(df4),dtype=float);
    df4['to_kv'] = np.zeros(len(df4),dtype=float);
    for i in range(len(df4)):
        df4['from_bus_num'][i] = busID2num[df4['from_bus_id'].iloc[i]];
        df4['from_kv'][i] = NE_bus['baseKV'].iloc[df4['from_bus_num'][i]];
        
        df4['to_bus_num'][i] = busID2num[df4['to_bus_id'].iloc[i]];
        df4['to_kv'][i] = NE_bus['baseKV'].iloc[df4['to_bus_num'][i]];
        
        loc1 = [NE_bus['lat'].iloc[df4['from_bus_num'].iloc[i]],NE_bus['lon'].iloc[df4['from_bus_num'].iloc[i]]];
        loc2 = [NE_bus['lat'].iloc[df4['to_bus_num'].iloc[i]],NE_bus['lon'].iloc[df4['to_bus_num'].iloc[i]]];
        df4['length'][i] = np.linalg.norm(np.array(loc1)-np.array(loc2));
        
    return df4;


def Branch_data(df_ne_bus, df_ne_br,df_ne_plt,path):
        
    plant2sym = {1:'ng',2:'dfo',3:'solar',4:'wind',5:'wind_offshore', 6:'hydro',7:'coal',8:'nuclear'};
    sym2plant = {'ng':1,'dfo':2,'solar':3,'wind':4,'wind_offshore':5,'hydro':6,'coal':7,'nuclear':8};
    
    dfb = df_ne_br;
    # get the existing branches
    AdjNodes = []; is_exist = [];
    Br = [];LnDist = [];Sus = [];maxF = [];
    for i in range(len(df_ne_bus)):
        AdjNodes.append([]); is_exist.append([]);
        Br.append([]);LnDist.append([]);Sus.append([]);maxF.append([]);
    
    for i in range(len(df_ne_br)):
        fb = df_ne_br['from_bus_num'].iloc[i];
        tb = df_ne_br['to_bus_num'].iloc[i];
        s1 = min(fb,tb);
        tb = max(fb,tb);
        fb = s1;
        
        Br[fb].append(tb); 
        LnDist[fb].append(df_ne_br['length'].iloc[i]);
        is_exist[fb].append(1);
        Sus[fb].append(df_ne_br['suscept'].iloc[i]);
        maxF[fb].append(df_ne_br['rateA'].iloc[i]);  
        AdjNodes[fb].append(tb);
        AdjNodes[tb].append(fb);        
    
    # get the max dist for each bus
    dfb2 = dfb[['from_bus_num','length']];
    dfb2.rename(columns = {'from_bus_num':'bus'}, inplace= True)
    dfb3 = dfb[['to_bus_num','length']];
    dfb3.rename(columns = {'to_bus_num':'bus'}, inplace= True)
    dfb4 = [dfb2,dfb3];
    dfb5 = pd.concat(dfb4)
    dfb6 = dfb5.groupby(['bus']).agg(mean=('length','mean'),Max=('length','max'));
    dfb6 = dfb6.reset_index()    
    
    # remove the following three lines
    # bus number 3814 is not connected. Add it manually
    # row = pd.DataFrame({'bus':3814,'mean':0,'Max':70},index=[3814]);     
    # dfb6 = pd.concat([dfb6.iloc[:3814],row,dfb6.iloc[3814:]]).reset_index(drop=True);
    # del dfb2,dfb3,dfb4,dfb5;
    
    #% apply a regression model for susceptance and max-flow
    dfr1 = dfb[dfb['suscept']!=0];
    from sklearn import linear_model;
    X = dfr1[['length','from_kv','to_kv']];
    Ys = dfr1['suscept'];
    Yf = dfr1['rateA'];
    regS = linear_model.LinearRegression();
    regS.fit(X,Ys);
    regF = linear_model.LinearRegression();
    regF.fit(X,Yf);
    
    #% create Jn and Dist for Jn
    for i in range(len(df_ne_bus)):  
    
        bn = df_ne_bus['bus_num'].iloc[i];
        lat = df_ne_bus['lat'][bn];lon = df_ne_bus['lon'][bn];
        
        if dfb6['Max'].iloc[i]==0:
            mxm = 100.0;
        s1 = (abs(df_ne_bus['lat']-lat)+abs(df_ne_bus['lon']-lon))<(mxm/30)+0.6;
        # if i==101:
        #         s1 = (abs(df_ne_bus['lat']-lat)+abs(df_ne_bus['lon']-lon))<(dfb6['Max'].iloc[i]/30)+0.5;
    
        s2 = df_ne_bus[s1];s2=s2.reset_index();    
        dist = np.zeros(len(s2),dtype=float);   
        for j in range(len(s2)):
            dist[j] = np.linalg.norm(np.array([lat,lon])-np.array([s2['lat'].iloc[j],s2['lon'].iloc[j]]));
        
        # find the nearest bus with dist>0
        a1 = min(dist[dist>0]);
        ind2 = np.where(dist==a1);
        ind2 = ind2[0][0];
        a2 = s2['bus_num'].iloc[ind2];
        if a2 not in AdjNodes[i]:
            m1 = min(i,a2);
            m2 = max(i,a2);
            AdjNodes[m1].append(m2);
            AdjNodes[m2].append(m1);
            Br[m1].append(m2);
            is_exist[m1].append(0); # check if it is already existed! luckily there is no candidate branch who was existed
            LnDist[m1].append(np.round(a1,1));
            Sus[m1].append(np.round(regS.predict([[a1,df_ne_bus['baseKV'][i],s2['baseKV'][ind2]]])[0],1));
            maxF[m1].append(np.round(regF.predict([[a1,df_ne_bus['baseKV'][i],s2['baseKV'][ind2]]])[0],1));
        
    BRs = list([]);
    for i in range(len(Br)):
        for j in range(len(Br[i])):
            BRs.append([i,Br[i][j], is_exist[i][j], maxF[i][j],Sus[i][j], LnDist[i][j]]);
    df_BRs = pd.DataFrame(BRs);
    nBus = len(AdjNodes);
    for i in range(nBus):
        b1 = AdjNodes[i];
        b1 = np.unique(b1);
        AdjNodes[i] = b1;
    df_adj = pd.DataFrame(AdjNodes);
    #% export data as a csv files
    # df_adj = pd.DataFrame(AdjNodes);
    # df_Br = pd.DataFrame(Br);
    # df_LnDist = pd.DataFrame(LnDist);
    # df_Sus = pd.DataFrame(Sus);
    # df_mf = pd.DataFrame(maxF);
    # df_ie = pd.DataFrame(is_exist);
    # cols=[''];
    # path2=os.path.join(path,'bus_adj_Nodes.csv');
    # df_adj.to_csv(path2,encoding='utf-8',index = True);
    
    # path2=os.path.join(path,'b2b_br.csv');
    # df_Br.to_csv(path2,encoding='utf-8',index = True);
    
    # path2=os.path.join(path,'b2b_br_dist.csv');
    # df_LnDist.to_csv(path2,encoding='utf-8',index = True);
    
    # path2=os.path.join(path,'b2b_br_Suscept.csv');
    # df_Sus.to_csv(path2,encoding='utf-8',index = True);
    
    # path2=os.path.join(path,'b2b_br_maxFlow.csv');
    # df_mf.to_csv(path2,encoding='utf-8',index = True);
    
    # path2=os.path.join(path,'b2b_br_is_existing.csv');
    # df_ie.to_csv(path2,encoding='utf-8',index = True);
    return df_BRs,df_adj;

def get_zonal_solar_profile(df_zone, df_bus, df_plt,df_solar,path,States):
    
    Reg_zone_ids = df_zone[df_zone['state'].isin(States)];
    Reg_zone_ids=(Reg_zone_ids['zone_id']);
    zonal_profile = np.zeros((len(df_solar),len(Reg_zone_ids)));
    
    for i in range(len(Reg_zone_ids)):
        Reg_bus= df_bus[df_bus['zone_id'].isin(pd.Series(Reg_zone_ids.iloc[i]))];    
        df_reg_plt = df_plt[df_plt['bus_id'].isin(Reg_bus['bus_id'])];
        reg_plt_id = np.array(df_reg_plt['plant_id']);
        aa= np.array(df_solar.columns[1:].tolist());
        aa = aa.astype(int);
        a2 =list(set(aa) & set(reg_plt_id));
        dfs_Tr = df_solar.transpose();
        dfs_Tr = dfs_Tr.iloc[1:,:];
        dfs_Tr.insert(0,'plant_id',aa);
        a2 = dfs_Tr[dfs_Tr['plant_id'].isin(reg_plt_id)];
        if len(a2)==0: continue; 
        for p in range(len(a2)):
            pl1 = a2['plant_id'][p];
            pmax = df_plt['Pmax'][pl1];
            if pmax<0.5:
                pmax = np.max(a2.iloc[p,1:])+0.01;
                #pmax = np.max(pmax,1);
            a2.iloc[p,1:] = a2.iloc[p,1:]/pmax;                
        a3 = a2.transpose(); 
        a3 = a3.iloc[1:,:];    
        zonal_profile[:,i] = a3.mean(axis=1);
    
    solar_profile = pd.DataFrame(zonal_profile);
    path2=os.path.join(path,'profile_solar_hourly.csv');
    solar_profile.to_csv(path2,encoding='utf-8',index = False);
    
    
def get_zonal_hydro_profile(df_zone, df_bus, df_plt,df_hydro,path,States):
    
    Reg_zone_ids = df_zone[df_zone['state'].isin(States)];
    Reg_zone_ids=(Reg_zone_ids['zone_id']);
    zonal_profile = np.zeros((len(df_hydro),len(Reg_zone_ids)));
    
    for i in range(len(Reg_zone_ids)):
        Reg_bus= df_bus[df_bus['zone_id'].isin(pd.Series(Reg_zone_ids.iloc[i]))];    
        df_reg_plt = df_plt[df_plt['bus_id'].isin(Reg_bus['bus_id'])];
        reg_plt_id = np.array(df_reg_plt['plant_id']);
        aa= np.array(df_hydro.columns[1:].tolist());
        aa = aa.astype(int);
        a2 =list(set(aa) & set(reg_plt_id));
        dfs_Tr = df_hydro.transpose();
        dfs_Tr = dfs_Tr.iloc[1:,:];
        dfs_Tr.insert(0,'plant_id',aa);
        a2 = dfs_Tr[dfs_Tr['plant_id'].isin(reg_plt_id)];
        if len(a2)==0: continue; 
        for p in range(len(a2)):
            pl1 = a2['plant_id'][p];
            pmax = df_plt['Pmax'][pl1];
            if pmax<0.5:
                pmax = np.max(a2.iloc[p,1:])+0.01;
                #pmax = np.max(pmax,1);
            a2.iloc[p,1:] = a2.iloc[p,1:]/pmax;                
        a3 = a2.transpose(); 
        a3 = a3.iloc[1:,:];    
        zonal_profile[:,i] = a3.mean(axis=1);
    
    solar_profile = pd.DataFrame(zonal_profile);
    path2=os.path.join(path,'profile_hydro_hourly.csv');
    solar_profile.to_csv(path2,encoding='utf-8',index = False);
   

def get_zonal_wind_profile(df_zone, df_bus, df_plt,df_wind,path,States):
    
    Reg_zone_ids = df_zone[df_zone['state'].isin(States)];
    Reg_zone_ids=(Reg_zone_ids['zone_id']);
    zonal_profile = np.zeros((len(df_wind),len(Reg_zone_ids)));
    
    for i in range(len(Reg_zone_ids)):
        Reg_bus= df_bus[df_bus['zone_id'].isin(pd.Series(Reg_zone_ids.iloc[i]))];    
        df_reg_plt = df_plt[df_plt['bus_id'].isin(Reg_bus['bus_id'])];
        reg_plt_id = np.array(df_reg_plt['plant_id']);
        aa= np.array(df_wind.columns[1:].tolist());
        aa = aa.astype(int);
        a2 =list(set(aa) & set(reg_plt_id));
        dfs_Tr = df_wind.transpose();
        dfs_Tr = dfs_Tr.iloc[1:,:];
        dfs_Tr.insert(0,'plant_id',aa);
        a2 = dfs_Tr[dfs_Tr['plant_id'].isin(reg_plt_id)];
        if len(a2)==0: continue; 
        for p in range(len(a2)):
            pl1 = a2['plant_id'][p];
            pmax = df_plt['Pmax'][pl1];
            if pmax<0.5:
                pmax = np.max(a2.iloc[p,1:])+0.01;
                #pmax = np.max(pmax,1);
            a2.iloc[p,1:] = a2.iloc[p,1:]/pmax;                
        a3 = a2.transpose(); 
        a3 = a3.iloc[1:,:];    
        zonal_profile[:,i] = a3.mean(axis=1);
    
    solar_profile = pd.DataFrame(zonal_profile);
    path2=os.path.join(path,'profile_wind_hourly.csv');
    solar_profile.to_csv(path2,encoding='utf-8',index = False);
    
   
def get_solar_profile(df_zone, df_bus, df_plt,df_solar,path,States):
        
    NE_zone_ids = df_zone[df_zone['state'].isin(States)];
    NE_zone_ids=NE_zone_ids['zone_id'];
    NE_bus= df_bus[df_bus['zone_id'].isin(NE_zone_ids)];
    
    df_ne_plt = df_plt[df_plt['bus_id'].isin(NE_bus['bus_id'])];
    ne_plt_id = np.array(df_ne_plt['plant_id']);
    
    aa= np.array(df_solar.columns[1:].tolist());
    aa = aa.astype(int)
    a2 =list( set(aa) & set(ne_plt_id));
    a2 = np.sort(a2);
    
    # find the index of plants
    a3 = np.zeros(len(a2),dtype=int);
    for i in range(len(a2)):
        a3[i] = np.where(aa==a2[i])[0]+1;
        
    aa = df_solar.iloc[:,a3];
    #aa = aa.drop('UTC',1);
    
    # take average to get a profile
    solar_profile = aa.mean(axis=1);
    
    path2=os.path.join(path,'profile_solar_hourly.csv');
    solar_profile.to_csv(path2,encoding='utf-8',index = False);
    
def get_hydro_profile(df_zone, df_bus, df_plt,df_hydro,path,States):
        
    NE_zone_ids = df_zone[df_zone['state'].isin(States)];
    NE_zone_ids=NE_zone_ids['zone_id'];
    NE_bus= df_bus[df_bus['zone_id'].isin(NE_zone_ids)];
    
    df_ne_plt = df_plt[df_plt['bus_id'].isin(NE_bus['bus_id'])];
    ne_plt_id = np.array(df_ne_plt['plant_id']);
    
    aa= np.array(df_hydro.columns[1:].tolist());
    aa = aa.astype(int)
    a2 =list( set(aa) & set(ne_plt_id));
    a2 = np.sort(a2);
    
    # find the index of NE plants
    a3 = np.zeros(len(a2),dtype=int);
    for i in range(len(a2)):
        a3[i] = np.where(aa==a2[i])[0];
        
    
    aa = df_hydro.iloc[:,a3];
    aa = aa.drop('UTC',1);
    
    # take average to get a profile
    hydro_profile = aa.mean(axis=1);
    
    path2=os.path.join(path,'profile_hydro_hourly.csv');
    hydro_profile.to_csv(path2,encoding='utf-8',index = False);
    

      
def get_wind_profile(df_zone, df_bus, df_plt,df_wind,path,States):
        
    NE_zone_ids = df_zone[df_zone['state'].isin(States)];
    NE_zone_ids=NE_zone_ids['zone_id'];
    NE_bus= df_bus[df_bus['zone_id'].isin(NE_zone_ids)];
    
    df_ne_plt = df_plt[df_plt['bus_id'].isin(NE_bus['bus_id'])];
    ne_plt_id = np.array(df_ne_plt['plant_id']);
    
    aa= np.array(df_wind.columns[1:].tolist());
    aa = aa.astype(int)
    a2 =list( set(aa) & set(ne_plt_id));
    a2 = np.sort(a2);
    
    # find the index of NE plants
    a3 = np.zeros(len(a2),dtype=int);
    for i in range(len(a2)):
        a3[i] = np.where(aa==a2[i])[0];
        
    
    aa = df_wind.iloc[:,a3];
    aa = aa.drop('UTC',1);
    
    # take average to get a profile
    hydro_profile = aa.mean(axis=1);
    
    path2=os.path.join(path,'profile_wind_hourly.csv');
    hydro_profile.to_csv(path2,encoding='utf-8',index = False);
    
    
    

