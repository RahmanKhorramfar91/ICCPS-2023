# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:29:25 2022

@author: Rahman Khorramfar
Allah sene tevekkul
"""
# from IPython import get_ipython;
# get_ipython().magic('reset -f') # to clear the namespace
# get_ipython().magic('clear');
import numpy as np;
from Setting import Setting,EV,GV;
import pandas as pd
import os;
# from geopy.distance import distance;
def time_weights(nD,fName):   
    #df = pd.read_csv(os.getcwd()+r'\Rep_Days\num_rep_days='+str(nD)+'.csv');
    df = pd.read_csv(os.getcwd() + '/temporal_cluster/temporal_cluster_days=' +str(nD)+'_'+fName+'.csv');
    # df = pd.read_csv(os.getcwd() + '/temporal_cluster/'+fName+'_days='+str(nD)+'.csv');

    # s1 = df[df['Medoid']==1];
    # s2 = s1.sort_values(['Cluster']);
    # g_rep_days = np.array(s2.index);    
    # e_time_weight = np.zeros(24*nD,dtype=int);
    # e_rep_hrs = np.zeros(24*nD,dtype=int);
    # g_time_weight = np.zeros(nD,dtype=int);    
    # days_in_cluster = dict();
    # days2Medoid = dict();
    
    # # sort g_rep_days and g_weights
    # g_rep_days = sorted(g_rep_days);
    # for i in range(nD):
    #     s1 = np.array(df['Cluster']);
    #     i2 = df['Cluster'].iloc[g_rep_days[i]];
    #     s2 = np.where(s1==i2);
    #     g_time_weight[i]= len(s2[0]);
    #     for d in s2[0]:
    #         days2Medoid[d] = g_rep_days[i];
    #     days_in_cluster[i] = s2[0];
    #     e_time_weight[i*24:(i+1)*24] += g_time_weight[i];
    #     for j in range(24):
    #          e_rep_hrs[i*24+j] =g_rep_days[i]*24+j;
    # return e_time_weight,g_time_weight,g_rep_days,e_rep_hrs,days_in_cluster,days2Medoid;
    e_time_weight = np.zeros(24*nD,dtype=int);
    e_rep_hrs = np.zeros(24*nD,dtype=int);
    g_time_weight = np.zeros(nD,dtype=int);
    g_rep_days = np.array(df['Day of Year']);
    for i in range(nD):
        e_time_weight[i*24:(i+1)*24] += int(df['Weight'][i]);
        for j in range(24):
            e_rep_hrs[i*24+j] =g_rep_days[i]*24+j;
        g_time_weight[i]= int(df['Weight'][i]);
    return e_time_weight,g_time_weight,g_rep_days,e_rep_hrs;

# # Pow_dir = os.getcwd()+r'\Power_System_Data';
# # Gas_dir = os.getcwd()+r'\Gas_System_Data';

# Pow_dir = os.getcwd()+'/Power_System_Data';
# Gas_dir = os.getcwd()+'/Gas_System_Data';

# plant2sym = {0:'ng',1:'solar',2:'wind', 3:'hydro',4:'nuclear',5:'CT',6:'CC',7:'CC-CCS',8:'solar-UPV',9:'wind-new',10:'wind-offshore-new',11:'nuclear-new'};
# sym2plant = {'ng':0,'solar':1,'wind':2,'hydro':3,'nuclear':4,'CT':5,'CC':6,'CC-CCS':7,'solar-UPV':8,'wind-new':9,'wind-offshore-new':10,'nuclear-new':11};

# # this is the order in Breakthrough data
# zone_id2state = {0:'Maine', 1:'New Hampshire',2:'Vermont',
#                3:'Massachusetts',4:'Rhode Island',5:'Connecticut'};
# state2zone_id = {'Maine':0, 'New Hampshire':1,'Vermont':2,
#                'Massachusetts':3,'Rhode Island':4,'Connecticut':5};

# # df_ebus = pd.read_csv(Pow_dir+r'\Region_buses.csv');
# # df_eAdj = pd.read_csv(Pow_dir+r'\bus_adj_Nodes.csv');
# # df_ePlt = pd.read_csv(Pow_dir+r'\Region_plant.csv');
# # df_eDem = pd.read_csv(Pow_dir+r'\zonal_demand.csv');
# # df_op =  pd.read_csv(Pow_dir+r'\Other_params.csv');
# # df_solar = pd.read_csv(Pow_dir+r'\profile_solar_hourly.csv'); 
# # df_hydro = pd.read_csv(Pow_dir+r'\profile_hydro_hourly.csv'); 
# # df_wind = pd.read_csv(Pow_dir+r'\profile_wind_hourly.csv'); 
# # df_wind_offshore = pd.read_csv(Pow_dir+r'\profile_wind_offshore_hourly.csv'); 
# # df_reg_mult = pd.read_csv(Pow_dir+r'\Regional_multipliers.csv');
# # df_plt = pd.read_csv(Pow_dir+r'\Plants.csv');
# # df_br = pd.read_csv(Pow_dir+r'\Branches.csv');
# # df_str = pd.read_csv(Pow_dir+r'\storage.csv');


# # df_gnode = pd.read_csv(Gas_dir+r'/ng_nodes.csv');
# # df_ng_Lexp = pd.read_csv(Gas_dir+r'/ng_L_exp.csv');
# # df_ng_Limp = pd.read_csv(Gas_dir+r'/ng_L_imp.csv');
# # df_ng_adjE = pd.read_csv(Gas_dir+r'/ng_adjE.csv');
# # df_g2g_br = pd.read_csv(Gas_dir+r'/g2g_br.csv');
# # df_ng_dem = pd.read_csv(Gas_dir+r'/ng_daily_demand2016.csv');
# # df_exis_svl = pd.read_csv(Gas_dir+r'/ng_SVL_facilities.csv');
# # df_svl_params = pd.read_csv(Gas_dir+r'/SVL_params.csv');
# # df_ng_other = pd.read_csv(Gas_dir+r'/Other_params.csv');

# # total production capacity for various plant types (E3 EFI New England Net-zero Study)
# # Utility scale solar, wind, wind-offshore, nuclear (for all other types there is no limit)
# type_prod_lim = [22e3, 10e3, 280e3,3.5e3]; # 10e6 means no limit

# df_ebus = pd.read_csv(Pow_dir+'/Region_Nodes.csv');
# #df_eAdj = pd.read_csv(Pow_dir+'/bus_adj_Nodes-6-nodes.csv');
# df_ePlt = pd.read_csv(Pow_dir+'/Region_plants.csv');
# df_br = pd.read_csv(Pow_dir+'/Branches.csv');
# df_eDem = pd.read_csv(Pow_dir+'/zonal_load-RM.csv');
# df_eDem = df_eDem.iloc[:8760,:];
# df_op =  pd.read_csv(Pow_dir+'/Other_params.csv');
# df_solar = pd.read_csv(Pow_dir+'/profile_solar_hourly.csv'); 
# df_solar = df_solar.iloc[:8760,:];
# df_hydro = pd.read_csv(Pow_dir+'/profile_hydro_hourly.csv');
# df_hydro = df_hydro.iloc[:8760,:]; 
# df_wind = pd.read_csv(Pow_dir+'/profile_wind_hourly.csv'); 
# df_wind = df_wind.iloc[:8760,:];
# df_wind_offshore = pd.read_csv(Pow_dir+'/profile_wind_offshore_hourly.csv'); 
# df_wind_offshore = df_wind_offshore.iloc[:8760,:];
# df_reg_mult = pd.read_csv(Pow_dir+'/Regional_multipliers.csv');
# df_plt = pd.read_csv(Pow_dir+'/Plants.csv');
# df_str = pd.read_csv(Pow_dir+'/storage.csv');
# df_ccs = pd.read_csv(Pow_dir+'/CCS.csv');
# node2zone = np.arange(len(df_ebus));

# # if Setting.Power_network_size==88:
# #     df_ebus = pd.read_csv(Pow_dir+'/Region_buses-188-nodes.csv');
# #     #df_eAdj = pd.read_csv(Pow_dir+'/bus_adj_Nodes-188-nodes.csv');
# #     df_ePlt = pd.read_csv(Pow_dir+'/Region_plant-188-nodes.csv');
# #     df_br = pd.read_csv(Pow_dir+'/Branches-188-nodes.csv');
# #     df_eDem = pd.read_csv(Pow_dir+'/bus_demand.csv');
# #     node2zone = np.array(df_ebus['zone_id'])-1;
# df_gnode = pd.read_csv(Gas_dir+'/ng_nodes.csv');
# df_ng_Lexp = pd.read_csv(Gas_dir+'/ng_L_exp.csv');
# df_ng_Limp = pd.read_csv(Gas_dir+'/ng_L_imp.csv');      
# df_ng_adjE = pd.read_csv(Gas_dir+'/ng_adjE-6-nodes.csv');
# df_g2g_br = pd.read_csv(Gas_dir+'/g2g_br.csv');
# df_ng_dem = pd.read_csv(Gas_dir+'/ng_daily_load2050_'+Setting.electrification_scenario+'.csv');
# df_exis_svl = pd.read_csv(Gas_dir+'/ng_SVL_facilities.csv');
# df_svl_params = pd.read_csv(Gas_dir+'/SVL_params.csv');
# df_ng_other = pd.read_csv(Gas_dir+'/Other_params.csv');

# if Setting.Power_network_size==188:
#     df_ng_adjE = pd.read_csv(Gas_dir+'/ng_adjE-188-nodes.csv');
#%% other parameters
class other_prm():
    SVL_lifetime = 30;
    pipeline_lifetime = 30;
    WACC = 0.07;
    RNG_price = 20;
    trans_unit_cost = 3500;
    trans_line_lifespan = 30;
    NG_price = 5.45;
    Nuclear_price = 0.72;
    pipe_per_mile = 7e+5;
    pi = 3.141592;
    NG_emission = 0.053; # ton/MW



#%% enode
class enode:    
    # def __init__(self,n)    
    num = [];  # scalar
    adj_buses = np.array([]);
    Init_plt_type = np.array([]);
    Init_plt_count = np.array([]);
    demand  = np.array([]); 
    cap_factors = np.array([]);
    arcs = [];
    arc_sign= [];
    # @staticmethod
    # def typeThis():
    #     print("Hi,this tutorial is about static class in python");
        


    
#%% plants
class plant:
    Type = '';
    num = [];
    is_exist = [];
    capex = [];
    VOM = [];
    FOM = [];
    co2_capture_rate = [];
    heat_rate=[];
    lifetime= 0;
    decom_cost = [];
    nameplate_cap = [];
    min_output = [];
    ramp_rate= [];
    startup_cost=[];
    startup_fuel=[];
    est_cost_coef = np.array([]);
    emission=[];
    
    regional_mult = np.array([]);



#%% branch
class branch:
    from_node = [];
    to_node= [];
    suscept = [];
    maxFlow = [];
    length = [];
    is_exist = [];
    est_coef = [];


#%% power storage
class eStorage:
    energy_capex= [];
    power_capex = [];
    eff_ch = [];
    eff_disCh = [];
    eFOM = [];
    pFOM = [];
    lifetime = int();
    est_coef = [];
    


    
#%% CC-CCS data
class CCS:
    pipe_capex = float();
    str_capex = []; #float
    elec_req_pipe = [];
    elec_req_pump = [];
    comp_dis = [];
    str_loc = np.zeros(2);  #1x2 array
    node2str_dis = []; #1xlen(nodes) array



#%% Gas system: 
    
class gnode:
    num = [];
    fips = [];
    out_dem = [];
    demand = np.array([]);
    injU = [];
    L_exp = np.array([]);
    L_imp = np.array([]);
    adjE = np.array([]);
    adjS = np.array([]);



class pipe:
    from_node = [];
    to_node = [];
    is_exist = [];
    length = [];
    Cap = [];
    inv_coef =[];
    



class exist_SVL:
    num = [];
    str_cap = [];
    vap_cap = [];
    liq_cap = [];


class SVL:
    capex = [];
    FOM = [];
    eff_ch = [];
    eff_disCh = [];
    BOG = [];
    inv_coef = [];
    
class data:    
    Enodes=[];Gnodes=[];Branches=[];node2zone=[];
    PipeLines=[];Plants=[];eStore=[];Other_input=[];state2zone_id=[];plant2sym=[];
    sym2plant=[]; time_weights=[];zone_id2state=[];Exist_SVL=[];SVLs=[];
    type_prod_lim=[];


def Run_problem_data():
    Other_input = other_prm();
    Pow_dir = os.getcwd()+'/Power_System_Data';
    Gas_dir = os.getcwd()+'/Gas_System_Data';
    
    plant2sym = {0:'ng',1:'solar',2:'wind', 3:'hydro',4:'nuclear',5:'CT',6:'CC',7:'CC-CCS',8:'solar-UPV',9:'wind-new',10:'wind-offshore-new',11:'nuclear-new'};
    sym2plant = {'ng':0,'solar':1,'wind':2,'hydro':3,'nuclear':4,'CT':5,'CC':6,'CC-CCS':7,'solar-UPV':8,'wind-new':9,'wind-offshore-new':10,'nuclear-new':11};
    
    # this is the order in Breakthrough data
    zone_id2state = {0:'Maine', 1:'New Hampshire',2:'Vermont',
               3:'Massachusetts',4:'Rhode Island',5:'Connecticut'};
    state2zone_id = {'Maine':0, 'New Hampshire':1,'Vermont':2,
               'Massachusetts':3,'Rhode Island':4,'Connecticut':5};
    
    type_prod_lim = [22e3, 10e3, 280e3,3.5e3]; # 10e6 means no limit
    
    df_ebus = pd.read_csv(Pow_dir+'/Region_Nodes.csv');
    #df_eAdj = pd.read_csv(Pow_dir+'/bus_adj_Nodes-6-nodes.csv');
    df_ePlt = pd.read_csv(Pow_dir+'/Region_plants.csv');
    df_br = pd.read_csv(Pow_dir+'/Branches.csv');
    df_eDem = pd.read_csv(Pow_dir+'/zonal_load-RM.csv');
    df_eDem = df_eDem.iloc[:8760,:];
    df_op =  pd.read_csv(Pow_dir+'/Other_params.csv');
    df_solar = pd.read_csv(Pow_dir+'/profile_solar_hourly.csv'); 
    df_solar = df_solar.iloc[:8760,:];
    df_hydro = pd.read_csv(Pow_dir+'/profile_hydro_hourly.csv');
    df_hydro = df_hydro.iloc[:8760,:]; 
    df_wind = pd.read_csv(Pow_dir+'/profile_wind_hourly.csv'); 
    df_wind = df_wind.iloc[:8760,:];
    df_wind_offshore = pd.read_csv(Pow_dir+'/profile_wind_offshore_hourly.csv'); 
    df_wind_offshore = df_wind_offshore.iloc[:8760,:];
    df_reg_mult = pd.read_csv(Pow_dir+'/Regional_multipliers.csv');
    df_plt = pd.read_csv(Pow_dir+'/Plants.csv');
    df_str = pd.read_csv(Pow_dir+'/storage.csv');
    df_ccs = pd.read_csv(Pow_dir+'/CCS.csv');
    node2zone = np.arange(len(df_ebus));
    
    # if Setting.Power_network_size==88:
    #     df_ebus = pd.read_csv(Pow_dir+'/Region_buses-188-nodes.csv');
    #     #df_eAdj = pd.read_csv(Pow_dir+'/bus_adj_Nodes-188-nodes.csv');
    #     df_ePlt = pd.read_csv(Pow_dir+'/Region_plant-188-nodes.csv');
    #     df_br = pd.read_csv(Pow_dir+'/Branches-188-nodes.csv');
    #     df_eDem = pd.read_csv(Pow_dir+'/bus_demand.csv');
    #     node2zone = np.array(df_ebus['zone_id'])-1;
    df_gnode = pd.read_csv(Gas_dir+'/ng_nodes.csv');
    df_ng_Lexp = pd.read_csv(Gas_dir+'/ng_L_exp.csv');
    df_ng_Limp = pd.read_csv(Gas_dir+'/ng_L_imp.csv');      
    df_ng_adjE = pd.read_csv(Gas_dir+'/ng_adjE-6-nodes.csv');
    df_g2g_br = pd.read_csv(Gas_dir+'/g2g_br.csv');
    df_ng_dem = pd.read_csv(Gas_dir+'/ng_daily_load2050_'+Setting.electrification_scenario+'.csv');
    df_exis_svl = pd.read_csv(Gas_dir+'/ng_SVL_facilities.csv');
    df_svl_params = pd.read_csv(Gas_dir+'/SVL_params.csv');
    df_ng_other = pd.read_csv(Gas_dir+'/Other_params.csv');
    

    Enodes = list();
    
    # for zonal demand
    # growth_rate =  df_op['demand_growth_rate'][0];
    # yrs = df_op['target_year'][0]-df_op['base_year'][0];
    # coef = (1+growth_rate)**yrs;
    
    for i in range(len(df_ebus)):
        en = enode();
        en.num = i;
        #s1 = np.array(df_eAdj.iloc[i,:]);
        #s1 = s1[np.logical_not(np.isnan(s1))];
        #en.adj_buses = s1.astype(int);
        en.Init_plt_type = sym2plant.keys();
        en.Init_plt_count = np.zeros(len(sym2plant.keys()),dtype=int);    
        s1 = np.array(df_eDem.iloc[:,i+1]); # directly using 2050 load
        en.demand = s1;
        en.cap_factors = np.ones((8760,len(sym2plant.keys())));
    
        
        en.cap_factors[:,sym2plant['solar']] = df_solar.iloc[:,i+1];
        en.cap_factors[:,sym2plant['solar-UPV']] = df_solar.iloc[:,i+1];
    
        en.cap_factors[:,sym2plant['wind']] = df_wind.iloc[:,i+1];
        en.cap_factors[:,sym2plant['wind-new']] = df_wind.iloc[:,i+1];
        #en.cap_factors[:,sym2plant['hydro']] = df_hydro.iloc[:,i];
        en.cap_factors[:,sym2plant['wind-offshore-new']] = df_wind_offshore.iloc[:,0];
        # else:
        #     s1 = int(df_ebus['zone_id'][i])-1;
        #     en.cap_factors[:,sym2plant['solar']] = df_solar.iloc[:,s1];
        #     en.cap_factors[:,sym2plant['solar-UPV']] = df_solar.iloc[:,s1];
        #     en.cap_factors[:,sym2plant['wind']] = df_wind.iloc[:,s1];
        #     en.cap_factors[:,sym2plant['wind-new']] = df_wind.iloc[:,s1];
        #     #en.cap_factors[:,sym2plant['hydro']] = df_hydro.iloc[:,i];
        #     en.cap_factors[:,sym2plant['wind-offshore-new']] = df_wind_offshore.iloc[:,0];
        
        Enodes.append(en);
        
    for i in range(len(df_ePlt)):
        s1  = df_ePlt['plant_type'][i];
        # if Setting.Power_network_size!=88:
        #     s2 = df_ePlt['zone_id'][i]-1;
        # else:
        s2 = int(df_ePlt['zone_id'][i])-1;
        if s1 in sym2plant.keys():
            if Setting.Power_network_size!=88:
                Enodes[s2].Init_plt_count[sym2plant[s1]] = int(df_ePlt['adjusted_count'][i]);
            else:
                Enodes[s2].Init_plt_count[sym2plant[s1]] = int(df_ePlt['count'][i]);
                
    
    Plants = list();
    for i in range(len(sym2plant.keys())):
        plt = plant();
        plt.Type = plant2sym[i];
        plt.num = i;
        plt.is_exist = df_plt['is existing'][i];
        plt.capex = df_plt['CAPEX per plant'][i];
        plt.VOM = df_plt['VOM ($/MWh)'][i];
        plt.co2_capture_rate = df_plt['Carbon capture rate'][i];
        plt.heat_rate = df_plt['Heat Rate  (MMBtu/MWh)'][i];
        plt.lifetime = df_plt['Lifetime (year)'][i];
        plt.decom_cost = df_plt['Decom. cost ($) per plant'][i];
        plt.nameplate_cap = df_plt['Nameplate capacity (MW)'][i];
        plt.FOM = df_plt['FOM ($/kW-yr)'][i]*plt.nameplate_cap*1000;
        plt.min_output = df_plt['Minimum stable output (%)'][i];
        plt.ramp_rate = df_plt['Hourly Ramp rate (%)'][i];
        plt.startup_cost = df_plt['Startup Cost (per  plant)'][i];
        plt.startup_fuel = df_plt['Startup Fuel (MMBtu)'][i];
        plt.emission = (1-plt.co2_capture_rate)*Other_input.NG_emission;
        if plt.is_exist==0:
            plt.regional_mult = df_reg_mult.iloc[:,i-4];
            s1 = (1/(1+Other_input.WACC)**plt.lifetime);
            plt.est_cost_coef = (Other_input.WACC/(1-s1))* plt.regional_mult;
        # else:
        #     plt.regional_mult = np.zeros(len(state2zone_id.keys()));
        #     plt.est_cost_coef = np.zeros(len(state2zone_id.keys()));
        
        
        Plants.append(plt);
        
    Branches = [];
    
    arcs = [[] for x in range(len(Enodes))];
    arc_sign = [[] for x in range(len(Enodes))];
    
    for b in range(len(df_br)):
        br = branch();
        br.from_node = int(df_br['from bus'][b]);
        br.to_node = int(df_br['to bus'][b]);
        
        arcs[br.from_node].append(b);
        arcs[br.to_node].append(b);
        if br.from_node>br.to_node:
            arc_sign[br.from_node].append(-1);
            arc_sign[br.to_node].append(1);
        else:
            arc_sign[br.from_node].append(1);
            arc_sign[br.to_node].append(-1);
        
        br.suscept = df_br['susceptance'][b];
        br.maxFlow = df_br['maxFlow'][b];
        br.length = df_br['distance'][b];
        br.is_exist = int(df_br['is_exist'][b]);
        s1 = (1/(1+Other_input.WACC)**Other_input.trans_line_lifespan);
        br.est_coef = (Other_input.WACC/(1-s1));
        
        Branches.append(br);
        
    for i in range(len(arcs)):
        Enodes[i].arcs = np.array(arcs[i]);
        Enodes[i].arc_sign = np.array(arc_sign[i]);
    
    eStore = list();
    
    for i in range(len(df_str)):
        st = eStorage();
        st.energy_capex = df_str['energy capex'][i];
        st.power_capex = df_str['power capex'][i];    
        st.eff_ch = df_str['charging efficiency'][i];
        st.eff_disCh = df_str['discharging efficiency'][i];
        st.eFOM = df_str['energy FOM'][i];
        st.pFOM = df_str['power FOM'][i];
        st.lifetime = int(df_str['lifetime'][i]);
        s1 = (1/(1+Other_input.WACC)**st.lifetime);
        st.est_coef = (Other_input.WACC/(1-s1));
        eStore.append(st);
        
    CC_CCS = CCS();
    CC_CCS.str_loc = [df_ccs['str_loc_lat'][0],df_ccs['str_loc_lon'][0]];
    CC_CCS.pipe_capex = df_ccs['lev_inv_pipe'][0];
    CC_CCS.str_capex = df_ccs['lev_inv_str'][0];
    CC_CCS.elec_req_pipe = df_ccs['E_pipe'][0];
    CC_CCS.elec_req_pump = df_ccs['E_pump'][0];
    CC_CCS.comp_dis = df_ccs['comp_dis'][0];
    
    CC_CCS.node2str_dis = np.zeros(len(df_ebus));
    CC_CCS.node2str_dis=np.array([433,338,294,302,297,232]);
    # for i in range(len(df_ebus)):
    #     lat_lon = [df_ebus['lat'][i],df_ebus['lon'][i]];
    #     CC_CCS.node2str_dis[i] = distance(CC_CCS.str_loc,lat_lon).miles;
    
    Gnodes = list();
    # for zonal demand
    # growth_rate =  df_ng_other['demand growth rate'][0];
    # yrs = df_op['target_year'][0]-df_op['base_year'][0];
    # coef = (1+growth_rate)**yrs;
    ng_dem=0;
    for i in range(len(df_gnode)):
        gas = gnode();
        gas.num = i;
        gas.fips = int(df_gnode['FIPS'][i]);
        gas.out_dem = df_gnode['out_of_state_demand'][i];
        gas.injU = df_gnode['inj_capacity'][i];
        gas.adjS =np.array([int(df_gnode['SVL1'][i]),int(df_gnode['SVL2'][i])]);
        
        s1 = np.array(df_ng_Lexp.iloc[i,:]);
        s1 = s1[np.logical_not(np.isnan(s1))];
        gas.L_exp = s1.astype(int);
        
        s1 = np.array(df_ng_Limp.iloc[i,:]);
        s1 = s1[np.logical_not(np.isnan(s1))];
        gas.L_imp = s1.astype(int);
        
        s1 = np.array(df_ng_adjE.iloc[i,:]);
        s1 = s1[np.logical_not(np.isnan(s1))];
        gas.adjE = s1.astype(int);
        gas.demand = df_ng_dem.iloc[:,i];
        ng_dem += np.sum(gas.demand);
        Gnodes.append(gas);
    
    
    PipeLines = list();
    for i in range(len(df_g2g_br)):
        pp = pipe();
        pp.from_node = df_g2g_br['from node'][i];
        pp.to_node = df_g2g_br['to node'][i];
        pp.is_exist = df_g2g_br['is exist'][i];
        pp.length = df_g2g_br['length'][i];
        pp.Cap = df_g2g_br['max capacity'][i];
        
        s1 = (1/(1+Other_input.WACC)**Other_input.pipeline_lifetime);
        pp.inv_coef = (Other_input.WACC/(1-s1));
        PipeLines.append(pp);
        
    Exist_SVL = list();
    for i in range(len(df_exis_svl)):
        eSvl = exist_SVL();
        eSvl.num = i;
        eSvl.str_cap = df_exis_svl['Storage-cap'][i];
        eSvl.vap_cap = df_exis_svl['Vap-cap'][i];
        eSvl.liq_cap = df_exis_svl['Liq-cap'][i];
        
        Exist_SVL.append(eSvl);
    
        
    
    SVLs = list();
    for i in range(len(df_svl_params)):
        svl = SVL();
        svl.capex = df_svl_params['capex'][i];
        svl.FOM = df_svl_params['FOM'][i];
        svl.eff_ch = df_svl_params['eff_ch'][i];
        svl.eff_disCh = df_svl_params['eff_disCh'][i];
        svl.BOG = df_svl_params['BOG'][i];
        s1 = (1/(1+Other_input.WACC)**Other_input.SVL_lifetime);
        svl.inv_coef = (Other_input.WACC/(1-s1));
    
        SVLs.append(svl);
    Data = data();
    Data.Enodes= Enodes;
    Data.Gnodes = Gnodes;
    Data.Branches = Branches;
    Data.node2zone= node2zone;
    Data.PipeLines = PipeLines;
    Data.Plants = Plants;
    Data.eStore = eStore;
    Data.Other_input = Other_input;
    Data.state2zone_id = state2zone_id;
    Data.plant2sym = plant2sym;
    Data.sym2plant = sym2plant;
    Data.zone_id2state = zone_id2state;
    Data.Exist_SVL = Exist_SVL;
    Data.SVLs = SVLs;
    Data.type_prod_lim = type_prod_lim;
    
    return Data;

#%% delete unnecessary data from the workspace
# del df_br,df_ebus,df_eDem,df_ePlt,df_exis_svl,df_g2g_br,df_gnode;
# del df_hydro,df_ng_adjE,df_ng_dem,df_ng_Lexp,df_ng_Limp,df_ng_other,df_op;
# del df_plt,df_reg_mult,df_solar,df_str,df_svl_params,df_wind,df_wind_offshore;
# del s1,s2,svl,pp,Pow_dir,Gas_dir,plt,i,st,en;

# del br,gas,enode;
# del sym2plant,plant2sym,zone_id2state,eSvl
