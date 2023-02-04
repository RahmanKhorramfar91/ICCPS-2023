# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 14:29:48 2022

@author: Rahman Khorramfar
"""
from Setting import Setting,EV,GV;
import os;import time;import csv;
from ProblemData_UB import Rub_problem_data, time_weights;
# from ProblemData import Enodes,Gnodes,Branches,node2zone;
# from ProblemData import PipeLines,Plants,eStore,Other_input,state2zone_id,plant2sym;
# from ProblemData import sym2plant,time_weights,zone_id2state,Exist_SVL,SVLs;
# from ProblemData import type_prod_lim;

import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
import numpy as np;
import pandas as pd;


e_time_weight,g_time_weight,g_rep_days,e_rep_hrs,days_in_cluster,days2Medoid = time_weights(Setting.num_rep_days,Setting.rep_day_folder);

# nE = len(Enodes); 
# nG = len(Gnodes);
# nBr = len(Branches);
# nPipe = len(PipeLines);
# nPlt = len(Plants);
# neSt = len(eStore);
# Te = range(len(e_time_weight));
# Tg = range(len(g_time_weight));
# nPipe = len(PipeLines);
# nSVL = len(Exist_SVL);
# nG = len(Gnodes);

thermal_units = ["ng","CT","CC","CC-CCS","nuclear","nuclear-new"];
   # ng_units = ['ng','solar','wind','hydro','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore-new'];
    # #ng_units = ['ng','solar','wind','hydro','CT','CC','CC-CCS','solar-UPV','wind-new','nuclear','wind-offshore-new'];


def Power_System_Model(Model,Data): 
# fetch data     
    Enodes=  Data.Enodes;
    Gnodes = Data.Gnodes;
    Branches = Data.Branches;
    node2zone= Data.node2zone;
    PipeLines = Data.PipeLines;
    Plants = Data.Plants;
    eStore = Data.eStore;
    Other_input = Data.Other_input;
    state2zone_id = Data.state2zone_id;
    plant2sym = Data.plant2sym;
    sym2plant = Data.sym2plant;
    zone_id2state = Data.zone_id2state;
    Exist_SVL = Data.Exist_SVL;
    SVLs = Data.SVLs;
    type_prod_lim = Data.type_prod_lim;
    nE = len(Enodes); 
    nG = len(Gnodes);
    nBr = len(Branches);
    nPipe = len(PipeLines);
    nPlt = len(Plants);
    neSt = len(eStore);
    Te = range(len(e_time_weight));
    Tg = range(len(g_time_weight));
    nPipe = len(PipeLines);
    nSVL = len(Exist_SVL);
    nG = len(Gnodes);
    
    
    
    #% define decision variables 
    EV.Xop = Model.addVars(nE,nPlt,vtype=GRB.INTEGER);
    EV.Xest = Model.addVars(nE,nPlt,vtype=GRB.INTEGER);
    EV.Xdec = Model.addVars(nE,nPlt,vtype=GRB.INTEGER);
    EV.Ze = Model.addVars(nBr,vtype=GRB.BINARY);
    
    if Setting.UC_active:
        EV.X = Model.addVars(nE,len(Te), nPlt,vtype=GRB.INTEGER);
        EV.Xup= Model.addVars(nE,len(Te), nPlt,vtype=GRB.INTEGER);
        EV.Xdown = Model.addVars(nE,len(Te), nPlt,vtype=GRB.INTEGER);
    
    if Setting.relax_int_vars:
        EV.Xop = Model.addVars(nE,nPlt,vtype=GRB.CONTINUOUS);
        EV.Xest = Model.addVars(nE,nPlt,vtype=GRB.CONTINUOUS);
        EV.Xdec = Model.addVars(nE,nPlt,vtype=GRB.CONTINUOUS);
        EV.Ze = Model.addVars(nBr,vtype=GRB.CONTINUOUS);
        
    if (Setting.relax_UC_vars or Setting.relax_int_vars) and (Setting.UC_active):
        EV.X = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
        EV.Xup= Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
        EV.Xdown = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    
    
    
    EV.prod= Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    EV.theta = Model.addVars(nE,len(Te),lb=np.zeros((nE,len(Te)))-GRB.INFINITY,ub=np.zeros((nE,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    EV.Shed = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    
    EV.YeCD = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    EV.YeLev = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    #EV.YeStr = Model.addVars(nE,neSt,vtype=GRB.BINARY);
    EV.kappa_capt = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    EV.kappa_pipe = Model.addVars(nE,vtype=GRB.CONTINUOUS);
    
    EV.eSch =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSdis =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSlev =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    
    EV.flowE =  Model.addVars(nBr,len(Te),lb=np.zeros((nBr,len(Te)))-GRB.INFINITY,ub=np.zeros((nBr,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    
    
    EV.est_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.est_trans_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.decom_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.FOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.startup_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.VOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.nuc_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.gas_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.shedding_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.elec_storage_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);

    EV.e_system_cost=Model.addVar(vtype=GRB.CONTINUOUS);
    
    
    #% Set some variables
    #1) existing types can not be established because there are new equivalent types
    #2) new types cannot be decommissioned
    #3) offshore only allowed in Massachusetts and Connecticut
    
    s1 = [state2zone_id['Massachusetts'],state2zone_id['Connecticut']];
    #Model.addConstrs(EV.Xest[n,i]==0 for n in range(nE) for i in range(nPlt) if n not in s1 if Plants[i].Type=='wind-offshore-new');
    Model.addConstrs(EV.Xest[n,i]==0 for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==1)
    Model.addConstrs(EV.Xdec[n,i]==0 for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==0)
    # for n in range(nE):
    #     #s1 = int(node2zone[n]);
    #     #if not(zone_id2state[s1]=='Massachusetts' or zone_id2state[s1]=='Connecticut'):
    #     Model.addConstr(EV.Xest[n,sym2plant['wind-offshore-new']]==0);


    #% Electricity System Objective Function    
    # cost components
    est_cost = LinExpr(quicksum(Plants[i].capex*EV.Xest[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==0));
    dec_cost = LinExpr(quicksum(Plants[i].decom_cost*EV.Xdec[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==1));
    fom_cost = LinExpr(quicksum(Plants[i].FOM*EV.Xop[n,i] for n in range(nE) for i in range(nPlt)));
    vom_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].VOM*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    nuc_fuel_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].heat_rate *Other_input.Nuclear_price*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if (Plants[i].Type=="nuclear" or Plants[i].Type=="nuclear-new") ));
    NG_units = ["ng","CT","CC","CC-CCS"];
    gas_fuel_cost = 0;
    gas_fuel_cost = LinExpr(quicksum(Other_input.NG_price*e_time_weight[t]*Plants[i].heat_rate*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if plant2sym[i] in NG_units))
    startup_cost = 0;
    if Setting.UC_active:
        startup_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].startup_cost*EV.Xup[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    shed_cost = LinExpr(quicksum(e_time_weight[t]*Setting.e_shed_penalty*EV.Shed[n,t] for n in range(nE) for t in Te));
    strg_cost = LinExpr(quicksum(eStore[r].est_coef*(eStore[r].power_capex*EV.YeCD[n,r]+eStore[r].energy_capex*EV.YeLev[n,r]) for n in range(nE) for r in range(neSt)));
    strg_cost += LinExpr(quicksum((eStore[r].pFOM*EV.YeCD[n,r]+eStore[r].eFOM*EV.YeLev[n,r]) for n in range(nE) for r in range(neSt)));
    tran_cost = LinExpr(quicksum(Branches[b].est_coef*Other_input.trans_unit_cost*Branches[b].maxFlow*Branches[b].length*EV.Ze[b] for b in range(nBr)));
    
    # total power system cost function
    if Setting.emis_case==1: # add gas fuel cost only if Case=1
        e_total_cost = gas_fuel_cost+est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+startup_cost+shed_cost+strg_cost+tran_cost;
    else:
        e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+startup_cost+shed_cost+strg_cost+tran_cost;

    #e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+shed_cost;
    
    Model.addConstr(EV.est_cost == est_cost);
    Model.addConstr(EV.decom_cost == dec_cost);
    Model.addConstr(EV.FOM_cost == fom_cost);
    Model.addConstr(EV.VOM_cost == vom_cost);
    Model.addConstr(EV.nuc_fuel_cost == nuc_fuel_cost);
    Model.addConstr(EV.startup_cost == startup_cost);
    Model.addConstr(EV.shedding_cost == shed_cost);
    Model.addConstr(EV.elec_storage_cost == strg_cost);
    Model.addConstr(EV.est_trans_cost==tran_cost);
    Model.addConstr(EV.gas_fuel_cost == gas_fuel_cost);
    Model.addConstr(EV.e_system_cost==e_total_cost);    
    
    #% Electricity System Constraints
    # C1: number of generation units at each node
    Model.addConstrs(EV.Xop[n,i] == Enodes[n].Init_plt_count[i]+EV.Xest[n,i]-EV.Xdec[n,i] for n in range(nE) for i in range(nPlt));
        
    # C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)    
    Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
    if Setting.UC_active:
        # UC:
        Model.addConstrs(EV.X[n,t,i]-EV.X[n,t-1,i]==EV.Xup[n,t,i]-EV.Xdown[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units));
        # ramping:
        Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*(EV.X[n,t,i]-EV.Xup[n,t,i])+max(Plants[i].min_output, Plants[i].ramp_rate)*Plants[i].nameplate_cap*EV.Xup[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*(EV.X[n,t,i]-EV.Xup[n,t,i])+max(Plants[i].min_output, Plants[i].ramp_rate)*Plants[i].nameplate_cap*EV.Xup[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        # production limit
        Model.addConstrs(EV.prod[n,t,i]>= Plants[i].min_output*Plants[i].nameplate_cap*EV.X[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
        Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap*EV.X[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
        
        Model.addConstrs(EV.X[n,t,i] <= EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
    else:
        Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
    
    # C5, C6: flow limit for electricity
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1))
        
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow*EV.Ze[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow*EV.Ze[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0))
    
    # C7: power balance      
    if Setting.copper_plate_approx:
        Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt) for n in range(nE))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt) for n in range(nE))+quicksum(EV.Shed[n,t] for n in range(nE)) == quicksum(Enodes[n].demand[e_rep_hrs[t]] for n in range(nE)) for t in Te);
    else:        
        Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt))-quicksum(Enodes[n].arc_sign[b]*EV.flowE[Enodes[n].arcs[b],t] for b in range(len(Enodes[n].arcs)))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt))+EV.Shed[n,t] == Enodes[n].demand[e_rep_hrs[t]] for n in range(nE) for t in Te);

    
    # C8: flow equation
    # Model.addConstrs(EV.flowE[b,t]==Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t]) for b in range(nBr) for t in Te if Branches[b].is_exist==1);
    # Model.addConstrs(EV.flowE[b,t]-Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1-EV.Ze[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
    # Model.addConstrs(-EV.flowE[b,t]+Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1-EV.Ze[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
        
    # # C9: phase angle (theta) limits. already applied in the definition of the variable
    # Model.addConstrs(EV.theta[n,t]<=Other_input.pi for n in range(nE) for t in Te);
    # Model.addConstrs(-EV.theta[n,t]<=Other_input.pi for n in range(nE) for t in Te);
    # Model.addConstrs(EV.theta[0,t]==0 for t in Te);
    
    # C10: VRE production according to capacity factors
    Model.addConstrs(EV.prod[n,t,i]<=Enodes[n].cap_factors[e_rep_hrs[t],i]* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te);
    
    # C11: demand curtainlment constraint
    Model.addConstrs(EV.Shed[n,t]<= Enodes[n].demand[e_rep_hrs[t]] for  n in range(nE) for t in Te);
    
    # C12: RPS constraints
    VRE = ["solar","wind","wind-offshore-new","solar-UPV","wind-new"];# hydro not included
    Model.addConstr(quicksum(e_time_weight[t]*EV.prod[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if(plant2sym[i] in VRE)) >= Setting.VRE_share*quicksum(e_time_weight[t]*Enodes[n].demand[e_rep_hrs[t]] for n in range(nE) for t in Te));
    
    # C14,C15,C16 storage constraints
    Model.addConstrs(EV.eSlev[n,0,r]==eStore[r].eff_ch*EV.eSch[n,0,r]-EV.eSdis[n,0,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt));
    Model.addConstrs(EV.eSlev[n,t,r]-EV.eSlev[n,t-1,r]==eStore[r].eff_ch*EV.eSch[n,t,r]-EV.eSdis[n,t,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt) for t in Te if t>0);
    
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSdis[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSch[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeLev[n,r]>= EV.eSlev[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    
    # start and ending of storage should be the same in case of using rep. days
    if len(g_rep_days)<364:
        for k in range(len(g_rep_days)):
            s1 = k*24; s2=(k+1)*24-1;
            Model.addConstrs(EV.eSlev[n,s1,r]==EV.eSlev[n,s2,r] for n in range(nE) for r in range(neSt));
    
    #sym2plant = {'ng':0,'solar':1,'wind':2,'hydro':3,'nuclear':4,'CT':5,'CC':6,'CC-CCS':7,'solar-UPV':8,'wind-new':9,'wind-offshore-new':10,'nuclear-new':11};
    # Utility scale solar, wind, wind-offshore, nuclear (for all other types there is no limit)
    #type_prod_lim = [22e3, 10e3, 280e3,3.5e3]; # 10e6 means no limit
    # productin capacity for each new type    
    Model.addConstr(quicksum(Plants[sym2plant['solar-UPV']].nameplate_cap*EV.Xop[n,sym2plant['solar-UPV']]+Plants[sym2plant['solar']].nameplate_cap*EV.Xop[n,sym2plant['solar']] for n in range(nE) )<= type_prod_lim[0]);
    Model.addConstr(quicksum(Plants[sym2plant['wind-new']].nameplate_cap*EV.Xop[n,sym2plant['wind-new']]+Plants[sym2plant['wind']].nameplate_cap*EV.Xop[n,sym2plant['wind']] for n in range(nE) )<= type_prod_lim[1]);
    Model.addConstr(quicksum(Plants[sym2plant['wind-offshore-new']].nameplate_cap*EV.Xop[n,sym2plant['wind-offshore-new']] for n in range(nE) )<= type_prod_lim[2]);
    Model.addConstr(quicksum(Plants[sym2plant['nuclear-new']].nameplate_cap*EV.Xop[n,sym2plant['nuclear-new']]+Plants[sym2plant['nuclear']].nameplate_cap*EV.Xop[n,sym2plant['nuclear']] for n in range(nE) )<= type_prod_lim[3]);


    # CSS constraints
    #Model.addConstrs(EV.kappa_capt[n,t]==Other_input.NG_emission*(Plants[sym2plant['CC-CCS']].co2_capture_rate)*Plants[sym2plant['CC-CCS']].heat_rate*EV.prod[n,t,sym2plant['CC-CCS']] for n in range(nE) for t in Te);    
    #Model.addConstrs(EV.kappa_pipe[n] >= EV.kappa_capt[n,t] for n in range(nE) for t in Te);
    

def NG_System_Model(Model,Data):
    # fetch data     
    Enodes=  Data.Enodes;
    Gnodes = Data.Gnodes;
    Branches = Data.Branches;
    node2zone= Data.node2zone;
    PipeLines = Data.PipeLines;
    Plants = Data.Plants;
    eStore = Data.eStore;
    Other_input = Data.Other_input;
    state2zone_id = Data.state2zone_id;
    plant2sym = Data.plant2sym;
    sym2plant = Data.sym2plant;
    zone_id2state = Data.zone_id2state;
    Exist_SVL = Data.Exist_SVL;
    SVLs = Data.SVLs;
    type_prod_lim = Data.type_prod_lim;
    nE = len(Enodes); 
    nG = len(Gnodes);
    nBr = len(Branches);
    nPipe = len(PipeLines);
    nPlt = len(Plants);
    neSt = len(eStore);
    Te = range(len(e_time_weight));
    Tg = range(len(g_time_weight));
    nPipe = len(PipeLines);
    nSVL = len(Exist_SVL);
    nG = len(Gnodes);
    
    
    GV.Zg = Model.addVars(nPipe,vtype=GRB.BINARY);
    if Setting.relax_int_vars:
        GV.Zg = Model.addVars(nPipe,vtype=GRB.CONTINUOUS);
    GV.Xvpr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Xstr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Sstr = Model.addVars(nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.Svpr = Model.addVars(nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.Sliq = Model.addVars(nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.supply = Model.addVars(nG,len(Tg),vtype = GRB.CONTINUOUS);
    GV.Shed =  Model.addVars(nG,len(Tg),vtype = GRB.CONTINUOUS);
    GV.RNG_use = Model.addVars(nG,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowGG =  Model.addVars(nPipe,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowGE =  Model.addVars(nG,nE,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowGL =  Model.addVars(nG,nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowVG =  Model.addVars(nSVL,nG,len(Tg),vtype = GRB.CONTINUOUS);
    
    GV.inv_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.inv_pipe_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.shed_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.RNG_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.fom_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.import_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.g_system_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    
    # NG System Objective Function
    s1 = LinExpr(quicksum(PipeLines[b].inv_coef*PipeLines[b].length*Other_input.pipe_per_mile*GV.Zg[b] for b in range(nPipe)));   
    s2 = LinExpr(quicksum(g_time_weight[tau]*GV.supply[k,tau]*Other_input.NG_price for k in range(nG) for tau in Tg));
    s3 = LinExpr(quicksum(SVLs[0].inv_coef*(SVLs[0].capex*GV.Xstr[j]+SVLs[1].capex*GV.Xvpr[j]) for j in range(nSVL)));
    s4 = LinExpr(quicksum(SVLs[1].FOM*(Exist_SVL[j].vap_cap+GV.Xvpr[j])+SVLs[0].FOM*(Exist_SVL[j].str_cap+GV.Xstr[j]) for j in range(nSVL)));
    s5 = LinExpr(quicksum(g_time_weight[tau]*Setting.g_shed_penalty*GV.Shed[k,tau] for k in range(nG) for tau in Tg));
    s6 = LinExpr(quicksum(g_time_weight[tau]*Other_input.RNG_price*GV.RNG_use[k,tau]for k in range(nG) for tau in Tg));
    
    
    s7 = s1+s2+s3+s4+s5+s6;
    Model.addConstr(GV.inv_pipe_cost==s1);
    Model.addConstr(GV.import_cost==s2);
    Model.addConstr(GV.inv_str_cost==s3);
    Model.addConstr(GV.fom_str_cost==s4);
    Model.addConstr(GV.shed_cost==s5);
    Model.addConstr(GV.RNG_cost==s6);
    Model.addConstr(GV.g_system_cost==s7);
    
    # NG System Constraints
    #C1, C2: flow limit for NG
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap for i in range(nPipe) for tau in Tg if PipeLines[i].is_exist==1);
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap*GV.Zg[i] for i in range(nPipe) for tau in Tg if PipeLines[i].is_exist==0);
    
    # C3: flow balance, NG node (no shedding allowed as RNG is considered)
    # Model.addConstrs(GV.supply[k,tau]
    #                   -quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_exp)
    #                   +quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_imp)
    #                   -quicksum(GV.flowGE[k,n,tau] for n in Gnodes[k].adjE)
    #                   +quicksum((GV.flowVG[j,k,tau]-GV.flowGL[k,j,tau]) for j in Gnodes[k].adjS)
    #                   +GV.Shed[k,tau]
    #                   +GV.RNG_use[k,tau]
    #                  == Gnodes[k].demand[g_rep_days[tau]] for k in range(nG) for tau in Tg);                             
    Model.addConstrs(GV.supply[k,tau]
                      -quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_exp)
                      +quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_imp)
                      -quicksum(GV.flowGE[k,n,tau] for n in Gnodes[k].adjE)                      
                      +GV.Shed[k,tau]
                      +GV.RNG_use[k,tau]
                     == Gnodes[k].demand[g_rep_days[tau]] for k in range(nG) for tau in Tg);                             

    # C3,C4: injection (supply) limit and curtailment limit
    Model.addConstrs(GV.supply[k,tau]<= Gnodes[k].injU for k in range(nG) for tau in Tg);    
    Model.addConstrs(GV.Shed[k,tau]+GV.RNG_use[k,tau]<= Gnodes[k].demand[g_rep_days[tau]] for k in range(nG) for tau in Tg);
    
    # no constrain on RNG use
    #Model.addConstr(quicksum(g_time_weight[tau]*GV.RNG_use[k,tau] for k in range(nG) for tau in Tg)<= Setting.RNG_cap*Setting.poss_gas_consump);
    
    # C5: storage balance (storage strats empty)
    Model.addConstrs(GV.Sstr[j,tau]==Exist_SVL[j].str_cap*0+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in Tg if tau==0);
    Model.addConstrs(GV.Sstr[j,tau]==(1-SVLs[0].BOG)*GV.Sstr[j,tau-1]+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in Tg if tau>0);
    
    # C6,8: calculate Sliq, Svpr
    # for j in range(nSVL):
    #     for tau in Tg:
    #         NG_adj = [];
    #         for k in range(nG): 
    #             for j2 in Gnodes[k].adjS:
    #                 if j2==j:
    #                     NG_adj.append(k);
                        
    #         Model.addConstr(GV.Sliq[j,tau]==quicksum(GV.flowGL[k,j,tau] for k in NG_adj));
    #         Model.addConstr(GV.Svpr[j,tau]==quicksum(GV.flowVG[j,k,tau] for k in NG_adj));
                
    # # C6: Sliq limit
    # Model.addConstrs(GV.Sliq[j,tau]<=Exist_SVL[j].liq_cap for j in range(nSVL) for tau in Tg);
    
    # # C9: Svpr limit
    # Model.addConstrs(GV.Svpr[j,tau]<=Exist_SVL[j].vap_cap+GV.Xvpr[j] for j in range(nSVL) for tau in Tg);
    
    # # C10: Sstr limit
    # Model.addConstrs(GV.Sstr[j,tau]<=Exist_SVL[j].str_cap+GV.Xstr[j] for j in range(nSVL) for tau in Tg);


    # no new pipeline 
    #Model.addConstrs(GV.Zg[b]==0 for b in range(nPipe));

def Coupling_constraints(Model,Data):
    # fetch data     
    Enodes=  Data.Enodes;
    Gnodes = Data.Gnodes;
    Branches = Data.Branches;
    node2zone= Data.node2zone;
    PipeLines = Data.PipeLines;
    Plants = Data.Plants;
    eStore = Data.eStore;
    Other_input = Data.Other_input;
    state2zone_id = Data.state2zone_id;
    plant2sym = Data.plant2sym;
    sym2plant = Data.sym2plant;
    zone_id2state = Data.zone_id2state;
    Exist_SVL = Data.Exist_SVL;
    SVLs = Data.SVLs;
    type_prod_lim = Data.type_prod_lim;
    nE = len(Enodes); 
    nG = len(Gnodes);
    nBr = len(Branches);
    nPipe = len(PipeLines);
    nPlt = len(Plants);
    neSt = len(eStore);
    Te = range(len(e_time_weight));
    Tg = range(len(g_time_weight));
    nPipe = len(PipeLines);
    nSVL = len(Exist_SVL);
    nG = len(Gnodes);
    
    
    NG_units = ["ng","CT","CC","CC-CCS"];
    if (Setting.emis_case !=1): #power system only
        #Model.addConstrs(GV.flowGE[k,n,tau]==quicksum(Plants[i].heat_rate*EV.prod[n,t,i] for t in range(tau*24,(tau+1)*24) for i in range(nPlt) if plant2sym[i] in NG_units) for k in range(nG) for tau in Tg for n in Gnodes[k].adjE);
        Model.addConstrs(quicksum(GV.flowGE[k,n,tau] for k in range(nG) if n in Gnodes[k].adjE)==quicksum(Plants[i].heat_rate*EV.prod[n,t,i] for t in range(tau*24,(tau+1)*24) for i in range(nPlt) if plant2sym[i] in NG_units) for n in range(nE) for tau in Tg);

    
    ex_xi = LinExpr(quicksum(g_time_weight[tau]*GV.flowGE[k,n,tau] for k in range(nG) for tau in Tg for n in Gnodes[k].adjE));
    e_emis = LinExpr(quicksum(e_time_weight[t]*Plants[i].emission*Plants[i].heat_rate*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if plant2sym[i] in NG_units));
    #g_emis = LinExpr(quicksum(g_time_weight[tau]*Other_input.NG_emission*GV.supply[k,tau] for k in range(nG) for tau in Tg))-Other_input.NG_emission*ex_xi;
    
    g_emis = LinExpr(quicksum(g_time_weight[tau]*Other_input.NG_emission*(Gnodes[k].demand[g_rep_days[tau]]-GV.RNG_use[k,tau]-GV.Shed[k,tau]) for k in range(nG) for tau in Tg));
    Model.addConstr(EV.emis_amount==e_emis);
    Model.addConstr(GV.emis_amount==g_emis);
    
   # if Setting.emis_case==4:  #global emission limit
    Model.addConstr(EV.emis_amount+GV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.CO2_emission_1990);
    # if Setting.emis_case==3:#global emission limit with no RNG
    #     Model.addConstr(EV.emis_amount+GV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.CO2_emission_1990);
    #     Model.addConstrs(GV.RNG_use[k,tau]==0 for k in range(nG) for tau in Tg);

    # if Setting.emis_case==2: # JPoNG with no RNG and emis const on power system    
    #     Model.addConstr(EV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.e_emis_lim);
    #     Model.addConstrs(GV.RNG_use[k,tau]==0 for k in range(nG) for tau in Tg);
    # if Setting.emis_case==1:
    #     Model.addConstr(EV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.e_emis_lim);
    # if Setting.emis_case==0: # no coupling constraint
    #     pass;


def Get_var_vals(Model,Data):
    # fetch data     
    Enodes=  Data.Enodes;
    Gnodes = Data.Gnodes;
    Branches = Data.Branches;
    node2zone= Data.node2zone;
    PipeLines = Data.PipeLines;
    Plants = Data.Plants;
    eStore = Data.eStore;
    Other_input = Data.Other_input;
    state2zone_id = Data.state2zone_id;
    plant2sym = Data.plant2sym;
    sym2plant = Data.sym2plant;
    zone_id2state = Data.zone_id2state;
    Exist_SVL = Data.Exist_SVL;
    SVLs = Data.SVLs;
    type_prod_lim = Data.type_prod_lim;
    nE = len(Enodes); 
    nG = len(Gnodes);
    nBr = len(Branches);
    nPipe = len(PipeLines);
    nPlt = len(Plants);
    neSt = len(eStore);
    Te = range(len(e_time_weight));
    Tg = range(len(g_time_weight));
    nPipe = len(PipeLines);
    nSVL = len(Exist_SVL);
    nG = len(Gnodes);
    
    EV.Xop_val = Model.getAttr('x',EV.Xop);
    EV.Xest_val = Model.getAttr('x',EV.Xest);
    EV.Xdec_val = Model.getAttr('x',EV.Xdec);
    #EV.X_val = Model.getAttr('x',EV.X);
    #EV.Xup_val = Model.getAttr('x',EV.Xup);
    
    EV.prod_val = Model.getAttr('x',EV.prod);
    EV.flowE_val = Model.getAttr('x',EV.flowE);
    EV.YeCD_val = Model.getAttr('x',EV.YeCD);
    EV.YeLev_val = Model.getAttr('x',EV.YeLev);
    EV.eSlev_val = Model.getAttr('x',EV.eSlev);
    EV.Shed_val = Model.getAttr('x',EV.Shed);
    EV.Ze_val = Model.getAttr('x',EV.Ze);
    EV.eSdis_val = Model.getAttr('x',EV.eSdis);
    EV.eSch_val = Model.getAttr('x',EV.eSch);
    EV.kappa_capt_val = Model.getAttr('x',EV.kappa_capt);
    EV.kappa_pipe_val =  Model.getAttr('x',EV.kappa_pipe);    
    EV.est_cost_val = EV.est_cost.X;
    EV.decom_cost_val = EV.decom_cost.X;
    EV.FOM_cost_val = EV.FOM_cost.X;
    EV.VOM_cost_val = EV.VOM_cost.X;
    EV.nuc_fuel_cost_val = EV.nuc_fuel_cost.X;
    EV.startup_cost_val = EV.startup_cost.X;
    EV.shedding_cost_val = EV.shedding_cost.X;
    EV.elec_storage_cost_val = EV.elec_storage_cost.X;
    EV.est_trans_cost_val = EV.est_trans_cost.X;
    EV.emis_amount_val = EV.emis_amount.X;
    EV.gas_fuel_cost_val = EV.gas_fuel_cost.X;
    EV.e_system_cost_val = EV.e_system_cost.X;
    #if Setting.emis_case==1:
    EV.e_system_cost_val = EV.e_system_cost.X;#-EV.gas_fuel_cost.X;
    
    # print(Model.getAttr('x',GV.Sstr));
    # print(Model.getAttr('x',GV.Svpr));
    # print(Model.getAttr('x',GV.flowGG));
    
    GV.Zg_val = Model.getAttr('x',GV.Zg);
    GV.supply_val = Model.getAttr('x',GV.supply);
    #GV.RNG_supply_val = Model.getAttr('x',GV.supply);
    GV.Shed_val = Model.getAttr('x',GV.Shed);
    GV.RNG_use_val = Model.getAttr('x',GV.RNG_use);
    GV.flowGE_val = Model.getAttr('x',GV.flowGE);
    GV.g_system_cost_val = GV.g_system_cost.X;
    
    # note this one
    GV.import_cost_val = GV.import_cost.X;
    if Setting.emis_case==1:
        GV.import_cost_val = GV.import_cost.X-EV.gas_fuel_cost.X;

    
    
    GV.RNG_cost_val = GV.RNG_cost.X;
    GV.fom_str_cost_val = GV.fom_str_cost.X;
    GV.shed_cost_val = GV.shed_cost.X;
    GV.inv_pipe_cost_val = GV.inv_pipe_cost.X;
    GV.emis_amount_val = GV.emis_amount.X;
    GV.inv_str_cost_val = GV.inv_str_cost.X;
    
    
    
#(EV.prod[n,t,i]<=Enodes[n].cap_factors[e_rep_hrs[t],i]* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te);

    # for n in range(nE):
    #     for t in Te:
    #         for i in range(nPlt):
    #             if plant2sym[i] =='solar-UPV':
    #                 print(f"prod[{n},{t},{i}] = {EV.prod_val[n,t,i]}");
    #                 print(f"cap: {Enodes[n].cap_factors[e_rep_hrs[t],i]}, pmax = {Plants[i].nameplate_cap}, xop:{EV.Xop_val[n,i]}")


def Publish_results(s_time,MIP_gap,Data):
    # fetch data     
    Enodes=  Data.Enodes;
    Gnodes = Data.Gnodes;
    Branches = Data.Branches;
    node2zone= Data.node2zone;
    PipeLines = Data.PipeLines;
    Plants = Data.Plants;
    eStore = Data.eStore;
    Other_input = Data.Other_input;
    state2zone_id = Data.state2zone_id;
    plant2sym = Data.plant2sym;
    sym2plant = Data.sym2plant;
    zone_id2state = Data.zone_id2state;
    Exist_SVL = Data.Exist_SVL;
    SVLs = Data.SVLs;
    type_prod_lim = Data.type_prod_lim;
    nE = len(Enodes); 
    nG = len(Gnodes);
    nBr = len(Branches);
    nPipe = len(PipeLines);
    nPlt = len(Plants);
    neSt = len(eStore);
    Te = range(len(e_time_weight));
    Tg = range(len(g_time_weight));
    nPipe = len(PipeLines);
    nSVL = len(Exist_SVL);
    nG = len(Gnodes);
    
    num_e_str = 0; str_lev = 0;str_cap=0;num_ze = 0;total_shed=0;total_flow=0;
    [num_e_str := num_e_str+1 if EV.YeLev_val[n,r]>0 else num_e_str+0 for r in range(neSt) for n in range(nE)];
    [str_lev := str_lev+EV.YeLev_val[n,r]for r in range(neSt) for n in range(nE)];
    [str_cap := str_cap +EV.YeCD_val[n,r] for r in range(neSt) for n in range(nE)];
    [num_ze := num_ze+1 if EV.Ze_val[b]>0 else num_ze+0 for b in range(nBr)];
    [total_shed := total_shed+EV.Shed_val[n,t] for n in range(nE) for t in Te];
    [total_flow:= total_flow+abs(EV.flowE_val[b,t]) for b in range(nBr) for t in Te];
    total_flow = total_flow/2;
    
    # co2_capt_total = 0;co2_pipe_cap = 0
    # [co2_capt_total:= co2_capt_total+e_time_weight[t]*EV.kappa_capt_val[n,t] for n in range(nE) for t in Te];
    # [co2_pipe_cap:= co2_pipe_cap+EV.kappa_pipe_val[n] for n in range(nE)];
    # print(co2_capt_total);
    # print(co2_pipe_cap);
    
    pr = np.zeros(nPlt);
    est = np.zeros(nPlt);
    dec = np.zeros(nPlt);
    
    for i in range(nPlt):
        s1 = 0;
        [s1:= s1+EV.prod_val[n,t,i]*e_time_weight[t] for n in range(nE) for t in Te];
        pr[i] = s1;s1 = 0;
        [s1:= s1+EV.Xest_val[n,i] for n in range(nE)];
        est[i] = s1;s1 = 0;
        [s1:= s1+EV.Xdec_val[n,i] for n in range(nE)];
        dec[i] = s1;
    
    num_est_pipe = 0;total_ng_shed = 0;total_rng=0;total_fge=0;
    [num_est_pipe := num_est_pipe+1 if GV.Zg_val[p]>0 else num_est_pipe+0 for p in range(nPipe)];
    [total_ng_shed := total_ng_shed+GV.Shed_val[n,tau] for n in range(nG) for tau in Tg];
    [total_rng := total_rng+GV.RNG_use_val[n,tau] for n in range(nG) for tau in Tg];
    #total_rng = total_rng/Setting.poss_gas_consump;
    [total_fge := total_fge+GV.flowGE_val[j,n,tau] for j in range(nG) for n in range(nE) for tau in Tg];
    elapsed = time.time()-s_time;
    header0 = ['Power_network_size','cluster_method','Rep-Days','Emis-case','Elec_scenario', 'reduc-goal','RPS','UC-active?',
              'UC-rlx?','int-vars-rlx?','MI-gap(%)', 'Run time(sec)','Total-cost','',
              'Power-cost','est-cost','decom-cost','FOM','VOM','nuc-fuel-cost','gas-fuel-cost','startup-cost','shed-cost',
              'storage-cost','tran-est-cost','CCS-cost','emission_e','num-est-tran',
              'total-str-lev','total-str-cap','num-est-str','total-flow','',
              'NG-cost','NG-import-cost','RNG-import-cost','inv-storage','FOM-storage','shed-cost','pipe-est-cost',
              'emission_g','num-est-pipe','total-ng-shed','total-RNG-import','total-flowGE','','Production:',
              'ng','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new','','established:', # production
              'ng','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new','','decommissioned:', # established
              'ng','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new','', # decommissioned
              ];
    row = []; row.append(Setting.Power_network_size); row.append(Setting.rep_day_folder);
    row.append(Setting.num_rep_days);row.append(Setting.emis_case);
    row.append(Setting.electrification_scenario);
    row.append(Setting.emis_reduc_goal);
    row.append(Setting.VRE_share);row.append(Setting.UC_active);
    row.append(Setting.relax_UC_vars);row.append(Setting.relax_int_vars);row.append(MIP_gap);
    row.append(elapsed);
    row.append(EV.e_system_cost_val+GV.g_system_cost_val);
    row.append(''); 
    row.append(EV.e_system_cost_val); row.append(EV.est_cost_val);
    row.append(EV.decom_cost_val); row.append(EV.FOM_cost_val);row.append(EV.VOM_cost_val);
    row.append(EV.nuc_fuel_cost_val);row.append(EV.gas_fuel_cost_val);
    row.append(EV.startup_cost_val);
    row.append(EV.shedding_cost_val);
    row.append(EV.elec_storage_cost_val);row.append(EV.est_trans_cost_val);
    row.append(0);
    row.append(EV.emis_amount_val);row.append(num_ze);row.append(str_lev);row.append(str_cap);
    row.append(num_e_str); row.append(total_flow);
    row.append('');
    row.append(GV.g_system_cost_val);row.append(GV.import_cost_val);row.append(GV.RNG_cost_val);
    row.append(GV.inv_str_cost_val);
    row.append(GV.fom_str_cost_val);row.append(GV.shed_cost_val);row.append(GV.inv_pipe_cost_val);
    row.append(GV.emis_amount_val);row.append(num_est_pipe);row.append(total_ng_shed);
    row.append(total_rng);row.append(total_fge);row.append('');row.append('Production:');
    
    for i in range(nPlt):
        row.append(pr[i]);
    row.append('');row.append('established:');
    for i in range(nPlt):
        row.append(est[i]);
    row.append('');row.append('decommissioned:');
    for i in range(nPlt):
        row.append(dec[i]);
    
    
    with open(os.getcwd()+'/JPoNG_Results.csv','a',encoding='UTF8',newline='') as f:
        writer = csv.writer(f);
        if Setting.print_result_header:
            writer.writerow(header0);
        writer.writerow(row);
        f.close();
        
    name =os.getcwd()+'/'+ str(Setting.Power_network_size)+'-'+ str(Setting.num_rep_days)+'-'+Setting.electrification_scenario+'-'+str(Setting.emis_case)+'-'+str(Setting.emis_reduc_goal)+'-'+str(Setting.VRE_share)+'.csv';

    
    if Setting.print_all_vars:
        header = ['ng','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new'
                   ,'Charge','Discharge',
                   'zg','from-g','to-g','dist-g','ze','from-e','to-e','dist-e',
                   'ng-supply(1-18)-(day by node)','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18',
                   'RNG-supply(1-18)-(day by node)','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18',
                   'ng-prod','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new',
                   'ng-est','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new',
                   'ng-dec','solar','wind','hydro','nuclear','CT','CC','CC-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new',
                   'co2_pipe_cap-1','co2_pipe_cap-2','co2_pipe_cap-3','co2_pipe_cap-4','co2_pipe_cap-5','co2_pipe_cap-6',
                   'co2-capt-1','co2-capt-2','co2-capt-3','co2-capt-4','co2-capt-5','co2-capt-6'];
        #'ng-marginal-price(1-18)-(day by node)','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'
       
        row0 = [['']*140 for i in range(max(368,Setting.num_rep_days*24+3))];
        #row0 = np.zeros((max(96,Setting.num_rep_days*24+3),77))-1;
        #row0[:] = np.nan;
        col = 0; 
        #row0 = [[] for t in range(len(Te))];
        for t in Te:
            for i in range(nPlt):
                r0 = [EV.prod_val[n,t,i]  for n in range(nE)];
                #row0[t].append(np.round(sum(r0)));
                row0[t][i+col] = np.round(sum(r0));
            r0 = [EV.eSch_val[n,t,0] for n in range(nE)];
            #row0[t].append(np.round(sum(r0)));
            row0[t][nPlt] = np.round(sum(r0));
            r0 = [EV.eSdis_val[n,t,0] for n in range(nE)];
            #row0[t].append(np.round(sum(r0)));
            row0[t][nPlt+1] = np.round(sum(r0));
            
        col = nPlt+2;#charge and discharge
        for p in range(nPipe):
            # row0[p].append(GV.Zg_val[p]);
            # row0[p].append(PipeLines[p].from_node);
            # row0[p].append(PipeLines[p].to_node);
            # row0[p].append(PipeLines[p].length);
            row0[p][col] = GV.Zg_val[p];
            row0[p][col+1] = PipeLines[p].from_node;
            row0[p][col+2] = PipeLines[p].to_node;
            row0[p][col+3] = PipeLines[p].length;
            
        col += 4;
        for br in range(nBr):
            # row0[br].append(EV.Ze_val[br]);
            # row0[br].append(Branches[br].from_node);
            # row0[br].append(Branches[br].to_node);
            # row0[br].append(Branches[br].length);
            row0[br][col] = EV.Ze_val[br];
            row0[br][col+1] = Branches[br].from_node;
            row0[br][col+2] = Branches[br].to_node;
            row0[br][col+3] = Branches[br].length;
        col += 4;
        for tau in Tg:
            for k in range(nG):
                #row0[tau].append(GV.supply_val[k,tau]);
                row0[tau][k+col] = GV.supply_val[k,tau];
            
        col += nG;
        for tau in Tg:
            for k in range(nG):
                #row0[tau].append(GV.supply_val[k,tau]);
                row0[tau][k+col] = GV.RNG_use_val[k,tau];
            
        col += nG;       
        
        
        for n in range(nE):
            for i in range(nPlt):
                r0 = [EV.prod_val[n,t,i]*e_time_weight[t] for t in Te];
                #row0[n].append(np.round(sum(r0)));
                row0[n][i+col] = np.round(sum(r0));
        col += nPlt;    
        for n in range(nE):
            for i in range(nPlt):                
                #row0[n].append(EV.Xest_val[n,i]);
                row0[n][i+col] = EV.Xest_val[n,i];
        col += nPlt;
        for n in range(nE):
            for i in range(nPlt):                
                #row0[n].append(EV.Xdec_val[n,i]);
                row0[n][i+col] = EV.Xdec_val[n,i];
        
        col += nPlt;
       
        # for n in range(nE):
        #     row0[0][col+n]= np.round(EV.kappa_pipe_val[n],1);

        # col += nE;
        # for n in range(nE):
        #     for t in Te:
        #         row0[t][col+n] = np.round(EV.kappa_capt_val[n,t],1);
        
        # # col += nE;
        # for k in range(nG):
        #     for t in Tg:
        #         row0[t][col+k] = GV.marginal_prices_val[k,t];
        with open(name,'a',encoding='UTF8',newline='') as fid:
            writer=csv.writer(fid);
            writer.writerow(header0);
            writer.writerow(row);
            writer.writerow(header);
            for r in range(len(row0)):
                writer.writerow(row0[r]);
            fid.close();


   
def JPoNG_full_year_LP(Model,Data,days,hours): # used in the AAAI paper
# fetch data     

    Enodes=  Data.Enodes;
    Gnodes = Data.Gnodes;
    Branches = Data.Branches;
    node2zone= Data.node2zone;
    PipeLines = Data.PipeLines;
    Plants = Data.Plants;
    eStore = Data.eStore;
    Other_input = Data.Other_input;
    state2zone_id = Data.state2zone_id;
    plant2sym = Data.plant2sym;
    sym2plant = Data.sym2plant;
    zone_id2state = Data.zone_id2state;
    Exist_SVL = Data.Exist_SVL;
    SVLs = Data.SVLs;
    type_prod_lim = Data.type_prod_lim;
    nE = len(Enodes); 
    nG = len(Gnodes);
    nBr = len(Branches);
    nPipe = len(PipeLines);
    nPlt = len(Plants);
    neSt = len(eStore);
    Te = range(len(hours));
    Tg = range(len(days));
    nPipe = len(PipeLines);
    nSVL = len(Exist_SVL);
    nG = len(Gnodes);
     #% define decision variables 
    EV.prod= Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    EV.theta = Model.addVars(nE,len(Te),lb=np.zeros((nE,len(Te)))-GRB.INFINITY,ub=np.zeros((nE,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    EV.Shed = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    
    EV.YeCD = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    EV.YeLev = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    #EV.YeStr = Model.addVars(nE,neSt,vtype=GRB.BINARY);
    
    
    EV.eSch =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSdis =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSlev =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    
    EV.flowE =  Model.addVars(nBr,len(Te),lb=np.zeros((nBr,len(Te)))-GRB.INFINITY,ub=np.zeros((nBr,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    
    
    EV.est_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.est_trans_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.decom_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.FOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.startup_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.VOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.nuc_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.gas_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.shedding_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.elec_storage_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.e_system_cost=Model.addVar(vtype=GRB.CONTINUOUS);


    #% Electricity System Objective Function    
    est_cost = 0;
    dec_cost = 0;
    fom_cost = 0;
    vom_cost = LinExpr(quicksum(Plants[i].VOM*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    nuc_fuel_cost = LinExpr(quicksum(Plants[i].heat_rate *Other_input.Nuclear_price*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if (Plants[i].Type=="nuclear" or Plants[i].Type=="nuclear-new") ));
    
    startup_cost = 0;
    shed_cost = LinExpr(quicksum(Setting.e_shed_penalty*EV.Shed[n,t] for n in range(nE) for t in Te));
    strg_cost = LinExpr(quicksum(eStore[r].est_coef*(eStore[r].power_capex*EV.YeCD[n,r]+eStore[r].energy_capex*EV.YeLev[n,r]) for n in range(nE) for r in range(neSt)));
    strg_cost += LinExpr(quicksum((eStore[r].pFOM*EV.YeCD[n,r]+eStore[r].eFOM*EV.YeLev[n,r]) for n in range(nE) for r in range(neSt)));
    tran_cost = 0;
    
    # total power system cost function
    e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+startup_cost+shed_cost+strg_cost+tran_cost;
    #e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+shed_cost;
    
    Model.addConstr(EV.est_cost == est_cost);
    Model.addConstr(EV.decom_cost == dec_cost);
    Model.addConstr(EV.FOM_cost == fom_cost);
    Model.addConstr(EV.VOM_cost == vom_cost);
    Model.addConstr(EV.nuc_fuel_cost == nuc_fuel_cost);
    Model.addConstr(EV.startup_cost == startup_cost);
    Model.addConstr(EV.shedding_cost == shed_cost);
    Model.addConstr(EV.elec_storage_cost == strg_cost);
    Model.addConstr(EV.est_trans_cost==tran_cost);
    Model.addConstr(EV.e_system_cost==e_total_cost);            
    #% Electricity System Constraints
    # C1: number of generation units at each node
         
    # C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)    
    Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap*EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
  
    Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
    Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
    
    # C5, C6: flow limit for electricity
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1))
        
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow*EV.Ze_val[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow*EV.Ze_val[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0))
    
    # C7: power balance  
    Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt))+quicksum(Enodes[n].arc_sign[b]*EV.flowE[Enodes[n].arcs[b],t] for b in range(len(Enodes[n].arcs)))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt))+EV.Shed[n,t] == Enodes[n].demand[hours[t]] for n in range(nE) for t in Te);
    
    # # C8: flow equation
    # Model.addConstrs(EV.flowE[b,t]==Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t]) for b in range(nBr) for t in Te if Branches[b].is_exist==1);
    # Model.addConstrs(EV.flowE[b,t]-Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1-EV.Ze_val[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
    # Model.addConstrs(-EV.flowE[b,t]+Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1-EV.Ze_val[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
        
    # # # C9: phase angle (theta) limits. already applied in the definition of the variable
    # Model.addConstrs(EV.theta[n,t]<=Other_input.pi for n in range(nE) for t in Te);
    # Model.addConstrs(-EV.theta[n,t]<=Other_input.pi for n in range(nE) for t in Te);
    # Model.addConstrs(EV.theta[0,t]==0 for t in Te);
    
    # C10: VRE production according to capacity factors
    Model.addConstrs(EV.prod[n,t,i]<=Enodes[n].cap_factors[hours[t],i]* Plants[i].nameplate_cap*EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te);
    
    # C11: demand curtainlment constraint
    Model.addConstrs(EV.Shed[n,t]<= Enodes[n].demand[hours[t]] for  n in range(nE) for t in Te);
    
    # C12: RPS constraints
    VRE = ["solar","wind","wind-offshore-new","solar-UPV","wind-new"];# hydro not included
    Model.addConstr(quicksum(EV.prod[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if(plant2sym[i] in VRE)) >= Setting.VRE_share*quicksum(Enodes[n].demand[hours[t]] for n in range(nE) for t in Te));
    
    # C14,C15,C16 storage constraints
    Model.addConstrs(EV.eSlev[n,0,r]==eStore[r].eff_ch*EV.eSch[n,0,r]-EV.eSdis[n,0,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt));
    Model.addConstrs(EV.eSlev[n,t,r]-EV.eSlev[n,t-1,r]==eStore[r].eff_ch*EV.eSch[n,t,r]-EV.eSdis[n,t,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt) for t in Te if t>0);
    
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSdis[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSch[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeLev[n,r]>= EV.eSlev[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    
    # start and ending of storage should be the same in case of using rep. days
    if len(days)<364:
        for k in range(len(days)):
            s1 = k*24; s2=(k+1)*24-1;
            Model.addConstrs(EV.eSlev[n,s1,r]==EV.eSlev[n,s2,r] for n in range(nE) for r in range(neSt));

 


 # NG system
    GV.Xvpr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Xstr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Sstr = Model.addVars(nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.Svpr = Model.addVars(nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.Sliq = Model.addVars(nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.supply = Model.addVars(nG,len(Tg),vtype = GRB.CONTINUOUS);
    GV.Shed =  Model.addVars(nG,len(Tg),vtype = GRB.CONTINUOUS);
    GV.RNG_use = Model.addVars(nG,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowGG =  Model.addVars(nPipe,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowGE =  Model.addVars(nG,nE,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowGL =  Model.addVars(nG,nSVL,len(Tg),vtype = GRB.CONTINUOUS);
    GV.flowVG =  Model.addVars(nSVL,nG,len(Tg),vtype = GRB.CONTINUOUS);
    
    GV.inv_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.inv_pipe_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.shed_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.RNG_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.fom_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.import_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.g_system_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    
    # NG System Objective Function
    s1 = 0;
    s2 = LinExpr(quicksum(GV.supply[k,tau]*Other_input.NG_price for k in range(nG) for tau in Tg));
    s3 = LinExpr(quicksum(SVLs[0].inv_coef*(SVLs[0].capex*GV.Xstr[j]+SVLs[1].capex*GV.Xvpr[j]) for j in range(nSVL)));
    s4 = LinExpr(quicksum(SVLs[1].FOM*(Exist_SVL[j].vap_cap+GV.Xvpr[j])+SVLs[0].FOM*(Exist_SVL[j].str_cap+GV.Xstr[j]) for j in range(nSVL)));
    s5 = LinExpr(quicksum(Setting.g_shed_penalty*GV.Shed[k,tau] for k in range(nG) for tau in Tg));
    s6 = LinExpr(quicksum(Other_input.RNG_price*GV.RNG_use[k,tau]for k in range(nG) for tau in Tg));
    
    
    s7 = s1+s2+s3+s4+s5+s6;
    Model.addConstr(GV.inv_pipe_cost==s1);
    Model.addConstr(GV.import_cost==s2);
    Model.addConstr(GV.inv_str_cost==s3);
    Model.addConstr(GV.fom_str_cost==s4);
    Model.addConstr(GV.shed_cost==s5);
    Model.addConstr(GV.RNG_cost==s6);
    Model.addConstr(GV.g_system_cost==s7);
    
    # NG System Constraints
    #C1, C2: flow limit for NG
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap for i in range(nPipe) for tau in Tg if PipeLines[i].is_exist==1);
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap*GV.Zg_val[i] for i in range(nPipe) for tau in Tg if PipeLines[i].is_exist==0);
    
    # C3: flow balance, NG node
    Model.addConstrs(GV.supply[k,tau]
                      -quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_exp)
                      +quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_imp)
                      -quicksum(GV.flowGE[k,n,tau] for n in Gnodes[k].adjE)                      
                      +GV.Shed[k,tau]+GV.RNG_use[k,tau]
                     == Gnodes[k].demand[days[tau]] for k in range(nG) for tau in Tg);                             
    
    # C3,C4: injection (supply) limit and curtailment limit
    Model.addConstrs(GV.supply[k,tau]<= Gnodes[k].injU for k in range(nG) for tau in Tg);
    Model.addConstrs(GV.Shed[k,tau]+GV.RNG_use[k,tau]<= Gnodes[k].demand[days[tau]] for k in range(nG) for tau in Tg);
    
    # # C5: storage balance
    # Model.addConstrs(GV.Sstr[j,tau]==Exist_SVL[j].str_cap*0+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in Tg if tau==0);
    # Model.addConstrs(GV.Sstr[j,tau]==(1-SVLs[0].BOG)*GV.Sstr[j,tau-1]+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in Tg if tau>0);
    
    # # C6,8: calculate Sliq, Svpr
    # for j in range(nSVL):
    #     for tau in Tg:
    #         NG_adj = [];
    #         for k in range(nG): 
    #             for j2 in Gnodes[k].adjS:
    #                 if j2==j:
    #                     NG_adj.append(k);
                        
    #         Model.addConstr(GV.Sliq[j,tau]==quicksum(GV.flowGL[k,j,tau] for k in NG_adj));
    #         Model.addConstr(GV.Svpr[j,tau]==quicksum(GV.flowVG[j,k,tau] for k in NG_adj));
                
    # # C6: Sliq limit
    # Model.addConstrs(GV.Sliq[j,tau]<=Exist_SVL[j].liq_cap for j in range(nSVL) for tau in Tg);
    
    # # C9: Svpr limit
    # Model.addConstrs(GV.Svpr[j,tau]<=Exist_SVL[j].vap_cap+GV.Xvpr[j] for j in range(nSVL) for tau in Tg);
    
    # # C10: Sstr limit
    # Model.addConstrs(GV.Sstr[j,tau]<=Exist_SVL[j].str_cap+GV.Xstr[j] for j in range(nSVL) for tau in Tg);

 # coupling constraints
    Coupling_constraints(Model,Data);   
    
   
    
def enforce_inv_decisions_cluster(Original_network_size,sp, Model,nPlt):    
    dfs = pd.read_csv(str(Original_network_size)+'_cluster_assignments.csv');
    dff = pd.read_csv(str(sp)+'_cluster_assignments.csv'); # full network (88)
    for c in range(Original_network_size):
        s1 = dfs[dfs['Cluster']==c];
        s2 = np.array(s1['Node']);
        s1 = dff[dff['Node'].isin(s2)];
        s2 = np.array(s1['Cluster'].unique());
        Model.addConstrs(quicksum(EV.Xest[n,i] for n in s2) >= (EV.Xest_val[c,i]-1) for  i in range(nPlt));    
        Model.addConstrs(quicksum(EV.Xdec[n,i] for n in s2) <= (EV.Xdec_val[c,i]+1) for  i in range(nPlt));
        
    #Model.addConstr(quicksum(EV.Xest[n,i] for n in s2 for i in range(nPlt)) <= quicksum(int(EV.Xest_val[c,i]) for  i in range(nPlt)));    
    #Model.addConstr(quicksum(EV.Xdec[n,i] for n in s2 for i in range(nPlt)) <= quicksum(int(EV.Xdec_val[c,i]) for  i in range(nPlt)));    

    