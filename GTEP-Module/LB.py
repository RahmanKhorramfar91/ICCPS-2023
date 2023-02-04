# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 21:22:46 2022

@author: Rahman Khorramfar
"""

import time;
import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
from Setting import Setting,EV,GV;
from ProblemData_UB import Rub_problem_data, time_weights;
import sys; 
import numpy as np;
#%% Set Default Setting for the Porblem 
Setting.Power_network_size = 88; #88 or 6 
import reconfig4segmentation;
reconfig4segmentation.Cluster_reconfig(Setting.Power_network_size);


Setting.rep_day_folder = 'joint_CF_with_extreme_days';
Setting.num_rep_days = 2;
Setting.emis_case = 4;
Setting.electrification_scenario = 'RM';   # currently HM and RM, but can be MM, MS,MR etc.
Setting.emis_reduc_goal = 0.8; # %80
Setting.VRE_share = 0.5;
#Setting.RNG_cap = 0.2;
Setting.solver_gap = 0.01;
Setting.wall_clock_time_lim = 1; #hour
Setting.UC_active = False;
Setting.relax_UC_vars = True;
Setting.relax_int_vars = 0;
Setting.solver_thread_num = 8;
Setting.e_shed_penalty = 1e4;
Setting.g_shed_penalty = 1e4;


if len(sys.argv)>1:
    print(str(sys.argv));
    Setting.rep_day_folder = sys.argv[1];
    Setting.Power_network_size = int(sys.argv[2]);
    Setting.num_rep_days = int(sys.argv[3]);
    Setting.emis_case = int(sys.argv[4]);
    Setting.electrification_scenario = sys.argv[5];
    Setting.emis_reduc_goal = float(sys.argv[6]); # %80
    Setting.VRE_share = float(sys.argv[7]);
    #Setting.RNG_cap = float(sys.argv[6]);
    Setting.solver_gap = float(sys.argv[8]);
    Setting.wall_clock_time_lim = int(sys.argv[9]); # hour
    Setting.UC_active = bool(int(sys.argv[10]));
    Setting.relax_UC_vars = bool(int(sys.argv[11]));
    Setting.relax_int_vars = bool(int(sys.argv[12]));
    Setting.solver_thread_num = int(sys.argv[13]);

Setting.wall_clock_time_lim = Setting.wall_clock_time_lim*3600; # convert to second for Gurobi
Setting.print_result_header = 1;
Setting.copper_plate_approx = 1; 
Setting.print_all_vars = 0;
s_time = time.time();

Original_network_size = Setting.Power_network_size;
print(sys.argv);
# from ProblemData import Enodes,Gnodes,Branches;
# from ProblemData import PipeLines,Plants,eStore,Other_input,state2zone_id,plant2sym;
# from ProblemData import plant2sym,sym2plant,time_weights,Exist_SVL,SVLs;

#%% run for the rep day model 
# step 1: solve the model for rep days 


#%% solve the full network with a few rep days (2) and 
# impose the sum of investment decisions for the nodes of each cluster 
# to be equal to the investent deicision from step 1
Model = gp.Model();

Setting.Power_network_size = 88; #88 or 6 

reconfig4segmentation.Cluster_reconfig(Setting.Power_network_size);

Setting.rep_day_folder = 'joint_CF_with_extreme_days';

Data = Rub_problem_data();
import Modules_UB;

#%% solve the full problem for the full year on a rolling-horizon bases
# slided the full year into 12 parts (one month), and set the final-day storage level to 0. 

part_length = 365;
ndY = 365;
UB_val = 0;
PipeLines = Data.PipeLines
nE = len(Data.Enodes); 
nBr = len(Data.Branches);
nPlt = len(Data.Plants);
nPipe = len(Data.PipeLines);
    


days =  np.arange(365);

hours = np.array(np.arange(8760));

Model = gp.Model();
Modules_UB.JPoNG_full_year_LP(Model,Data,days,hours)

Model.modelSense = GRB.MINIMIZE;
Model.setObjective(UB_val+GV.g_system_cost+ EV.e_system_cost);
#Model.setParam('OutputFlag', 0);
Model.setParam('MIPGap', Setting.solver_gap);
Model.setParam('Timelimit', Setting.wall_clock_time_lim);
Model.optimize();
UB_val = Model.objVal;
    
Setting.num_rep_days = 365;    
EV.prod_val = Model.getAttr('x',EV.prod);
EV.flowE_val = Model.getAttr('x',EV.flowE);
EV.YeCD_val = Model.getAttr('x',EV.YeCD);
EV.YeLev_val = Model.getAttr('x',EV.YeLev);
EV.eSlev_val = Model.getAttr('x',EV.eSlev);
EV.Shed_val = Model.getAttr('x',EV.Shed);
EV.eSdis_val = Model.getAttr('x',EV.eSdis);
EV.eSch_val = Model.getAttr('x',EV.eSch);    
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
GV.shed_cost_val = UB_val;  # note this. NG shedding is usually zero, so set it to UB
GV.inv_pipe_cost_val = GV.inv_pipe_cost.X;
GV.emis_amount_val = GV.emis_amount.X;
GV.inv_str_cost_val = GV.inv_str_cost.X;
Setting.print_all_vars = 1;
Modules_UB.Publish_results(s_time,0,Data);
print(f"lower Bound is: {UB_val}");
#    s1 = LinExpr(quicksum(PipeLines[b].inv_coef*PipeLines[b].length*Other_input.pipe_per_mile*GV.Zg_val[b] for b in range(nPipe)));   
#    tran_cost = LinExpr(quicksum(Branches[b].est_coef*Other_input.trans_unit_cost*Branches[b].maxFlow*Branches[b].length*EV.Ze_val[b] for b in range(nBr)));
#    est_cost = LinExpr(quicksum(Plants[i].capex*EV.Xest_val[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==0));
#    dec_cost = LinExpr(quicksum(Plants[i].decom_cost*EV.Xdec_val[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==1));
#    fom_cost = LinExpr(quicksum(Plants[i].FOM*EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt)));

