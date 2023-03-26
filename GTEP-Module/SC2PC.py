# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:55:42 2022

@author: Rahman Khorramfar
"""
import multiprocessing
import sys
from subprocess import PIPE, Popen
import  os
import subprocess;
import numpy as np;
def system(command):
    process = Popen(
        args=command,
        stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
        bufsize=0,
        universal_newlines=True,
        shell=True
    )

    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        return -1#raise ProcessException(command, exitCode, output)

    return process.communicate()[0]


folder = 'extreme_days';
net_size = 6;
rep_days = [40];
    #rep_days = np.insert(rep_days,len(rep_days),365);
case = [3];
elec_scen = ['RM','HM'];
emis_reduc_goal = [0.8,0.85,0.9,0.95];
VRE_share = [0.5];
solver_gap = 0.01;
wall_clock_time_lim = 10; #hours;
UC_active = 1;
relax_UC_vars = 1;
relax_int_vars = 0;
solver_thread_num = 4;
param_list=[];
SuperClound_Thread = 96;
param_list=[];

for i2 in case:
    for i1 in rep_days:
        for i3 in elec_scen:
            for i4 in emis_reduc_goal:
                for i5 in VRE_share:
                    param = str(net_size)+'-'+ str(i1)+'-'+str(i3)+'-'+str(i2)+'-'+str(i4)+'-'+str(i5)+'.csv';
                    
                    param_list.append(param);    
                    
system('cd '+os.getcwd());
system('scp -r rkhorramfar@txe1-login.mit.edu:/home/gridsan/rkhorramfar/GTEP-Module/JPoNG_Results.csv ./')

for param in param_list:
    system('scp -r rkhorramfar@txe1-login.mit.edu:/home/gridsan/rkhorramfar/GTEP-Module/'+param+' ./');
    
    
