import multiprocessing
import sys
from subprocess import PIPE, Popen
import  os
import subprocess
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



if __name__ == "__main__":
    import time;
    from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor;
    from datetime import datetime;
    
    folder = ['NG=1_CF=1','NG=0_CF=0'];
    net_size = [6];
    rep_days = np.array([5,10,20]);
    cluster_method = ['GNN','GNN-0.5-0.5'];
    # rep_days = [10];
    # cluster_method = ['GNN-0-0'];

    #rep_days = np.insert(rep_days,len(rep_days),365);
    case = [4];
    elec_scen = ['RM'];
    emis_reduc_goal = [0.8];
    solver_gap = 0.01;
    wall_clock_time_lim = 4; #hours;
    solver_thread_num = 4;
    param_list=[];
    SuperClound_Thread = 4;
    for i1 in net_size:
        for i2 in rep_days:
            for i3 in cluster_method:
                for i4 in folder:                    
                    # if i2==10 and i3=='GNN': continue;
                    param = 'python UB.py '+' '+str(i1)+' '+ str(i2)+' '+str(i3)+' '+str(solver_gap)+' '+str(wall_clock_time_lim)+' '+str(solver_thread_num)+' '+i4;                        
                    param_list.append(param);     
                            
    del i1,i2,i3;
    #print(param_list)
    for i in range(int(np.ceil(solver_thread_num*len(param_list)/SuperClound_Thread))):
        j = int(SuperClound_Thread/solver_thread_num);
        with ThreadPoolExecutor() as executor:
            results = executor.map(system, param_list[i*j:(i+1)*j]);
    # for param in param_list:
    #      tmp = system(param)
    # print(tmp)

#%%
# import time
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
# from datetime import datetime

# now = datetime.now().time() # time object

# print("now =", now)
# def sleep_secs(seconds):
#   time.sleep(seconds)
#   print(f'{seconds} has been processed')

# secs_list = [2,4, 6, 8, 10, 12];
# with ThreadPoolExecutor() as executor:
#   results = executor.map(sleep_secs, secs_list)
#   print(results)

# now = datetime.now().time() # time object

# print("now =", now)