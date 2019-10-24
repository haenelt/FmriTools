"""
Analysis of behavioural data

Analysis of the behavioural performance of the dummy fixation task during scanning. The script 
essentially performs the same calculation as directly done automatically after single functional 
runs.

created by Daniel Haenelt
Date created: 01-10-2019             
Last modified: 01-10-2019  
"""
import os
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib import rc

input = [
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_1/logfiles/p4_GE_EPI3_Run1_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_2/logfiles/p4_GE_EPI3_Run2_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_3/logfiles/p4_GE_EPI3_Run3_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_4/logfiles/p4_GE_EPI3_Run4_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_5/logfiles/p4_GE_EPI3_Run5_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_6/logfiles/p4_GE_EPI3_Run6_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_7/logfiles/p4_GE_EPI3_Run7_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_8/logfiles/p4_GE_EPI3_Run8_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_9/logfiles/p4_GE_EPI3_Run9_odc.mat",
        "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI3/Run_10/logfiles/p4_GE_EPI3_Run10_odc.mat",
         ]

""" do not edit below """

# font parameters for plots
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

for i in range(len(input)):
    
    # get output path and filename
    path_output = os.path.dirname(input[i])
    
    # get mat-file
    FixationData = loadmat(input[i])["FixationData"]

    # number of responses
    response = 0
    change = 0
    for j in range(len(FixationData)):
        if FixationData[j,0] == 3:
            response += 1
        else:
            change += 1

    # number of hits, misses and rt
    change_miss = 0
    change_hit = 0
    rt = []
    for j in range(len(FixationData)-1):
        if FixationData[j,0] != 3 and FixationData[j+1,0] != 3:
            change_miss += 1
        elif FixationData[j,0] != 3 and FixationData[j+1,0] == 3:
            change_hit += 1
            rt.append((FixationData[j+1,1]- FixationData[j,1])*1000)

    # add a miss if the run does not end with a response
    if FixationData[-1,0] != 3:
        change_miss += 1

    # compute error rate
    error_rate = change_miss / change * 100

    # reaction time mean and std
    mean_rt = np.mean(rt)
    std_rt = np.std(rt)

    # output
    fileID = open(os.path.join(path_output,"fixation_task_summary.txt"),"w")
    fileID.write("Number of changes: %i\n" % change)
    fileID.write("Number of responses: %i\n" % response)
    fileID.write("Number of hits: %i\n" % change_hit)
    fileID.write("Number of misses: %i\n" % change_miss)
    fileID.write("Error rate: %.2f %%\n" % error_rate)
    fileID.write("Mean RT: %.2f ms\n" % mean_rt)
    fileID.write("Corresponding SD: %.2f ms" % std_rt)
    fileID.close()

    # hist plot
    fig, ax = plt.subplots()
    ax.hist(rt)
    ax.set_xlabel("RT in ms")
    ax.set_ylabel("Number of responses")
    ax.set_title("Dummy fixation task reaction times")
    fig.savefig(os.path.join(path_output,"fixation_task_hist.png"), format='png', bbox_inches='tight')
    #plt.show()