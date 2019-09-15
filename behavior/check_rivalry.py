"""
Check rivalry switch data

From the condition logfile of the binocular rivalry paradigm, the variable SwitchData is read and 
analyzed. The mean and standard deviation of green and red periods are estimated. The relative
amount of piecemeal events is calculated. This is just a quick overview of the condition file.
Onsets are not classified by their neighbours and not double button presses are considered.

created by Daniel Haenelt
Date created: 31-05-2019             
Last modified: 31-05-2019  
"""
import numpy as np
from scipy.io import loadmat

input = "/home/raid2/haenelt/V2STRIPES/p4/rivalry/GE_EPI1/Run_8/logfiles/p4_GE_EPI1_Run8_rivalry_Cond.mat"

""" do not edit below """

# load mat file
data = loadmat(input)

# load SwitchData key
response = data["SwitchData"][:,0]
time = data["SwitchData"][:,1]

# get time of red and green grating perception
red = []
green = []
for i in range(len(response)-2):
    if response[i] == 1 and response[i+1] == 2:
        red.append(time[i+2]-time[i+1])
    elif response[i] == 2 and response[i+1] == 1:
        green.append(time[i+2]-time[i+1])

# get last point
if response[-2] == 1 and response[-1] == 2:
    red.append(data["onsets"][0][0][0][1]-time[-1])
elif response[-2] == 2 and response[-1] == 1:
    green.append(data["onsets"][0][0][0][1]-time[-1])

# mean and std
green_mean = np.mean(green)
green_std = np.std(green)
red_mean = np.mean(red)
red_std = np.std(red)

print("Green: "+str(green))
print("Green (mean): "+str(green_mean))
print("Green (std): "+str(green_std))
print("Red: "+str(red))
print("Red (mean): "+str(red_mean))
print("Red (std): "+str(red_std))

# get relative amount of mixed perceptions
cnt = 0
for i in range(len(response)):
    if response[i] == 3:
        cnt += 1

print("Piecemeal events in percent: "+str(cnt/len(response)*100))