# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import pydicom


"""
Read bandwidth per pixel in phase encoding direction

This scripts reads the bandwidth per pixel phase encode from siemens dicom data. 
This parameter is needed for fieldmap undistortion.

created by Daniel Haenelt
Date created: 03-07-2019             
Last modified: 12-10-2020  
"""

# input data
input = "/data/pt_01880/dcm/disparity2/"

# output data
path_output = "/data/pt_01880" 
name_output = "test"

# do not edit below

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

# read dicom
count = 0
file_length = len(os.listdir(input))
with open(os.path.join(path_output,name_output+".txt"), 'w') as f:
    for file in os.listdir(input):
        data = pydicom.dcmread(os.path.join(input,file))
        print(file+": file "+str(count)+" of "+str(file_length)+" files.")
        count += 1
        if file.endswith(".ima") or file.endswith(".dcm") and ["0019","1028"] in data:
            f.write(str(data.SeriesNumber)+"\t" \
                    +str(data.SeriesDescription)+"\t" \
                    +str(data["0019","1028"]) \
                    +"\n")
