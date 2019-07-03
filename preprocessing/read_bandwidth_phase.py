"""
Read bandwidth per pixel in phase encoding direction

This scripts reads the bandwidth per pixel phase encode from siemens dicom data. This parameter is 
needed for fieldmap undistortion.

created by Daniel Haenelt
Date created: 03-07-2019             
Last modified: 03-07-2019  
"""
import os
import pydicom

# input data
input = "/home/daniel/Schreibtisch/epi-example-dicom/dicom"

# output data
path_output = "/home/daniel/Schreibtisch" 
name_output = "bla"

""" do not edit below """

with open(os.path.join(path_output,name_output+".txt"), 'w') as f:
    for file in os.listdir(input):
        if file.endswith(".dcm"):
            f.write(str(pydicom.dcmread(file).SeriesNumber)+"\t" \
                    +str(pydicom.dcmread(file).SeriesDescription)+"\t" \
                    +str(pydicom.dcmread(file)["0019","1028"]) \
                    +"\n")