def path2label(path_in, label_out):
    """
    This function converts a freesurfer path file into a freesurfer label file. Separate paths 
    within the same path file are put into one label file with different label numbers.
    Inputs:
        *path_in: filename of input path file.
        *label_out: filename of output label file.
        
    created by Daniel Haenelt
    Date created: 04-03-2020          
    Last modified: 04-03-2020
    """
    import numpy as np

    # read path file
    fileID = open(path_in, "r")

    # skip path header and get number of points for each path
    data = []
    num_points = []
    access = 0
    line = fileID.readline()
    while line:
        line = fileID.readline()
        if "ENDPATH" in line:
            access = 0
        elif access:
            data = np.append(data, line.split()[-1])
        elif "NUMPOINTS" in line:
            num_points = np.append(num_points, line.split()[-1])
            access = 1
    
    fileID.close()

    # convert to integer
    data = data.astype(int)
    num_points = num_points.astype(int)

    # write label file    
    fileID = open(label_out, "w")    
    fileID.write("#!ascii label , from subject daniel vox2ras=TkReg coords=white\n")
    fileID.write(str(len(data))+"\n")

    num_counter = 0    
    c = 0
    for i in range(len(data)):
        fileID.write(str(data[i])+" 0.000 0.000 0.000 "+str(num_counter+1)+"\n")
        
        c += 1
        if c == num_points[num_counter]:
            c = 0
            num_counter += 1
    
    fileID.close()