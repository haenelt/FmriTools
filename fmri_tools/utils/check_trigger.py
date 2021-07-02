# -*- coding: utf-8 -*-

# external inputs
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


def check_trigger(input_biopac):
    """Check trigger.
    
    This function loads the digital input from a saved biopac mat-file and 
    checks the number of sent triggers and the time difference between triggers.    

    Parameters
    ----------
    input_biopac : str
        Biopac *.mat file

    Returns
    -------
    None.
    
    """
    
    # load digital input from biopac mat-file
    data = sio.loadmat(input_biopac)
    trigger = data["data"][:, 3]

    # count triggers and get onset time of triggers
    trigger_count = 0
    trigger_time = []
    for i in range(len(trigger)-1):
        if trigger[i] == 0 and trigger[i+1] != 0:
            trigger_count += 1
            trigger_time.append(i)
    
    # get time difference between following triggers (digital input is saved in
    # ms)
    trigger_time = np.array(trigger_time)    
    trigger_time1 = trigger_time[:-1]
    trigger_time2 = trigger_time[1:]

    trigger_time_delta = (trigger_time2 - trigger_time1) / 1000
    trigger_time_delta = np.round(trigger_time_delta)  # round to the nearest second

    # plot time difference between following triggers
    plt.plot(trigger_time_delta)

    # print out number of triggers
    print("Number of sent triggers: "+str(trigger_count))
