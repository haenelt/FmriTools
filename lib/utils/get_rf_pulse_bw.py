def get_rf_pulse_bw(input, npad=1000, ninterp=1000000, threshold=0.1):
    """
    This function computes the bandwidth in Hz of an RF pulse. The RF pulseshape has to be exported
    from the POET simulation. This can be done by clicking on the wanted pulse in the pulse sequence 
    diagram and saving the event block as textfile. The bandwidth is calculated as the width at the
    threshold value of the normalized frequency magnitude.
    Inputs:
        *input: textfile of one event block.
        *npad: pad zeros before and after pulse.
        *ninterp: number of frequency steps for spline interpolation.
        *threshold: value at which bw is estimated.
    Outputs:
        *rf: pulse shape in time domain (microseconds).
        *rf_fft: pulse shape in frequency domain (cycles/second).
        *bw: bandwidth (Hz).
        
    created by Daniel Haenelt
    Date created: 15-05-2019
    Last modified: 15-05-2019
    """
    import re
    import numpy as np
    import matplotlib.pyplot as plt
    from numpy.fft import fft, fftshift

    # read file
    file=[]
    with open(input) as fp:  
        line = fp.readline()
        while line:
            line = fp.readline()
            file.append(line)

    # remove header from array
    file = file[9:]

    # get rf-pulse shape from array
    rf = []
    for i in range(len(file)):
        temp = [float(s) for s in re.findall(r'-?\d+\.?\d*', file[i])]
        if len(temp) > 0:
            rf.append(temp[0])
    
    # remove zeros at the beginning
    exit_loop = 0
    cnt = 0
    while exit_loop == 0:
        if rf[cnt] == 0:
            cnt += 1
        else:
            exit_loop = 1
    rf = rf[cnt:]

    # remove zeros at the end
    exit_loop = 0
    cnt = -1
    while exit_loop == 0:
        if rf[cnt] == 0:
            cnt -= 1
        else:
            exit_loop = 1
    rf = rf[:cnt+1]

    # pad beginning and ending with zeros
    rf = np.pad(rf,npad,'constant')

    # get magnitude fft
    rf_fft = np.abs(fft(rf))

    # get frequency axis 
    freq = np.fft.fftfreq(len(rf),1e-6) # timesteps in seconds

    # shift center frequency to center and normalize magnitude
    freq = fftshift(freq)
    rf_fft = fftshift(rf_fft)
    rf_fft = rf_fft / np.max(rf_fft)

    # get spline interpolation of fft
    freq_spline = np.linspace(np.min(freq), np.max(freq), ninterp)
    rf_fft_spline = np.interp(freq_spline, freq, rf_fft)

    threshold = 0.1
    for i in range(len(rf_fft_spline)-1):
        if rf_fft_spline[i] <= threshold and rf_fft_spline[i+1] > threshold:
            freq1 = freq_spline[i]
    
        if rf_fft_spline[i] > threshold and rf_fft_spline[i+1] <= threshold:
            freq2 = freq_spline[i]

    # bandwidth in Hz
    bw = freq2 - freq1
    print("Pulse bandwidth in Hz: "+str(bw))

    # plot pulse and corresponding fft
    fig, ax = plt.subplots()
    ax.plot(np.linspace(0,len(rf),len(rf)),rf)
    ax.set_xlabel("Time in microseconds")
    ax.set_ylabel("Amplitude in a.u.")
    ax.set_title("Pulse shape in time domain")
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(freq_spline,rf_fft_spline)
    ax.set_xlabel("frequency in cycles/second")
    ax.set_ylabel("Normalized amplitude in a.u.")
    ax.set_title("Pulse shape in frequency domain")
    plt.show()
    
    return rf, rf_fft, bw