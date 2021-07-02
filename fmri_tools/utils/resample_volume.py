# -*- coding: utf-8 -*-

# python standard library inputs
import os
import datetime
from shutil import copyfile

# external inputs
import numpy as np
from sh import gunzip

# local inputs
from ..io.get_filename import get_filename


def resample_volume(file_in, file_out, dxyz=[0.4, 0.4, 0.4], rmode="Cu"):
    """Resample volume.

    This function resamples a nifti volume using the afni function 3dresample. 
    Before running the function, set the afni environment by calling AFNI in 
    the terminal.    

    Parameters
    ----------
    file_in : str
        Nifti input filename.
    file_out : str
        Nifti output filename.
    dxyz : list, optional
        Array of target resolution in single dimensions. The default is 
        [0.4, 0.4, 0.4].
    rmode : str, optional
        Interpolation methods (Linear, NN, Cu, Bk). The default is "Cu".

    Returns
    -------
    None.
    
    """

    # get path and file extension of input file
    path_in, _, ext_in = get_filename(file_in)

    # make temporary copy of input file
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    file_tmp = os.path.join(path_in, "tmp_" + tmp_string + ext_in)

    if not os.path.exists(file_tmp) and not os.path.exists(file_tmp[:-3]):
        copyfile(file_in, file_tmp)
    else:
        raise FileExistsError("Temporary file already exists!")

    if os.path.splitext(file_tmp)[1] == ".gz":
        gunzip(file_tmp)
        file_tmp = os.path.splitext(file_tmp)[0]

    # resample volume
    os.system("3dresample " +
              "-dxyz " + str(dxyz[0]) + " " + str(dxyz[1]) + " " + str(dxyz[2]) + " " +
              "-rmode " + str(rmode) + " " +
              "-inset " + file_tmp + " " +
              "-prefix " + file_out)

    # remove temporary copy
    os.remove(file_tmp)
