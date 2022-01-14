# -*- coding: utf-8 -*-

# python standard library inputs
import os
import datetime
import subprocess
import functools
import shutil as sh
from shutil import copyfile
from abc import ABC, abstractmethod

# external inputs
import numpy as np

# local inputs
from .mesh import Mesh
from ..io.get_filename import get_filename

__all__ = ['HeatKernel', 'intracortical_smoothing']


class Filter(ABC):
    NSHUFFLE = 10
    FIT_RANGE = np.arange(1, 10)

    def __init__(self, vtx, fac):
        self.vtx = vtx
        self.fac = fac
        self._mesh = Mesh(vtx, fac)

    @abstractmethod
    def apply(self, data, n_iter):
        pass

    @property
    @functools.lru_cache
    def fit_fwhm(self):
        fwhm_mean = []
        fwhm_std = []
        for i in self.FIT_RANGE:
            fwhm = []
            for j in range(self.NSHUFFLE):
                fwhm.append(self._fwhm(i))

            fwhm_mean.append(np.mean(fwhm))
            fwhm_std.append(np.std(fwhm))

        return fwhm_mean, fwhm_std

    def _fwhm(self, n_iter):

        data = np.random.normal(0, 1, len(self.vtx))
        data_filt = self.apply(data, n_iter)
        edges = self._mesh.edges
        edge_length = self._mesh.avg_edge_length

        var_ds = np.var(data_filt[edges[:,0]] - data_filt[edges[:,1]])
        var_s = np.var(data_filt)

        hallo = 1 - var_ds / (2 * var_s)
        fwhm = edge_length * np.sqrt(-2*np.log(2) / np.log(hallo))
        return fwhm








class HeatKernel(Filter):

    def __init__(self, vtx, fac, sigma):
        super().__init__(vtx, fac)
        #self.vtx = vtx
        #self.fac = fac
        self.sigma = sigma
        #self._mesh = Mesh(vtx, fac)

    @property
    @functools.lru_cache
    def kernel(self):
        """Heat kernel smoothing.

        This function performs heat kernel smoothing [1,2,3] on a triangle mesh. The
        code is mainly adapted from the matlab code by Chung et al. [4]. The kernel
        bandwidth corresponds to diffusion time in the heat equation [3]. The FWHM
        follows 4*sqrt(log 2*n_smooth*sigma) with the natural log.

        If you use this code, please reference one of the following papers. The
        details on the mathematical basis of of the algorithm can be found in these
        papers.

        Parameters
        ----------
        vtx : ndarray
            Vertex points of surface mesh.
        data : ndarray
            Array of vertex-wise sampled data points.
        adjm : ndarray
            Adjacency matrix.
        sigma : float
            Kernel bandwidth.
        n_smooth : int
            Number of iterations.

        Returns
        -------
        res : ndarray
            Array of vertex-wise smoothed data points.

        References
        -------
        .. [1] Chung, MK, et al. Cortical thickness analysis in autism via heat
        kernel smoothing. Neuroimage 25(1), 1256--1265 (2005).
        .. [2] Chung, MK, et al. Unified statistical approach to cortical thickness
        analysis. Inf Process Med Imaging 19, 627--638 (2005).
        .. [3] Chung, MK, et al. Encoding cortical surface by spherical harmonics.
        Statistica Sinica 18, 1269--1291 (2008).
        .. [4] http://pages.stat.wisc.edu/~mchung/softwares/hk/hk_smooth.m

        """

        # heat kernel shape
        k_shape = lambda x, s: np.exp(-x / (2 * s**2)) / np.sum(
            np.exp(-x / (2 * s**2)))

        # parameters
        nvtx = len(self.vtx)  # number of vertices

        # heat kernel weight computation
        neighbor = np.zeros((nvtx, self._max_neighbor + 1), dtype=np.int64)  # including the current vertex
        weight = np.zeros((nvtx, self._max_neighbor + 1))  # including the current vertex

        for i in range(nvtx):

            # get vertex neighbors
            nn = self._mesh.neighborhood(i)
            degree = len(nn)

            # get distance to vertex neighbors
            distance = [np.sum((self.vtx[n, :]-self.vtx[i, :])**2) for n in nn]
            distance.insert(0, 0)
            distance = np.array(distance)

            # get heat kernel weighting for each neighbor
            weight[i, :1 + degree] = k_shape(distance, self.sigma)

            # get corresponding neighbor (add 1 because of dummy row)
            neighbor[i, :1 + degree] = np.append([i], nn) + 1

        return weight, neighbor

    @property
    @functools.lru_cache
    def _max_neighbor(self):
        return np.max(self._mesh.n_neighbors).astype(int)  # max number of first order neighbors

    def apply(self, data, n_iter=1):

        weight, neighbor = self.kernel
        data = np.append(1, data)  # add dummy row

        # iterative kernel smoothing
        for i in range(n_iter):

            # add weights
            res = np.zeros_like(data)
            for j in range(self._max_neighbor):
                res[1:] += data[neighbor[:, j]] * weight[:, j]

            # initialize new data array
            data = res.copy()

        # remove dummy row
        res = res[1:]

        return res

    def apply_inverse(self, data):
        return data - self.apply(data)




def intracortical_smoothing(file_surf, file_overlay, file_out, tan_size=0,
                            rad_start=0, rad_size=1, tan_weights="gauss",
                            cleanup=True):
    """Intracortical smoothing.

    This function applies the freesurfer function mris_smooth_intracortical
    (available since FreeSurfer 7) which enables simultaneous smoothing in
    radial and tangential direction of the cortical sheet. [1]

    Parameters
    ----------
    file_surf : list[str]
        List of input surfaces which must be sorted from wm to pial.
    file_overlay : list[str]
        List of input overlays which must be sorted from wm to pial.
    file_out : str
        File name of smoothed output overlay.
    tan_size : float, optional
        Tangential extent of the smoothing kernel which is defined as the order
        of the vertex neighborhood (number of vertices; max = 6). The default
        is 0.
    rad_start : int, optional
        Starting surface mesh of the intracortical smoothing kernel in radial
        direction (max = number of input surfaces). The default is 0 (white
        surface).
    rad_size : int, optional
        Radial extent of the intracortical smoothing kernel (number of adjacent
        meshes; max = number of input surfaces). The default is 1 (no
        smoothing).
    tan_weights : str, optional
        Weighting function for tangential smoothing. "gauss": tan_size = FWHM,
        "distance": tan_size = radius of the neighborhood around the central
        vertex of the smoothing kernel. The default is "gauss".
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Raises
    ------
    ValueError
        If `file_out` has an invalid file extension.
    FileExistsError
        If `path_temp` already exists.

    Returns
    -------
    None.

    References
    -------
    .. [1] Blazejewska, AI, et al., Intracortical smoothing of small-voxel fMRI
    data can provide increased detection power without spatial resolution
    losses compared to conventional large-voxel fMRI data, NeuroImage 189,
    601--614 (2019).

    """

    # check output file name
    path_output, name_output, ext_output = get_filename(file_out)
    if ext_output not in [".mgh", ".mgz"]:
        raise ValueError("Output file name is expected to have the file "
                         "extension mgh or mgz!")

    # make output folder
    create_folder = 0
    if not os.path.exists(path_output):
        create_folder = 1
        os.makedirs(path_output)

    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_" + tmp_string)
    path_surf = os.path.join(path_temp, "surf")
    path_overlay = os.path.join(path_temp, "overlay")

    # make temporary folder
    if not os.path.exists(path_temp):
        os.mkdir(path_temp)
        os.mkdir(path_surf)
        os.mkdir(path_overlay)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # copy input files to temporary folder
    for i, f in enumerate(file_surf):
        name_tmp = "surf_" + str(i).zfill(4)
        copyfile(f, os.path.join(path_surf, name_tmp))

    for i, f in enumerate(file_overlay):
        name_tmp = "overlay_" + str(i).zfill(4) + ext_output
        copyfile(f, os.path.join(path_overlay, name_tmp))

    # unix command
    command = [
        'mris_smooth_intracortical',
        '--surf_dir', path_surf,
        '--surf_name', 'surf_*',
        '--overlay_dir', path_overlay,
        '--overlay_name', 'overlay_*.mgh',
        '--output_dir', path_output,
        '--output_name', name_output + ext_output,
        '--tan-size', str(tan_size),
        '--rad-size', str(rad_size),
        '--rad-start', str(rad_start),
        '--tan-weights', str(tan_weights)]

    command = ' '.join(command)

    # run smoothing
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError:
        sh.rmtree(path_temp, ignore_errors=True)
        if create_folder:
            sh.rmtree(path_output, ignore_errors=True)
        print("Intracortical smoothing failed!")

    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
