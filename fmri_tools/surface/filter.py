# -*- coding: utf-8 -*-
"""Filter scalar values on surface mesh."""

import datetime
import functools
import os
import shutil as sh
import subprocess
import warnings
from abc import ABC, abstractmethod
from shutil import copyfile

import numpy as np
from scipy.sparse.linalg import expm, expm_multiply
from scipy.stats import norm

from .. import execute_command
from ..io.filename import get_filename
from .mesh import Mesh

__all__ = [
    "HeatKernel",
    "IterativeNN",
    "Gaussian",
    "LaplacianGaussian",
    "intracortical_smoothing",
    "mris_fwhm",
]


class Filter(ABC):
    """Abstract Base Class for data filtering.

    Base class for the application of a compact spatial filter kernel on scalar
    data defined on a triangular surface mesh. Compact means that the input
    data is weighted with a kernel that has all its weights concentrated in a
    compact local neighborhood.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.

    """

    def __init__(self, verts, faces):
        self.verts = verts
        self.faces = faces
        self._mesh = Mesh(verts, faces)

    @abstractmethod
    def kernel(self):
        """Abstract method for computing kernel weights and a corresponding
        neighbor list."""
        pass

    def apply(self, data, n_iter=1):
        """Iterative application of smoothing kernel (lowpass filter).

        Parameters
        ----------
        data : np.ndarray, shape=(N,)
            Data to be filtered.
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        res : np.ndarray, shape=(N,)
            Filtered data.

        """
        # parameters
        weight, neighbor = self.kernel
        res = np.zeros_like(data)

        # iterative kernel smoothing
        for _ in range(n_iter):
            # add weights
            res = np.zeros_like(data)
            for j in range(self._max_neighbor):
                res += data[neighbor[:, j]] * weight[:, j]

            # initialize new data array
            data = res.copy()

        return res

    def apply_inverse(self, data, n_iter=1):
        """Inverse application of smoothing kernel (highpass filter).

        Parameters
        ----------
        data : np.ndarray, shape=(N,)
            Data to be filtered.
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        np.ndarray, shape=(N,)
            Filtered data.

        """
        return data - self.apply(data, n_iter)

    def apply_noise(self, n_iter=1):
        """Apply smoothing kernel (lowpass filter) to random noise.

        Parameters
        ----------
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        np.ndarray, shape=(N,)
            Filtered random noise.

        """
        return self.apply(self._noise, n_iter)

    def apply_inverse_noise(self, n_iter=1):
        """Inverse application of smoothing kernel (highpass filter) to random
        noise.

        Parameters
        ----------
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        np.ndarray, shape=(N,)
            Filtered random noise.

        """
        data = self._noise

        return data - self.apply(data, n_iter)

    def fwhm(self, n_iter=1):
        """Estimate the full width at half maximum (FWHM) of the filter.

        The kernel size is described by its full width at half maximum (FWHM)
        and is estimated by the equation given in [1]_.

        Parameters
        ----------
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        float
            FWHM of the filter.

        References
        ----------
        .. [1] Hagler, DJ, et al. Smoothing and cluster thresholding for
        cortical surface-based group analysis of fMRI data. Neuroimage 33(4),
        1093--1103 (2005).

        """
        data_filt = self.apply_noise(n_iter)
        edges = self._mesh.edges
        edge_length = self._mesh.avg_edge_length

        var_ds = np.var(data_filt[edges[:, 0]] - data_filt[edges[:, 1]])
        var_s = np.var(data_filt)
        var_ratio = 1 - var_ds / (2 * var_s)

        return edge_length * np.sqrt(-2 * np.log(2) / np.log(var_ratio))

    @property
    def _noise(self):
        """Generate random gaussian noise."""
        return np.random.normal(0, 1, len(self.verts))

    @property
    @functools.lru_cache()
    def _max_neighbor(self):
        """Maximum number of first order neighbors."""
        return np.max(self._mesh.n_neighbors).astype(int)


class HeatKernel(Filter):
    """Heat kernel smoothing.

    Heat kernel smoothing [1]_, [2]_, [3]_ on scalar data defined on a
    triangular surface mesh. The code is mainly adapted from the matlab code by
    Chung et al. [4]_. The kernel bandwidth corresponds to the diffusion time in
    the heat equation [3]_. The FWHM follows 2*sqrt(log(4)*n_iter*sigma) with
    the natural log.

    As pointed out in [4]_, the heat diffusion model becomes less valid as the
    bandwidth increases, and for large bandwidths (which depends on the average
    edge length for the surface), this method becomes indistinguishable from the
    nearest-neighbor method.

    If you use this code, please reference one of the following papers. The
    details on the mathematical basis of of the algorithm can be found there.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.
    sigma : float
        Kernel bandwidth.

    References
    ----------
    .. [1] Chung, MK, et al. Cortical thickness analysis in autism via heat
    kernel smoothing. Neuroimage 25(1), 1256--1265 (2005).
    .. [2] Chung, MK, et al. Unified statistical approach to cortical thickness
    analysis. Inf Process Med Imaging 19, 627--638 (2005).
    .. [3] Chung, MK, et al. Encoding cortical surface by spherical harmonics.
    Statistica Sinica 18, 1269--1291 (2008).
    .. [4] http://pages.stat.wisc.edu/~mchung/softwares/hk/hk_smooth.m

    """

    def __init__(self, verts, faces, sigma):
        super().__init__(verts, faces)
        self.sigma = sigma

    @property
    @functools.lru_cache()
    def kernel(self):
        """Computation of smoothing kernel.

        Returns
        -------
        np.ndarray, shape=(N, self._max_neighbors + 1)
            Heat kernel weighting for each vertex neighbor.
        np.ndarray, shape=(N, self._max_neighbors + 1)
            Correspond array containing vertex neighbors.

        """
        # parameters
        # self._max_neighbor + 1 because the current index is included as well
        nverts = len(self.verts)
        neighbor = np.zeros((nverts, self._max_neighbor + 1), dtype=np.int64)
        weight = np.zeros((nverts, self._max_neighbor + 1))

        for i in range(nverts):
            # vertex neighbors
            nn = self._mesh.neighborhood(i)
            degree = len(nn)

            # get distance to vertex neighbors
            distance = [
                float(np.sum((self.verts[n, :] - self.verts[i, :]) ** 2)) for n in nn
            ]
            distance.insert(0, 0)
            distance = np.array(distance)

            # heat kernel weighting for each neighbor and corresponding
            # neighbors
            weight[i, : 1 + degree] = self._kernel_shape(distance, self.sigma)
            neighbor[i, : 1 + degree] = np.append([i], nn)

        return weight, neighbor

    @staticmethod
    def _kernel_shape(x, s):
        """Heat kernel shape."""
        return np.exp(-x / (2 * s**2)) / np.sum(np.exp(-x / (2 * s**2)))


class IterativeNN(Filter):
    """Iterative nearest neighbor smoothing.

    Implementation of a simple iterative smoothing kernel (averaging box) for
    scalar data defined on a triangular surface mesh. The code follows the
    description in [1]_. In brief, for each iteration, data at each vertex is
    averaged with data from all first order neighbors.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.

    References
    ----------
    .. [1] Hagler, DJ, et al. Smoothing and cluster thresholding for cortical
    surface-based group analysis of fMRI data. Neuroimage 33(4), 1093--1103
    (2005).

    """

    def __init__(self, verts, faces):
        super().__init__(verts, faces)

    @property
    @functools.lru_cache()
    def kernel(self):
        """Computation of smoothing kernel.

        Returns
        -------
        np.ndarray, shape=(N, self._max_neighbors + 1)
            Heat kernel weighting for each vertex neighbor.
        np.ndarray, shape=(N, self._max_neighbors + 1)
            Correspond array containing vertex neighbors.

        """
        # parameters
        # self._max_neighbor + 1 because the current index is included as well
        nverts = len(self.verts)
        neighbor = np.zeros((nverts, self._max_neighbor + 1), dtype=np.int64)
        weight = np.zeros((nverts, self._max_neighbor + 1))

        for i in range(nverts):
            # vertex neighbors
            nn = self._mesh.neighborhood(i)
            degree = len(nn)

            # weighting for each neighbor and corresponding neighbors
            weight[i, : 1 + degree] = 1 / (1 + degree)
            neighbor[i, : 1 + degree] = np.append([i], nn)

        return weight, neighbor


class Gaussian(Filter):
    """Gaussian filter.

    Implementation of a Gaussian filter which can be described by a heat
    diffusion process. The code follows the derivation in [1]_.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.
    t : float
        Gaussian filter size.
    full : bool
        If True, the matrix exponential is computed. If False, the matrix
        exponential will be applied without explicit determination which is
        computationally less demanding.

    References
    ----------
    .. [1] Chen, Yi, et al. Scale-specific analysis of fMRI data on the
    irregular cortical surface. Neuroimage 181, 370--381 (2018).

    """

    def __init__(self, verts, faces, t, full=False):
        super().__init__(verts, faces)
        self.t = t
        self.full = full

    @property
    @functools.lru_cache()
    def kernel(self):
        """Computation of smoothing kernel.

        Computation of the diffusion operator which is defined by the Laplace-
        Beltrami operator. The full kernel needs the computations of a matrix
        exponential which is computationally expensive. Therefore, only the
        exponent can be returned optionally. In this case, the matrix
        exponential of the exponent is applied without explicit determination.

        Returns
        -------
        scipy.sparse.csr.csr_matrix
            Heat diffusion operator.

        """
        exponent = -self.t * self._mesh.laplace_beltrami
        if not self.full:
            return exponent

        return expm(exponent)

    def apply(self, data, n_iter=1):
        """Iterative application of smoothing kernel (lowpass filter). The
        iteration parameter is only included here for consistency with the other
        filter classes. No iterative application was done in the original paper,
        hence n_iter=1.

        Parameters
        ----------
        data : np.ndarray, shape=(N,)
            Data to be filtered.
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        res : np.ndarray, shape=(N,)
            Filtered data.

        """
        res = data.copy()
        for _ in range(n_iter):
            res = self.kernel.dot(res) if self.full else expm_multiply(self.kernel, res)
        return res


class LaplacianGaussian(Gaussian):
    """Laplacian of Gaussian filter.

    Implementation of a Laplacian of Gaussian (LoG) filter which can be
    described by a heat diffusion process. The code follows the derivation in
    [1]_.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.
    t : float
        Gaussian filter size.
    full : bool
        If True, the matrix exponential is computed. If False, the matrix
        exponential will be applied without explicit determination which is
        computationally less demanding.

    References
    ----------
    .. [1] Chen, Yi, et al. Scale-specific analysis of fMRI data on the
    irregular cortical surface. Neuroimage 181, 370--381 (2018).

    """

    def __init__(self, verts, faces, t, full=False):
        super().__init__(verts, faces, t, full)

    def apply(self, data, n_iter=1):
        """Iterative application of bandpass kernel (bandpass filter). The
        iteration parameter is only included here for consistency with the other
        filter classes. No iterative application was done in the original paper,
        hence n_iter=1.

        Parameters
        ----------
        data : np.ndarray, shape=(N,)
            Data to be filtered.
        n_iter : int, optional
            Number of iterations.

        Returns
        -------
        res : np.ndarray, shape=(N,)
            Filtered data.

        """
        res = data.copy()
        laplacian = self._mesh.laplace_beltrami
        for _ in range(n_iter):
            if self.full:
                res = laplacian.dot(self.kernel.dot(res))
            else:
                res = laplacian.dot(expm_multiply(self.kernel, res))

        return res

    def spatial_scale(self, n_iter=1, n_max=None):
        """Estimation of filtered spatial frequency range.

        When applying the bandpass filter, the spatial frequency characterics
        of the filter should be known. Here, we estimate the main spatial
        frequency of the filter by applying the filter to gaussian noise. In the
        filtered noise map, all maxima are then located and their distance to
        the nearest minimum is computed by following the path with maximum
        descent. The distribution of min-max distances is described as a
        log-normal distribution (since distances cannot be negative) and its
        mean and variance are estimated as described in [1]_. The main spatial
        cycle period is then defined as two times the expected value of the
        min-max length distribution.

        Parameters
        ----------
        n_iter : int, optional
            Number of iterations.
        n_max : int, optional
            Number of maxima to be used. If None, all found maxima are taken. If
            less maxima are found than set, all found maxima are used.

        Returns
        -------
        dict
            Dictionary collecting the output under the following keys

            * period (float) : Mean spatial cycle period of the filter.
            * freq (float) : Mean spatial frequency of the filter.
            * mean (float) : Expected value of the log-normal distribution.
            * variance (float) : Variance of the log-normal distribution.
            * length (np.ndarray) : Min-Max distance distribution.

        References
        ----------
        .. [1] http://www.talkstats.com/threads/convert-mean-and-variance-of-lognormal-to-normal-distribution.2757/

        """
        indices = np.arange(len(self.verts))
        data_filt = self.apply_noise(n_iter)

        # find maxima
        ind_max = []
        for ind in indices:
            nn = self._mesh.neighborhood(ind)
            if data_filt[ind] >= np.max(data_filt[nn]):
                ind_max.append(ind)

        # select subset of found maxima
        n_ind_max = len(ind_max)
        if n_max and n_ind_max > n_max:
            ind_max = np.random.choice(ind_max, n_max, replace=False)
        elif n_max and n_ind_max < n_max:
            warnings.warn("Not enough maxima found, using all!")

        # compute min-max distances
        length = np.zeros(n_ind_max)
        for i, ind in enumerate(ind_max):
            l_tmp = 0
            while True:
                nn = self._mesh.neighborhood(ind)
                nn_min = nn[np.argmin(data_filt[nn])]
                if data_filt[nn_min] < data_filt[ind]:
                    l_tmp += self._euclidean_dist(self.verts[ind], self.verts[nn_min])
                    ind = nn_min
                else:
                    length[i] = l_tmp
                    break

        # describe distribution
        length = length[length != 0]  # remove singularities
        data = np.log(length)  # transform to log-normal
        mu, sigma = norm.fit(data)

        mean = np.exp(mu + sigma**2 / 2)
        variance = np.exp(sigma**2 + 2 * mu) * (np.exp(sigma**2) - 1)

        return {
            "period": 2 * mean,
            "freq": 1 / (2 * mean),
            "mean": mean,
            "variance": variance,
            "length": length,
        }

    @staticmethod
    def _euclidean_dist(a, b):
        """Calculate euclidean distance between two points."""
        return np.sqrt(np.sum((a - b) ** 2))


def intracortical_smoothing(
    file_surf,
    file_overlay,
    file_out,
    tan_size=0,
    rad_start=0,
    rad_size=1,
    tan_weights="gauss",
    cleanup=True,
):
    """This function applies the freesurfer function mris_smooth_intracortical
    (available since FreeSurfer 7) which enables simultaneous smoothing in radial and
    tangential direction of the cortical sheet. [1]

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
        raise ValueError(
            "Output file name is expected to have the file " "extension mgh or mgz!"
        )

    # make output folder
    create_folder = 0
    if not os.path.exists(path_output):
        create_folder = 1
        os.makedirs(path_output)

    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_" + tmp_string)
    path_surf = os.path.join(path_temp, "surf")
    path_overlay = os.path.join(path_temp, "overlay")

    # make temporary folder
    if os.path.exists(path_temp):
        raise FileExistsError("Temporary folder already exists!")

    os.mkdir(path_temp)
    os.mkdir(path_surf)
    os.mkdir(path_overlay)

    # copy input files to temporary folder
    for i, f in enumerate(file_surf):
        name_tmp = "surf_" + str(i).zfill(4)
        copyfile(f, os.path.join(path_surf, name_tmp))

    for i, f in enumerate(file_overlay):
        name_tmp = "overlay_" + str(i).zfill(4) + ext_output
        copyfile(f, os.path.join(path_overlay, name_tmp))

    # unix command
    command = [
        "mris_smooth_intracortical",
        "--surf_dir",
        path_surf,
        "--surf_name",
        "surf_*",
        "--overlay_dir",
        path_overlay,
        "--overlay_name",
        "overlay_*.mgh",
        "--output_dir",
        path_output,
        "--output_name",
        name_output + ext_output,
        "--tan-size",
        str(tan_size),
        "--rad-size",
        str(rad_size),
        "--rad-start",
        str(rad_start),
        "--tan-weights",
        str(tan_weights),
    ]

    command = " ".join(command)

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


def mris_fwhm(file_in, file_out, sub, fwhm):
    """Smooth surface overlay using FreeSurfer.

    Parameters
    ----------
    file_in : str
        File name of input file.
    file_out : str
        File name of output file.
    sub : str
        FreeSurfer subject name.
    fwhm : float
        Filter size.

    Raises
    ------
    ValueError
        Wrong file extension.
    ValueError
        Unknown hemisphere.
    """
    _, name_in, ext_in = get_filename(file_in)
    path_out, _, _ = get_filename(file_out)

    if ext_in != ".mgh":
        raise ValueError("File extension should be *.mgh!")

    hemi = name_in[:2]
    if hemi not in ["lh", "rh"]:
        raise ValueError("Unknown hemisphere!")

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "mris_fwhm"
    command += f" --s {sub}"
    command += f" --hemi {hemi}"
    command += " --smooth-only"
    command += f" --fwhm {fwhm}"
    command += f" --i {file_in}"
    command += f" --o {file_out}"

    # run
    execute_command(command)
