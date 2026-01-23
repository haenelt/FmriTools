# -*- coding: utf-8 -*-
"""Motion correction of fMRI time series using afni and optional application to
different echoes in ME-fMRI projects."""

import uuid
import shutil
from pathlib import Path
from joblib import Parallel, delayed
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from ..registration.afni import volreg, allineate, prepare_header, extract_ref

__all__ = ["MotionCorrection"]


class MotionCorrection:
    """Apply motion correction to a multi-run session.

    Parameters
    ----------
    fnames : tuple
        List of input time series file names. All time series are realigned to the first
        file.
    dir_out : str, optional
        Output directory for motion estimates.

    Examples
    --------
    Basic usage:

    >>> fnames = ("run1.nii", "run2.nii", "run3.nii")
    >>> mc = MotionCorrection(fnames)
    >>> mc()

    With application to other datasets, e.g. different echoes:

    >>> fnames = ("run1.nii", "run2.nii", "run3.nii")
    >>> paths_other = ("/path1", "/path2", "/path3")
    >>> mc = MotionCorrection(fnames)
    >>> mc.estimate()
    >>> for i, path in enumerate(paths_other):
    >>>     files = list(Path(path).glob("*.nii"))
    >>>     mc.apply_other(i, files)

    """

    # folder names for motion estimates and for realigned time series
    NAME_MOCO_SINGLE = "moco"
    NAME_MOCO = "moco"

    # number of cores for parallel application of motion parameters
    N_CORE = 4

    def __init__(self, fnames, dir_out=None):
        self.fnames = fnames
        self.dir_out = dir_out if dir_out else str(Path(self.fnames[0]).parent.parent)
        self.nrun = len(self.fnames)
        # Get extension from first file name
        self.ext = "".join(Path(self.fnames[0]).suffixes)

    def __call__(self):
        """Run the full motion correction pipeline."""
        self.estimate()
        self.apply()
        self.plot_summary()
        self.make_mean()
        self.make_vol1()

    @property
    def file_ref(self):
        """Get reference volume from first run. Equivalent to other time series, the
        data is deobliqued before extraction, which is applied to a copied temporary
        file to not change the original data file."""
        _file_ref = self.path_moco(0) / f"ref{self.ext}"
        if not _file_ref.exists():
            file_tmp = self.path_moco(0) / f"tmp{self.ext}"
            shutil.copyfile(self.fnames[0], file_tmp)
            prepare_header(file_tmp)
            extract_ref(file_tmp, _file_ref)
            file_tmp.unlink(missing_ok=True)
        return _file_ref

    @property
    def dim(self):
        """Get image dimensions of time series (nx, ny, nz, nt)."""
        return nb.load(self.fnames[0]).header["dim"][1:5]

    def file_res(self, run):
        """File name of volreg output time series."""
        return self.path_moco(run) / f"res{self.ext}"

    def file_moco(self, run):
        """File name of afni realignment matrix needed to apply motion estimates."""
        return self.path_moco(run) / "moco_matrix.1D"

    def file_param(self, run):
        """File name of afni motion parameters needed to plot motion estimates."""
        return self.path_moco(run) / "moco_params.1D"

    def path_in(self, run):
        """Path to single runs."""
        if run >= self.nrun or run < 0:
            raise ValueError("Invalid run!")
        return Path(self.fnames[run]).parent

    def path_out(self, run):
        """Path to output folder of single runs."""
        _path = self.path_in(run) / self.NAME_MOCO_SINGLE
        _path.mkdir(exist_ok=True, parents=True)
        return _path

    def path_moco(self, run):
        """Path to moco folder of single runs."""
        _path = Path(self.dir_out) / self.NAME_MOCO / f"Run_{run + 1:02d}"
        _path.mkdir(exist_ok=True, parents=True)
        return _path

    def path_summary(self):
        """Path to folder containing information about realignment procedure."""
        _path = Path(self.dir_out) / self.NAME_MOCO / "summary"
        _path.mkdir(exist_ok=True, parents=True)
        return _path

    def estimate(self):
        """Estimate motion parameters. The reference volume is taken from the first time
        series. All data sets are deobliqued before motion correction, which is applied
        to a temporay file to not change the original file."""
        for i, fname in enumerate(self.fnames):
            file_in = self.path_moco(i) / f"in{self.ext}"
            shutil.copyfile(fname, file_in)
            prepare_header(file_in)
            file_out = self.path_moco(i) / f"out{self.ext}"
            volreg(file_in, file_out, self.file_ref)
            allineate(file_in, self.file_res(i), self.file_ref, self.file_moco(i))

    def apply(self):
        """Apply estimated motion parameters to input time series with final
        interpolation."""
        for i, fname in enumerate(self.fnames):
            self._apply_transform(i, fname)

    def apply_other(self, run, fnames):
        """Apply estimated motion parameters to different datasets listed in fnames.
        These could be different echoes in a multi-echo fMRI acquisition."""
        Parallel(n_jobs=self.N_CORE)(
            delayed(self._apply_transform)(run, fname) for fname in fnames
        )

    def _apply_transform(self, run, filename):
        """Apply estimated motion parameters from a specific run to one file. The same
        preprocessing is done as for the motion estimation."""
        file_in = self.path_out(run) / f"tmp_{uuid.uuid4()}{self.ext}"
        file_out = self.path_out(run) / f"u{Path(filename).name}"
        shutil.copyfile(filename, file_in)
        prepare_header(file_in)
        allineate(file_in, file_out, self.file_ref, self.file_moco(run))
        file_in.unlink(missing_ok=True)

    def plot_summary(self):
        """Plot estimated motion parameters as line plots across runs."""
        motion_data = []
        for i, _ in enumerate(self.fnames):
            _trace = np.loadtxt(self.file_param(i))
            motion_data.extend(_trace)
        motion_data = np.array(motion_data)

        fig, ax = plt.subplots()
        ax.plot(motion_data[:, 0], label="roll")
        ax.plot(motion_data[:, 1], label="pitch")
        ax.plot(motion_data[:, 2], label="yaw")
        ax.set_xlabel("Time in TR")
        ax.set_ylabel("Rotation in deg")
        ax.legend()
        fig.savefig(self.path_summary() / "rotation.svg")

        fig, ax = plt.subplots()
        ax.plot(motion_data[:, 3], label="superior direction")
        ax.plot(motion_data[:, 4], label="left direction")
        ax.plot(motion_data[:, 5], label="posterior direction")
        ax.set_xlabel("Time in TR")
        ax.set_ylabel("Translation in mm")
        ax.legend()
        fig.savefig(self.path_summary() / "translation.svg")

    def make_mean(self):
        """Make mean volume across runs."""
        nx, ny, nz, _ = self.dim
        _mean = np.zeros((nx, ny, nz))
        for i, _ in enumerate(self.fnames):
            _data = nb.load(self.file_res(i))
            _mean += np.mean(_data.get_fdata(), axis=3)
        _mean /= self.nrun
        _data0 = nb.load(self.file_res(0))
        output = nb.Nifti1Image(_mean, _data0.affine, _data0.header)
        nb.save(output, self.path_summary() / f"mean{self.ext}")

    def make_vol1(self):
        """Make time series of first volumes of each run."""
        nx, ny, nz, _ = self.dim
        _vol1 = np.zeros((nx, ny, nz, self.nrun))
        for i, _ in enumerate(self.fnames):
            _data = nb.load(self.file_res(i))
            _vol1[:, :, :, i] = _data.get_fdata()[:, :, :, 0]
        _data0 = nb.load(self.file_res(0))
        output = nb.Nifti1Image(_vol1, _data0.affine, _data0.header)
        nb.save(output, self.path_summary() / f"vol1{self.ext}")
