# -*- coding: utf-8 -*-
"""Python package for for analysis of high-resolution fMRI data."""

import os


class MatlabCommand:
    """Class for running Matlab commands.

    Parameters
    ----------
    fname : str
        Name of the Matlab script to run.
    *args
        Variable length argument list to pass to the Matlab script.

    """

    def __init__(self, fname, *args):
        self.fname = fname
        self.args = [x for x in args]

        # basename of matlab file
        if self.fname.endswith(".m"):
            self.fname = self.fname[:-2]

        # path to matlab folder
        self.path = os.path.abspath(os.path.dirname(__file__))

    def run(self):
        """Run the Matlab command."""

        # check if matlab file exists and print matlab version
        self._check_exists()
        self._matlab_version()

        cmd_args = ""
        for arg in self.args:
            if isinstance(arg, str):
                cmd_args += "\'" + arg + "\',"
            else:
                cmd_args += str(arg) + ","

        # change to matlab folder and execute script
        path_current = os.path.abspath(os.getcwd())
        os.chdir(self.path)
        command = "matlab -nodisplay -nodesktop -r \"" + \
                  self.fname + "(" + cmd_args[:-1] + "); exit;\""
        print("Execute: " + command)
        os.system(command)
        os.chdir(path_current)

    def _check_exists(self):
        """Check if the Matlab file exists."""

        if not os.path.isfile(os.path.join(self.path, self.fname + ".m")):
            raise Exception("Matlab file does not exist")

    @staticmethod
    def _matlab_version():
        """Get the version of Matlab."""

        cmd = "echo $(basename $(matlab -n | grep ' MATLAB ' | cut -d= -f2))"
        version = os.popen(cmd).read().rstrip()
        print("Found Matlab version: " + version)
