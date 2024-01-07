# -*- coding: utf-8 -*-
"""Test filename utilities."""

import pytest

from fmri_tools.io.filename import get_filename, natsort


@pytest.fixture
def filenames():
    """Test data for get_filename."""
    return [
        ("/home/test/test_data.nii", ".nii"),
        ("/home/test/test_data.nii.gz", ".nii.gz"),
        ("/home/test/test_data.mgh", ".mgh"),
        ("/home/test/test_data.mgz", ".mgz"),
    ]


@pytest.fixture
def name_list():
    """Test data for natsort."""
    return [
        (["Run_4", "Run_10", "Run_0"], ["Run_0", "Run_4", "Run_10"]),
        (["Ru4n", "Ru0n", "Ru10n"], ["Ru0n", "Ru4n", "Ru10n"]),
        (["4run", "8run", "1run"], ["1run", "4run", "8run"]),
        (["4run2", "8run0", "1run7"], ["1run7", "4run2", "8run0"]),
    ]


def test_get_filename(filenames):
    """Asserts the retrieval of correct file extensions."""
    for f_ in filenames:
        file_name = f_[0]
        ext_expected = f_[1]
        _, _, ext = get_filename(file_name)
        assert ext == ext_expected


def test_natsort(name_list):
    """Assert natural sorting."""
    for n_ in name_list:
        lst_in = n_[0]
        lst_sorted = n_[1]
        lst_out = natsort(lst_in)
        assert lst_out == lst_sorted
