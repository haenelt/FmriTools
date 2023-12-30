# -*- coding: utf-8 -*-
"""Test filename utilities."""

import pytest

from fmri_tools.io.filename import get_filename


@pytest.fixture
def filenames():
    """Test data for get_filename."""
    return [
        ("/home/test/test_data.nii", ".nii"),
        ("/home/test/test_data.nii.gz", ".nii.gz"),
        ("/home/test/test_data.mgh", ".mgh"),
        ("/home/test/test_data.mgz", ".mgz"),
    ]


def test_get_filename(filenames):
    """Asserts the retrieval of correct file extensions."""
    for f_ in filenames:
        file_name = f_[0]
        ext_expected = f_[1]
        _, _, ext = get_filename(file_name)
        assert ext == ext_expected
