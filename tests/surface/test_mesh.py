# -*- coding: utf-8 -*-
"""Test mesh utilities."""

import os

import pytest

from fmri_tools.surface.mesh import Mesh

# file name of sphere in data folder
FILE_SPHERE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "sphere")


def test_n_neighbors():
    """Assert number of vertex neighbors."""
    mesh = Mesh.from_file(FILE_SPHERE)
    n_neighbors = mesh.n_neighbors
    n_neighbors_expected = [
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        5.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
        6.0,
    ]
    for i, j in zip(n_neighbors, n_neighbors_expected):
        assert i == j


def test_neighborhood():
    """Assert neighborhood around vertex."""
    mesh = Mesh.from_file(FILE_SPHERE)
    neighborhood = mesh.neighborhood(0)
    neighborhood_expected = [12, 14, 16, 18, 20]
    for i, j in zip(neighborhood, neighborhood_expected):
        assert i == j
