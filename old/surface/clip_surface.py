# -*- coding: utf-8 -*-

import numpy as np
from nibabel.freesurfer.io import read_label

from ..io.surf import read_mgh, write_mgh

__all__ = ["clip_surface", "clip_mgh"]


def clip_surface(verts, faces, ind_keep):
    """Clip surface.

    Clip a triangle surface mesh based on a separate index list. Indices should
    form a connected region on the surface. Vertices which are not in the index
    list will be removed and the corresponding face array will be updated.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.
    ind_keep : list
        Index list of vertices to keep, which should be a connected region.

    Returns
    -------
    verts : np.ndarray, shape=(N,3)
        Remaining vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of remaining triangles.

    """

    # get new vertices
    verts = verts[ind_keep, :]

    # get new faces
    faces_keep = np.zeros(len(faces))
    faces_keep += np.in1d(faces[:, 0], ind_keep)
    faces_keep += np.in1d(faces[:, 1], ind_keep)
    faces_keep += np.in1d(faces[:, 2], ind_keep)
    faces = faces[faces_keep == 3, :]

    # reindex faces
    for i, ind in enumerate(ind_keep):
        faces[faces == ind] = i

    return verts, faces


def clip_mgh(mgh_in, label_in, mgh_out):
    """Clip mgh overlay.

    Parameters
    ----------
    mgh_in : str
        Path to input mgh file.
    label_in : str
        Path to input label file.
    mgh_out : str
        Path to output mgh file.

    Returns
    -------
    None.

    """

    arr, affine, header = read_mgh(mgh_in)
    label = read_label(label_in)
    write_mgh(mgh_out, arr[label], affine, header)
