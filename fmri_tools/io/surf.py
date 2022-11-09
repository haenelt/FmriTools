# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.mghformat import MGHHeader
from nibabel.freesurfer.io import write_geometry, read_geometry, read_label, \
    read_morph_data, write_morph_data
from scipy.spatial import Delaunay
from gbb.neighbor.nn_2d import nn_2d

# local inputs
from .get_filename import get_filename

__all__ = ['write_mgh', 'read_mgh', 'write_label', 'read_patch', 'patch_as_mesh',
           'mgh_to_patch', 'curv_to_patch', 'label_to_patch', 'write_vector_field',
           'write_white2pial']


def write_mgh(file_out, arr, affine=None, header=None):
    """Write MGH.

    This function adds two empty dimensions to an array and saves it as a
    freesurfer mgh surface file.

    Parameters
    ----------
    file_out : str
        Filename of output file.
    arr : ndarray
        Image array.
    affine : ndarray, optional
        Affine transformation matrix. The default is None.
    header : MGHHeader, optional
        Image header. The default is None.

    Raises
    ------
    ValueError
        If `file_out` is not a string or has a file extension which is not
        supported.

    Returns
    -------
    None.

    """

    # check filename
    if not isinstance(file_out, str):
        raise ValueError("Filename must be a string!")

    if not file_out.endswith("mgh"):
        raise ValueError("Currently supported file format is mgh.")

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # add empty dimensions
    arr = np.expand_dims(arr, axis=1)
    arr = np.expand_dims(arr, axis=1)

    if affine is None:
        affine = np.eye(4)

    if header is None:
        header = MGHHeader()

    # write output
    output = nb.Nifti1Image(arr, affine, header)
    nb.save(output, file_out)


def read_mgh(file_in, read_affine=False, read_header=False):
    """Read MGH.

    This function reads a surface mgh file and removes empty dimensions from the
    data array.

    Parameters
    ----------
    file_in : str
        File name of input file.
    read_affine : bool
        If True, return affine transformation matrix.
    read_header : bool
        If True, return header information.

    Raises
    ------
    ValueError
        If `file_in` is not a string or has a file extension which is not
        supported.

    Returns
    -------
    arr : ndarray
        Image array.
    affine : ndarray
        Affine transformation matrix. Returned only if `read_affine` is True.
    header : MGHHeader
        Image header. Returned only if `read_header` is True.

    """

    # check filename
    if not isinstance(file_in, str):
        raise ValueError("Filename must be a string!")

    if not file_in.endswith("mgh"):
        raise ValueError("Currently supported file format is mgh.")

    # get header
    header = nb.load(file_in).header
    affine = nb.load(file_in).affine

    # get data
    res = nb.load(file_in).get_fdata()
    res = np.squeeze(res)

    if read_affine:
        res += (affine,)

    if read_header:
        res += (header,)

    return res


def write_label(file_out, arr_label):
    """Write label.

    This function writes a textfile which can be read as label file in
    freesurfer.

    Parameters
    ----------
    file_out : str
        Filename of label file.
    arr_label : list
        List of label indices.

    Raises
    ------
    ValueError
        If `file_out` is not a string or has a file extension which is not
        supported.

    Returns
    -------
    None.

    """

    # check filename
    if not isinstance(file_out, str):
        raise ValueError("Filename must be a string!")

    if not file_out.endswith("label"):
        raise ValueError("Currently supported file format is txt.")

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # number of labels
    n_label = len(arr_label)

    with open(file_out, "w") as f:
        f.write('#!ascii label  , from subject  vox2ras=TkReg\n')
        f.write(str(n_label) + '\n')
        for i in range(n_label):
            f.write(str(arr_label[i]) + ' 0.000 0.000 0.000 0.000\n')


def read_patch(file_in):
    """Read patch.

    This function reads an freesurfer patch saved in binary format. It is an
    equivalent to the matlab function read_patch in the ./freesurfer/matlab
    folder. Data is read in big endian order.

    Parameters
    ----------
    file_in : str
        Full path of the input file.

    Returns
    -------
    x : ndarray
        x-coordinates of patch.
    y : ndarray
        y-coordinates of patch.
    z : ndarray
        z-coordinates of patch (if flattened, there should be only zeros).
    ind : ndarray
        index of each node.

    """

    # load data
    data_array_int = np.fromfile(file_in, np.dtype(">i"))
    data_array_float = np.fromfile(file_in, np.dtype(">f"))

    # check version
    ver = data_array_int[0]
    if ver != -1:
        sys.exit("incorrect version # " + str(ver) + " (not -1) found in file")

    # size of data array
    data_size = data_array_int[1]

    # reshape data into array
    data_array_int = data_array_int[2:]
    data_array_float = data_array_float[2:]
    data_array_int = data_array_int.reshape(data_size, 4)
    data_array_float = data_array_float.reshape(data_size, 4)

    # get vertex indices and coordinates
    ind = data_array_int[:, 0]
    ind[ind < 0] = -ind[ind < 0] - 1
    ind[ind >= 0] = ind[ind >= 0] - 1

    x = data_array_float[:, 1]
    y = data_array_float[:, 2]
    z = data_array_float[:, 3]

    return x, y, z, ind


def patch_as_mesh(file_out, file_patch):
    """Patch as Mesh.

    This function reads coordinates from a flattened freesurfer patch and creates a
    triangular surface mesh using Delaunay triangulation.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """

    x, y, _, ind = read_patch(file_patch)
    coords = np.zeros((len(x), 2))
    coords[:, 0] = x
    coords[:, 1] = y

    tri = Delaunay(coords)
    fac = tri.simplices

    vtx = np.zeros((len(x), 3))
    vtx[:, 0] = x
    vtx[:, 1] = y

    write_geometry(file_out, vtx, fac)


def mgh_to_patch(file_out, file_mgh, file_patch):
    """MGH to Patch.

    This function reads an MGH overlay and saves only coordinates that match with the
    corresponding freesurfer patch.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_mgh : str
        Filename of mgh overlay.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """

    _, _, _, ind = read_patch(file_patch)
    arr, affine, header = read_mgh(file_mgh)
    arr = arr[ind]
    write_mgh(file_out, arr, affine=affine, header=header)


def curv_to_patch(file_out, file_curv, file_patch):
    """Curv to Patch.

    This function reads a freesurfer curvature file and saves only coordinates that
    match with the corresponding freesurfer patch.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_curv : str
        Filename of freesurfer curvature file.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """

    _, _, _, ind = read_patch(file_patch)
    curv = read_morph_data(file_curv)
    curv = curv[ind]
    write_morph_data(file_out, curv)


def label_to_patch(file_out, file_label, file_patch):
    """Label to Patch.

    This function reads a freesurfer label file and saves only coordinates that match
    with the corresponding freesurfer patch.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_label : str
        Filename of freesurfer label file.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """

    _, _, _, ind = read_patch(file_patch)
    label = read_label(file_label)
    label_new = [np.where(ind == i)[0][0] for i in label]
    write_label(file_out, label_new)


def write_vector_field(vtx0, vtx1, adjm, file_out, meta_data=None, step_size=100,
                       shape="line"):
    """Write vector field.

    This function generates a surface mesh to visualize a vector field as a
    triangular mesh.

    Parameters
    ----------
    vtx0 : ndarray
        Array of vector start points.
    vtx1 : ndarray
        Array of vector end points.
    adjm : ndarray
        Adjacency matrix.
    file_out : str
        Filename of output surface mesh.
    meta_data : dict-like or None
        Header information. See documentation of argument volume_info in
        nibabel.freesurfer.io.write_geometry for more information. The default is None.
    step_size : int, optional
        Vector subset which will be visualized. The default is 100.
    shape : str, optional
        Line, triangle, prism. The default is "line".

    Returns
    -------
    None.

    """

    # array containing a list of considered vectors
    t = np.arange(0, len(vtx0), step_size)

    # initialise faces for specific shape
    if shape == "prism":
        f_new = [[0, 1, 2],
                 [3, 4, 5],
                 [0, 1, 4],
                 [0, 3, 4],
                 [1, 2, 5],
                 [1, 4, 5],
                 [0, 2, 5],
                 [0, 3, 5]]
        f_iter = 6
    elif shape == "triangle":
        f_new = [[0, 1, 2]]
        f_iter = 3
    elif shape == "line":
        f_new = [[0, 1, 0]]
        f_iter = 2

    v_res = []
    f_res = []
    for i in range(len(t)):

        # get index from nearest neighbour of a given vertex
        nn = nn_2d(t[i], adjm, 0)
        nn = nn[:2]

        # get all vertex points for specific shape
        if shape == "prism":
            A = list(vtx0[t[i]])
            B = list(vtx0[nn[0]])
            C = list(vtx0[nn[1]])
            D = list(vtx1[t[i]])
            E = list(vtx1[nn[0]])
            F = list(vtx1[nn[1]])
            v_new = [A, B, C, D, E, F]
        elif shape == "triangle":
            A = list(vtx0[t[i]])
            B = list(vtx0[nn[0]])
            C = list(vtx1[t[i]])
            v_new = [A, B, C]
        elif shape == "line":
            A = list(vtx0[t[i]])
            B = list(vtx1[t[i]])
            v_new = [A, B]

        # update faces
        if i > 0:
            for j in range(len(f_new)):
                f_new[j] = [x + f_iter for x in f_new[j]]

        # update resulting vertex and face list
        v_res.extend(v_new)
        f_res.extend(f_new)

    # vertices and faces as array
    v_res = np.array(v_res)
    f_res = np.array(f_res)

    # write output geometry
    write_geometry(file_out, v_res, f_res, volume_info=meta_data)


def write_white2pial(file_out, file_white, file_pial, adjm, step_size=100,
                     shape="line"):
    """Plot white to pial.

    This function generates lines between corresponding vertices at the white
    and pial surface to visualize the shift between matched vertices caused by
    realigning surfaces independently. You can either construct prisms,
    triangles or lines.

    Parameters
    ----------
    file_out : str
        Filename of output surface.
    file_white : str
        Filename of white surface.
    file_pial : str
        Filename of pial surface.
    adjm : obj
        Adjacency matrix.
    step_size : int, optional
        Subset of vertices.
    shape : str, optional
        line, triangle, prism.

    Returns
    -------
    None.

    """

    # read geometry
    vtx_white, fac_white, header_white = read_geometry(file_white,
                                                       read_metadata=True)
    vtx_pial, fac_pial = read_geometry(file_pial)

    # array containing a list of considered vertices
    t = np.arange(0, len(vtx_white), step_size)

    # initialise faces for specific shape
    if shape == "prism":
        fac_new = [[0, 1, 2],
                   [3, 4, 5],
                   [0, 1, 4],
                   [0, 3, 4],
                   [1, 2, 5],
                   [1, 4, 5],
                   [0, 2, 5],
                   [0, 3, 5]]
        fac_iter = 6
    elif shape == "triangle":
        fac_new = [[0, 1, 2]]
        fac_iter = 3
    elif shape == "line":
        fac_new = [[0, 1, 0]]
        fac_iter = 2

    vtx_res = []
    fac_res = []
    for i in range(len(t)):

        # get index from nearest neighbour of a given vertex
        nn = nn_2d(t[i], adjm, 0)
        nn = nn[:2]

        # get all vertex points for specific shape
        if shape == "prism":
            A = list(vtx_white[t[i]])
            B = list(vtx_white[nn[0]])
            C = list(vtx_white[nn[1]])
            D = list(vtx_pial[t[i]])
            E = list(vtx_pial[nn[0]])
            F = list(vtx_pial[nn[1]])
            vtx_new = [A, B, C, D, E, F]
        elif shape == "triangle":
            A = list(vtx_white[t[i]])
            B = list(vtx_white[nn[0]])
            C = list(vtx_pial[t[i]])
            vtx_new = [A, B, C]
        elif shape == "line":
            A = list(vtx_white[t[i]])
            B = list(vtx_pial[t[i]])
            vtx_new = [A, B]

        # update faces
        if i > 0:
            for j in range(len(fac_new)):
                fac_new[j] = [x + fac_iter for x in fac_new[j]]

        # update resulting vertex and face list
        vtx_res.extend(vtx_new)
        fac_res.extend(fac_new)

    # vertices and faces as array
    vtx_res = np.array(vtx_res)
    fac_res = np.array(fac_res)

    # write output geometry
    write_geometry(file_out, vtx_res, fac_res, volume_info=header_white)
