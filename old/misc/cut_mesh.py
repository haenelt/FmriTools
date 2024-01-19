# -*- coding: utf-8 -*-
"""
Surface mesh across cortical depth

This scripts cuts through a mesh to show the cortex in its depth. To do so,
cortical depth-dependent vertex points and faces are added after the cut. As a
last step, a new overlay is created to show depth-dependent values.

"""

import os
import sys

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import read_geometry, read_label, write_geometry
from numpy.linalg import norm

from ..surface.mesh import Mesh
from ..surface.smooth import mris_smooth

# input files
file_in = "/home/daniel/Schreibtisch/AAA_tmp/lh.mesh"
file_in2 = ""
file_contrast = [
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer0.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer1.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer2.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer3.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer4.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer5.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer6.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer7.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer8.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer9.mgh",
    "/home/daniel/Schreibtisch/AAA_tmp/lh.Z_all_left_right_GE_EPI2_upsampled_layer10.mgh",
]
file_label = ""  # "/home/daniel/Schreibtisch/movie/data/occ2.label"
path_output = "/home/daniel/Schreibtisch/test4"

# parameters
s = 1  # scale parameter for vertex shifts
n_smooth = 5  # final smoothing iterations

# cut coordinates
x_min = None
x_max = None
y_min = -18
y_max = None
z_min = -42
z_max = None
merge = "union"

# do not edit below


def cut_mesh(
    vertex,
    x_min=None,
    x_max=None,
    y_min=None,
    y_max=None,
    z_min=None,
    z_max=None,
    merge="union",
):
    """This function labels vertices which are outside of given vertex coordinates. If
    more than one threshold coordinate is used, the resulting label list can be
    calculating either with a union or with an intersection operation.

    Parameters
    ----------
    vertex: np.ndarray, shape=(N,3)
        Vertex array.
    x_min: float, optional
        Minimum x-coordiante.
    x_max: float, optional
        Maximum x-coordiante.
    y_min: float, optional
        Minimum y-coordiante.
    y_max: float, optional
        Maximum y-coordiante.
    z_min: float, optional
        Minimum z-coordiante.
    z_max: float, optional
        Maximum z-coordiante.
    merge: str, optional
        Merging operation (union or intersection).

    Returns
    -------
    label : list
        Resulting label list.

    """

    # initialize label list
    label = []

    # get single label arrays for each threshold coordinate
    if x_min is not None:
        label.append(np.argwhere(vertex[:, 0] < x_min)[:, 0])

    if x_max is not None:
        label.append(np.argwhere(vertex[:, 0] > x_max)[:, 0])

    if y_min is not None:
        label.append(np.argwhere(vertex[:, 1] < y_min)[:, 0])

    if y_max is not None:
        label.append(np.argwhere(vertex[:, 1] > y_max)[:, 0])

    if z_min is not None:
        label.append(np.argwhere(vertex[:, 2] < z_min)[:, 0])

    if z_max is not None:
        label.append(np.argwhere(vertex[:, 2] > z_max)[:, 0])

    # merge labels
    nlabel = len(label)
    if nlabel == 1:
        label = list(set(label))

    elif nlabel == 2:
        if merge == "union":
            label = list(set(label[0]).union(set(label[1])))
        elif merge == "intersection":
            label = list(set(label[0]).intersection(set(label[1])))
        else:
            sys.exit("No valid merge parameter!")

    elif nlabel == 3:
        if merge == "union":
            label_tmp = list(set(label[0]).union(set(label[1])))
            label = list(set(label_tmp).union(set(label[2])))
        elif merge == "intersection":
            label_tmp = list(set(label[0]).intersection(set(label[1])))
            label = list(set(label_tmp).intersection(set(label[2])))
        else:
            sys.exit("No valid merge parameter!")

    elif nlabel == 4:
        if merge == "union":
            label_tmp = list(set(label[0]).union(set(label[1])))
            label_tmp = list(set(label_tmp).union(set(label[2])))
            label = list(set(label_tmp).union(set(label[3])))
        elif merge == "intersection":
            label_tmp = list(set(label[0]).intersection(set(label[1])))
            label_tmp = list(set(label_tmp).intersection(set(label[2])))
            label = list(set(label_tmp).intersection(set(label[3])))
        else:
            sys.exit("No valid merge parameter!")

    elif nlabel == 5:
        if merge == "union":
            label_tmp = list(set(label[0]).union(set(label[1])))
            label_tmp = list(set(label_tmp).union(set(label[2])))
            label_tmp = list(set(label_tmp).union(set(label[3])))
            label = list(set(label_tmp).union(set(label[4])))
        elif merge == "intersection":
            label_tmp = list(set(label[0]).intersection(set(label[1])))
            label_tmp = list(set(label_tmp).intersection(set(label[2])))
            label_tmp = list(set(label_tmp).intersection(set(label[3])))
            label = list(set(label_tmp).intersection(set(label[4])))
        else:
            sys.exit("No valid merge parameter!")

    elif nlabel == 6:
        if merge == "union":
            label_tmp = list(set(label[0]).union(set(label[1])))
            label_tmp = list(set(label_tmp).union(set(label[2])))
            label_tmp = list(set(label_tmp).union(set(label[3])))
            label_tmp = list(set(label_tmp).union(set(label[4])))
            label = list(set(label_tmp).union(set(label[5])))
        elif merge == "intersection":
            label_tmp = list(set(label[0]).intersection(set(label[1])))
            label_tmp = list(set(label_tmp).intersection(set(label[2])))
            label_tmp = list(set(label_tmp).intersection(set(label[3])))
            label_tmp = list(set(label_tmp).intersection(set(label[4])))
            label = list(set(label_tmp).intersection(set(label[5])))
        else:
            sys.exit("No valid parameter!")

    print("Number of cuts: " + str(nlabel))
    print("Number of label points from cut: " + str(len(label)))

    return label


def remove_face(face, index_keep, index_remove):
    """Faces which contain vertex indices not in index_keep are removed from the face
    array. The resulting face array is then resorted.

    Parameters
    ----------
    face: np.ndarray, shape=(N,3)
        Face array.
    index_keep: np.ndarray, shape=(U,)
        Index array of vertex indices to keep.
    index_remove: np.ndarray, shape=(V,)
        Index array of vertex indices to remove.

    Returns
    -------
    face_new: np.ndarray, shape=(M,3)
        New face array.

    """

    # remove faces
    face_keep = np.zeros(len(face[:, 0]))
    face_keep += np.in1d(face[:, 0], index_keep)
    face_keep += np.in1d(face[:, 1], index_keep)
    face_keep += np.in1d(face[:, 2], index_keep)
    face_new = face[face_keep == 3, :]

    # reindex faces
    n = len(index_remove)
    for i in range(n):
        # current status
        print("sort faces: " + str(i) + " of " + str(n))

        tmp = face_new[face_new >= index_remove[i]] - 1
        face_new[face_new >= index_remove[i]] = tmp

    return face_new


def get_edge(index_keep, vtx, fac):
    """After the surface mesh is cut, vertex indices lying at the surface mesh border
    are listed.

    Parameters
    ----------
    index_keep: list
        Index array of vertex indices to keep.
    vtx: ndarray
        Vertex array.
    fac: ndarray
        Corresponding face array.

    Returns
    -------
    index_edge: list
        Index array of vertex indices at surface mesh border.

    """

    index_edge = []
    n = len(index_keep)
    for i in range(n):
        # current status
        print("search edges: " + str(i) + " of " + str(n))

        # get first order neighbors
        nn = Mesh(vtx, fac).neighborhood(index_keep[i])

        j = 0
        while True:
            if nn[j] not in index_keep:
                index_edge.append(i)
                break

            if j == len(nn) - 1:
                break

            j += 1

    return index_edge


def nn_3d(vtx0, vtx, r_size):
    """This function computes all nearest neighbors found within a sphere with a radius
    defined in ras coordinates. Note that the defined neighborhood does not have to be
    fully connected in this case. Vertex coordinates are in ras space.

    Parameters
    ----------
    vtx0 : ndarray
        Vertex point.
    vtx : ndarray
        Array of vertices.
    r_size : float
        Radius of sphere in ras coordinates.

    Returns
    -------
    nn : ndarray
        Array of neighbor indices.
    r[nn] : ndarray
        Euclidean distance to neighbors.

    """
    rx = (vtx[:, 0] - vtx0[0]) ** 2
    ry = (vtx[:, 1] - vtx0[1]) ** 2
    rz = (vtx[:, 2] - vtx0[2]) ** 2

    r = np.sqrt(rx + ry + rz)

    nn = np.where(r < r_size)[0]

    return nn, r[nn]


# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

# filenames of output
file_out0 = os.path.join(path_output, "surf0")
file_out1 = os.path.join(path_output, "surf1")
file_out0_cut = os.path.join(path_output, "surf0_cut")
file_out1_cut = os.path.join(path_output, "surf1_cut")
file_out_merge = os.path.join(path_output, "surf_merge")
file_out_merge2 = os.path.join(path_output, "surf_merge2")
file_out_merge3 = os.path.join(path_output, "surf_merge3")
file_out_merge4 = os.path.join(path_output, "surf_merge4")
file_edge = os.path.join(path_output, "edge.mgh")
file_edge_sorted = os.path.join(path_output, "edge_sorted.mgh")
file_overlay = os.path.join(path_output, "overlay.mgh")

# load input
vtx0, fac0 = read_geometry(file_in)
n = Mesh(vtx0, fac0).vertex_normals
nlayer = len(file_contrast)

if file_in2:
    vtx1, _ = read_geometry(file_in2)
else:
    # shift vertices
    vtx1 = vtx0 + s * n
    vtx0 = vtx0 - s * n

# get header for overlay
affine = nb.load(file_contrast[0]).affine
header = nb.load(file_contrast[0]).header

# define label to define vertices which will be removed
label = cut_mesh(
    vtx0,
    x_min=x_min,
    x_max=x_max,
    y_min=y_min,
    y_max=y_max,
    z_min=z_min,
    z_max=z_max,
    merge=merge,
)

if file_label:
    label_tmp = read_label(file_label)
    if merge == "union":
        label = list(set(label).union(set(label_tmp)))
    elif merge == "intersection":
        label = list(set(label).intersection(set(label_tmp)))
    else:
        sys.exit("No valid merge parameter")

label = np.unique(label)
label = np.sort(label)[::-1]

ind_keep = np.arange(0, len(vtx0))
ind_keep[label] = -1
ind_keep = ind_keep[ind_keep != -1]

# smooth surfaces
write_geometry(file_out0, vtx0, fac0)
write_geometry(file_out1, vtx1, fac0)
mris_smooth(file_out0, file_out0, 5)
mris_smooth(file_out1, file_out1, 5)

# reload smooth surfaces
vtx0, fac0 = read_geometry(file_out0)
vtx1, _ = read_geometry(file_out1)

# remove vertices
vtx0_new = vtx0[ind_keep, :]
vtx1_new = vtx1[ind_keep, :]

# sort faces
fac_new = remove_face(fac0, ind_keep, label)

write_geometry(file_out0_cut, vtx0_new, fac_new)
write_geometry(file_out1_cut, vtx1_new, fac_new)

# get edge vertices
ind_edge = get_edge(ind_keep, vtx0, fac0)

# label edges in overlay
arr_edge = np.zeros_like(vtx0_new[:, 0])
arr_edge[ind_edge] = 1

# write output
arr_edge = np.expand_dims(arr_edge, axis=1)
arr_edge = np.expand_dims(arr_edge, axis=1)
output = nb.Nifti1Image(arr_edge, affine, header)
nb.save(output, file_edge)

# sort edge indices
ind_edge_sorted = [ind_edge[0]]
ind_edge_rest = ind_edge.copy()
ind_edge_rest.remove(ind_edge[0])
nn_2d_found = True
hole = []
c = 0
while True:
    # find nearest neighbors of last sorted edge vertex. If no vertices are
    # found, edges are searched within a sphere. This is necessary, since the
    # edge boundary might not be fully connected. Thus, edge indices are sorted
    # here but neighboring indices might be far away and do not share a
    # boundary. This has to be accounted for in the cortical depth filling
    # procedure.
    if nn_2d_found:
        nn = Mesh(vtx0_new, fac_new).neighborhood(ind_edge_sorted[-1])
    else:
        hole.append(c)
        nn, _ = nn_3d(vtx0_new[ind_edge_sorted[-1]], vtx0_new, 100)
        nn_2d_found = True

    # get all neighboring vertices which are on the edge of the surface mesh and
    # are not already stored.
    nn_edge = list(set(nn).intersection(set(ind_edge)))
    nn_edge = list(set(nn_edge) - set(ind_edge_sorted))

    # get neighbor numbers from all left neighbor vertices. This is done because
    # the neighbor vertex with least own neighbors will be finally stored in the
    # sorted edge array. This ensures that no vertex will be missed when we move
    # along the edge.
    tmp = [
        len(Mesh(vtx0_new, fac_new).neighborhood(nn_edge[i]))
        for i in range(len(nn_edge))
    ]

    # if no neighbor is found, search again within the sphere.
    if not nn_edge:
        nn_2d_found = False
        continue

    c += 1

    # store chosen index and remove that index from the edge array.
    ind_edge_sorted.append(nn_edge[np.argmin(tmp)])
    ind_edge_rest.remove(ind_edge_sorted[-1])

    if not ind_edge_rest:
        break

# sample sorted edges
arr_sorted = nb.load(file_contrast[0]).get_fdata()
arr_sorted = arr_sorted[ind_keep]
arr_sorted[:] = 0
arr_sorted[ind_edge_sorted, 0, 0] = np.arange(1, len(ind_edge_sorted) + 1)
output = nb.Nifti1Image(arr_sorted, affine, header)
nb.save(output, file_edge_sorted)

# merge surfaces
vtx_merge = np.concatenate((vtx0_new, vtx1_new), axis=0)
fac_merge = np.concatenate((np.flip(fac_new, axis=1), fac_new + len(vtx0_new)), axis=0)

write_geometry(file_out_merge, vtx_merge, fac_merge)

# add cortical depth-dependent coordinates
ind_array = np.zeros((nlayer + 1, len(ind_edge_sorted))).astype(int)
ind_array[0, :] = ind_edge_sorted
ind_array[-1, :] = ind_edge_sorted
ind_array[-1, :] += len(vtx0_new)

j = 0
while j < nlayer:
    for i in range(np.shape(ind_array)[1]):
        p0 = vtx0_new[ind_edge_sorted[i]]
        p1 = vtx1_new[ind_edge_sorted[i]]

        line_curr = np.linspace(
            (0, 0, 0),
            (
                1,
                1,
                1,
            ),
            nlayer + 1,
            dtype=np.float,
        )
        line_curr = (p1 - p0) * line_curr + p0

        vtx_merge = np.vstack((vtx_merge, line_curr[j + 1, :]))

        if j != nlayer - 1:
            ind_array[j + 1, i] = len(vtx_merge) - 1

    j += 1

    # current status
    print("add indices for layer: " + str(j) + " of " + str(nlayer))

# add cortical depth-dependent faces
ind_edge_sorted_nn = []
for i in range(np.shape(ind_array)[0] - 1):
    # current status
    print("add faces for layer: " + str(i + 1) + " of " + str(nlayer))

    for j in range(np.shape(ind_array)[1] - 1):
        if j in hole:
            continue

        if i == 0:
            ind_edge_sorted_nn.append(ind_array[0, j + 1])

        a = ind_array[i, j]
        b = ind_array[i + 1, j]
        c = ind_array[i, j + 1]
        d = ind_array[i + 1, j + 1]

        fac_merge = np.vstack((fac_merge, [a, b, c]))
        fac_merge = np.vstack((fac_merge, [b, d, c]))

write_geometry(file_out_merge2, vtx_merge, fac_merge)

# final edge sorting
for i, _ in enumerate(hole):
    # Neighbors in a sphere around the current hole vertex are found. The hole
    # index itself is then removed from the array. Furthermore, only vertices
    # which are at the edge are kept.
    indd = ind_edge_sorted[hole[i]]
    nn, _ = nn_3d(vtx0_new[indd], vtx0_new, 10)
    nn = nn[nn != indd]
    nn = list(set(nn).intersection(ind_edge_sorted))

    # The euclidean distances between hole vertex and all kept neighbors are
    # calculated. The index of the neighbor with minimum distance is found. The
    # found vertex is then classified as right neighbor if the hole vertex is
    # within the first order neighborhood of the putative neighbor vertex on the
    # surface.
    while True:
        tmp = norm(vtx0_new[indd] - vtx0_new[nn], axis=1)
        index = np.where(tmp == np.min(tmp))[0][0]
        neighbor = np.where(ind_array[0, :] == nn[index])[0][0]
        nn_check = Mesh(vtx_merge, fac_merge).neighborhood(
            ind_array[int(nlayer / 2), neighbor]
        )

        if ind_array[int(nlayer / 2), hole[i]] in nn_check:
            nn.remove(nn[index])
        else:
            break

    # add face to fill the hole.
    for j in range(np.shape(ind_array)[0] - 1):
        a = ind_array[j, hole[i]]
        b = ind_array[j + 1, hole[i]]
        c = ind_array[j, neighbor]
        d = ind_array[j + 1, neighbor]

        tri1 = [a, b, c]
        tri2 = [b, d, c]

        fac_merge = np.vstack((fac_merge, [a, b, c]))
        fac_merge = np.vstack((fac_merge, [b, d, c]))

# write and smooth final surface
write_geometry(file_out_merge3, vtx_merge, fac_merge)
mris_smooth(file_out_merge3, file_out_merge4, n_smooth)

# new overlay
arr_overlay = np.zeros(len(vtx0_new))  # inner side
arr_overlay = np.concatenate(
    (arr_overlay, nb.load(file_contrast[int(nlayer / 2)]).get_fdata()[ind_keep, 0, 0])
)  # outer side

# add cortical depth-dependent overlay values
for i in range(nlayer):
    tmp = nb.load(file_contrast[i]).get_fdata()[ind_keep, 0, 0]
    arr_overlay = np.concatenate((arr_overlay, tmp[ind_array[0, :]]))

# add empty dimensions
arr_overlay = np.expand_dims(arr_overlay, axis=1)
arr_overlay = np.expand_dims(arr_overlay, axis=1)

# write output
output = nb.Nifti1Image(arr_overlay, affine, header)
nb.save(output, file_overlay)
