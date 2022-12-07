# load mesh
# load cmap
# load timeseries


# external inputs
import numpy as np
import nibabel as nb
from fmri_tools.surface.mesh import Mesh
















def linear_interpolation3d(x, y, z, arr_c):
    """Apply a linear interpolation of values in a 3D volume to an array of
    coordinates.
    Parameters
    ----------
    x : (N,) np.ndarray
        x-coordinates in voxel space.
    y : (N,) np.ndarray
        y-coordinates in voxel space.
    z : (N,) np.ndarray
        z-coordinates in voxel space.
    arr_c : (U,V,W) np.ndarray
        3D array with input values.
    Returns
    -------
    c : (N,) np.ndarray
        Interpolated values for [x,y,z].
    """

    # corner points
    x0 = np.floor(x).astype(int)
    x1 = np.ceil(x).astype(int)
    y0 = np.floor(y).astype(int)
    y1 = np.ceil(y).astype(int)
    z0 = np.floor(z).astype(int)
    z1 = np.ceil(z).astype(int)

    # distances to corner points
    xd = [_careful_divide(x[i], x0[i], x1[i]) for i, _ in enumerate(x)]
    yd = [_careful_divide(y[i], y0[i], y1[i]) for i, _ in enumerate(y)]
    zd = [_careful_divide(z[i], z0[i], z1[i]) for i, _ in enumerate(z)]

    xd = np.asarray(xd)
    yd = np.asarray(yd)
    zd = np.asarray(zd)

    # corner values
    c000 = arr_c[x0, y0, z0]
    c001 = arr_c[x0, y0, z1]
    c010 = arr_c[x0, y1, z0]
    c011 = arr_c[x0, y1, z1]
    c100 = arr_c[x1, y0, z0]
    c101 = arr_c[x1, y0, z1]
    c110 = arr_c[x1, y1, z0]
    c111 = arr_c[x1, y1, z1]

    # interpolation along x-axis
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    # interpolation along y-axis
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    # interpolation along z-axis
    c = c0 * (1 - zd) + c1 * zd

    return c


def nn_interpolation3d(x, y, z, arr_c):
    """Apply a nearest neighbor interpolation of values in a 3D volume to an
    array of coordinates.
    Parameters
    ----------
    x : (N,) np.ndarray
        x-coordinates in voxel space.
    y : (N,) np.ndarray
        y-coordinates in voxel space.
    z : (N,) np.ndarray
        z-coordinates in voxel space.
    arr_c : (U,V,W) np.ndarray
        3D array with input values.
    Returns
    -------
    c : (N,) np.ndarray
        Interpolated values for [x,y,z].
    """

    # get nearest neighbour grid points
    x0 = np.round(x).astype(int)
    y0 = np.round(y).astype(int)
    z0 = np.round(z).astype(int)

    c = arr_c[x0, y0, z0]

    return c


def _careful_divide(v, v0, v1):
    """Only divide if v0 and v1 are different from each other."""

    return (v - v0) / (v1 - v0) if v1 != v0 else v


if __name__ == "__main__":
    # mache ein bisschen schoener
    # choose interpolation method
    # cmap optional
    # label optional

    from fmri_tools.io.affine import read_vox2ras_tkr
    from fmri_tools.utils.apply_affine_chunked import apply_affine_chunked
    from nibabel.freesurfer.io import read_label
    from fmri_tools.io.surf import write_mgh

    file_surf = "/data/pt_01880/zdata2/lh.white_def1"
    file_cmap = "/data/pt_01880/zdata2/epi2ana.nii.gz"
    file_timeseries = "/data/pt_01880/zdata2/udata.nii"
    file_label = "/data/pt_01880/zdata2/lh.V1_exvivo.label"

    mesh = Mesh.from_file(file_surf)
    cmap = nb.load(file_cmap)
    arr_cmap = cmap.get_fdata()
    nx, ny, nz, _ = np.shape(arr_cmap)
    data = nb.load(file_timeseries)
    arr = data.get_fdata()
    label = read_label(file_label)

    _, ras2vox = read_vox2ras_tkr(file_cmap)
    vox2ras, _ = read_vox2ras_tkr(file_timeseries)
    vtx_vox = apply_affine_chunked(ras2vox, mesh.verts)

    mask = np.arange(len(mesh.verts))
    mask[vtx_vox[:, 0] < 0] = -1
    mask[vtx_vox[:, 1] < 0] = -1
    mask[vtx_vox[:, 2] < 0] = -1
    mask[vtx_vox[:, 0] > nx] = -1
    mask[vtx_vox[:, 1] > ny] = -1
    mask[vtx_vox[:, 2] > nz] = -1
    mask = mask[mask != -1]

    vtx_vox = vtx_vox[mask, :]
    x = linear_interpolation3d(vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2], arr_cmap[:,:,:,0])
    y = linear_interpolation3d(vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2], arr_cmap[:,:,:,1])
    z = linear_interpolation3d(vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2], arr_cmap[:,:,:,2])

    vtx_vox = np.array([x, y, z]).T
    nx, ny, nz, nt = np.shape(arr)

    mask[vtx_vox[:, 0] < 0] = -1
    mask[vtx_vox[:, 1] < 0] = -1
    mask[vtx_vox[:, 2] < 0] = -1
    mask[vtx_vox[:, 0] > nx] = -1
    mask[vtx_vox[:, 1] > ny] = -1
    mask[vtx_vox[:, 2] > nz] = -1
    mask = mask[mask != -1]

    #vtx_res = apply_affine_chunked(vox2ras, vtx_res)
    #from nibabel.freesurfer import write_geometry
    #write_geometry("/data/pt_01880/zzz", vtx_res, mesh.faces)

    vtx_vox = vtx_vox[mask, :]
    res = nn_interpolation3d(vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2], arr[:, :, :, 0])
    res_overlay = np.zeros(len(mesh.verts))
    res_overlay[mask] = res



    mask2 = np.ones(len(mesh.verts), dtype=int)
    mask2[label] = 0
    res_overlay[mask2 == 1] = 0

    write_mgh("/data/pt_01880/zzz2.mgh", res_overlay)
