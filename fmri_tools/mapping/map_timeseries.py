# choose interpolation method
# label optional?!?

# external inputs
import numpy as np
import nibabel as nb
from fmri_tools.surface.mesh import Mesh
from fmri_tools.io.affine import read_vox2ras_tkr
from fmri_tools.utils.apply_affine_chunked import apply_affine_chunked
from nibabel.freesurfer.io import read_label
from fmri_tools.io.surf import write_mgh
from fmri_tools.utils.interpolation import linear_interpolation3d, nn_interpolation3d




if __name__ == "__main__":
    file_surf = "/data/pt_01880/zdata2/lh.white_def1"
    file_cmap = "/data/pt_01880/zdata2/epi2ana.nii.gz"
    file_timeseries = "/data/pt_01880/zdata2/udata.nii"
    file_label = "/data/pt_01880/zdata2/lh.V1_exvivo.label"

    mesh = Mesh.from_file(file_surf)
    vtx_res = mesh.transform_coords(file_cmap, file_timeseries)

    #from nibabel.freesurfer import write_geometry
    #write_geometry("/data/pt_01880/zzz", vtx_res, mesh.faces)




    #nx, ny, nz, nt = np.shape(arr)

    #mask[vtx_vox[:, 0] < 0] = -1
    #mask[vtx_vox[:, 1] < 0] = -1
    #mask[vtx_vox[:, 2] < 0] = -1
    #mask[vtx_vox[:, 0] > nx] = -1
    #mask[vtx_vox[:, 1] > ny] = -1
    #mask[vtx_vox[:, 2] > nz] = -1
    #mask = mask[mask != -1]

    #vtx_vox = vtx_vox[mask, :]
    #res = nn_interpolation3d(vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2], arr[:, :, :, 0])
    #res_overlay = np.zeros(len(mesh.verts))
    #res_overlay[mask] = res

    #mask2 = np.ones(len(mesh.verts), dtype=int)
    #mask2[label] = 0
    #res_overlay[mask2 == 1] = 0

    #write_mgh("/data/pt_01880/zzz2.mgh", res_overlay)
