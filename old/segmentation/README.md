### HOWTO: manual white surface correction
- copy old `wm.mgz` in freesurfer trash folder and name like `wm_backup_201812061015.mgz`
- copy old white surfaces the same way
- edit `wm.mgz` (brush value: 255, eraser value: 1)
- best to overlay with orig and surfaces

### HOWTO: manual pial surface correction
- move old `brain.finalsurfs.mgz` to trash folder as in part 2 (remove from `mri` folder!)
- copy old pial and white surface the same way
- copy `orig.mgz` and save as `pial_edit.mgz` in the same folder
- `mri_convert pial_edit.mgz pial_edit.mgz -odt float`
- apply changes (brush value: 256, eraser value: -1)
- best to overlay with surfaces
- N.B. in freesurfer, manual changes of the pial surfaces are applied by creating the file `brain.finalsurfs.manedit.mgz` (copy of `brainmask.mgz`, brush value: 255, eraser: 1). This file will be created when the script is run.

### HOWTO: defining a patch for surface flattening
- open `tksurfer` and define manually the patch of the occipital pole
- if you define the patch onto the upsampled surface mesh, you have to open first the original surface in `tksurfer` and then load the upsampled surface in the opened GUI
- load the inflated surface
- rotate to the medial surface
- select points along the calcarine fissure and press the button "Cut line"
- select 3 points to define the cutting plane: 2 on medial side and 1 on lateral side
- choose a 4th points to specify which portion of surface to keep and press button "Cut plane"
- save file: File > Patch > Save as file `<hemi>.<namePATCH>.patch.3d`
- save the file in the dense subfolder
- after flattening you can visualise the patch by loading first the inflated surface in tksurfer, then File > Patch > Load patch ...