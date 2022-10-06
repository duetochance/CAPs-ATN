#Preliminary FLIRT registration

tmp_map = os.path.join(T1_on_DTI, 'tmp.mat')
cmd = f"flirt -v -in {b0mean_upsaml} -ref {TT_T1_brain_FS} -dof 6 -omat {tmp_map}"
os.system(cmd)

prereg_output_matrix = os.path.join(T1_on_DTI, 'diff_to_structural_bbr.mat')
schedule = os.path.join(schedule_path, 'bbr.sch')
segmentation_result = os.path.join(T1_on_DTI, 'diff_to_structural_bbr.nii.gz')

# refine the registration 
cmd = f"flirt -v -in {b0mean_upsaml} -ref {TT_T1_brain_FS} -dof 6 -cost bbr -wmseg {whitematter_fsl} -init {tmp_map} -omat {prereg_output_matrix} -out {segmentation_result} -schedule {schedule}"
os.system(cmd)

diff_to_structural_bbr_inverse = os.path.join(T1_on_DTI, 'diff_to_structural_bbr_inverse.mat')

# Invert the previous registration
cmd = f"convert_xfm -omat {diff_to_structural_bbr_inverse} -inverse {prereg_output_matrix}"
os.system(cmd)

T1_to_b0mean_FSLbbr = os.path.join(T1_on_DTI, "T1_to_b0mean_FSLbbr.nii.gz")

# Applying the inverse transformation to T1
cmd = f"flirt -v -applyxfm -init {diff_to_structural_bbr_inverse} -in {TT_T1_brain_FS} -ref {b0mean_upsaml} -out {T1_to_b0mean_FSLbbr}"
os.system(cmd)

WM_to_b0 = os.path.join(T1_on_DTI, "WM_to_b0.nii.gz")
#Applying the inverse transformation to T1_WM
cmd = f"flirt -v -applyxfm -init {diff_to_structural_bbr_inverse} -in {whitematter_fsl} -ref {b0mean_upsaml} -out {WM_to_b0}"
os.system(cmd)