function roi_mergeMasks(path_to_mask_1,path_to_mask_2, newName)
%
% roi_mergeMasks
% Merges two masks provided. Code from bidspm FAQ question: 
% "SPM: How do I merge 2 masks with SPM?"
% https://bidspm.readthedocs.io/en/latest/FAQ.html#spm-how-do-i-merge-2-masks-with-spm

% get header of Nifti images
header_1 = spm_vol(path_to_mask_1);
header_2 = spm_vol(path_to_mask_2);

% if you want to make sure that images are in the same space
% and have same resolution
masks = char({path_to_mask_1; path_to_mask_2});
spm_check_orientations(spm_vol(masks));

% get data of Nifti images
mask_1 = spm_read_vols(header_1);
mask_2 = spm_read_vols(header_2);

% concatenate data along the 4th dimension
merged_mask = cat(4, mask_1, mask_2);

% keep any voxel that has some value along the 4th dimension
merged_mask = any(merged_mask, 4);

% create a new header of the final mask
merged_mask_header = header_1;
merged_mask_header.fname = [newName, '.nii'];

spm_write_vol(merged_mask_header, merged_mask);
end

