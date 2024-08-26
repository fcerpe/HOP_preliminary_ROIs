
# HOPLAB fMRI pipeline - Preliminary ROIs extraction

This repository contains the scripts used in [HOPLAB](hoplab.be) to extract ROIs from fMRI data. 
It is currently under development, much is subject to change and may not work as a out-of-the-box tool. 

### Usage

Clone the repository

Set up your desired ROI in `hop_option()`: specify the paths and all the necessary paramters for each method. 

Run `hop_preliminary_rois` and collect your ROIs in the indicated output folder


### Methods implemented

We are trying to implement the most common methods used in HOPLAB.
Methods currently implemented: 
- extracttion of areas based on Brodmann atlas
- extraction of sphere from center coordinates
- intersection of sphere and localizer mask

Methods to implement:
- expansion within localizer mask from specific coordinates 
- extraction from different atlases

If you develop a method that you think it's useful for others, you can open a pull request to add it to the main repository.


### Workflow

- `hop_preliminary_rois`: 

  * loads the options specified for each desired ROI 
  * calls the right function to extract them

- `hop_option()`:

  Returns:
  * `opt` - default options relative to paths of the data structure, information over the experiment (MRI space, task name, default subjects to process) 
  * `roiMethod` - cell array containing information for each method requested

  Methods information
  1. Sphere: size of the sphere in radius (or radii); 
             nifti file for reference space;
             information over the area (name to assign in filename, coordinates for left and right hemisphere, information to merge or not the hemisphere into one ROI)

  2. Atlas: which atlas to use (currently available: Brodmann), or indication of custom atlas; 
            information over the area (name to assign to the ROI, label of the parcels to join). If the atlas is custom, provide a mask path instead of the parcel list

  3. Intersection: minimum voxel size 
                   information over the area (name to assign, path of the sphere to intersect, task and contrast information to find the correpsonding contrast to intersect)

  4. Expansion: TBD


- `hop_roi_atlas()`, `hop_roi_sphere()`, `hop_roi_intersection()`, `hop_roi_expansion()`:

  Call the specific method requested and store the output in the specified folder, with additonal subfolders if needed
  Details about the step-by-step workflow of each method can be found in the functions' documentations
