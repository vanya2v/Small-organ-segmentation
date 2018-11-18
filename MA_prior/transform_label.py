#!/usr/bin/env python
# 
# Script for segmenting cardiac images using multi-atlas label fusion.
# Note that these parameters may need to be adjusted for specific applications.


import sys
import imp
import os
import subprocess

sys.path.append('/usr/lib/python2.7/dist-packages:/vol/medic02/users/vvv214/sw')
import nibabel as nib


path = "malibo_C19"
listing = sorted(os.listdir(path))

def find_module(modulename, filename=None):
    """Finds a python module or package on the standard path.
       If a filename is specified, add its containing folder
       to the system path.
 
       Returns a string of the full path to the module/package."""

 
    full_path = []
    if filename:
        full_path.append(os.path.dirname(os.path.abspath(filename)))
    full_path += sys.path
    fname = imp.find_module(modulename, full_path)
    return fname[1]



def main():

    k = sys.argv[1] # testing
    
#    k = listing[int(n)]
      
    # Atlas images and label maps that are already warped or registered to the target image space.
    # Image registration can be performed using either IRTK or other registration tools (SyN, ITK etc).
#    n_target = 24
    #n_atlas = 1

    for i in range (0, 48, 2):
        l = listing[i] #training

        print  k,l
        
        nim1 = nib.load("/malibo_C19/" + str(k) +"/vol.nii.gz") 
        nim2 = nib.load("/malibo_C19/" + str(k) +"/label.nii.gz") 
        nim3 = nib.Nifti1Image(nim2.get_data(), nim1.get_affine())
        
        nib.save(nim3, "/MA_Prior/atlas/atlas_" + str(k) + "/" + str(k) + "_label.nii.gz")
        nib.save(nim1, "/MA_Prior/atlas/atlas_" + str(k) + "/" + str(k) + "_img.nii.gz")
        
        subprocess.call("MIRTK_bin/bin/mirtk transform-image /MA_Prior/atlas/atlas_" + str(k)+ "/" + str(k)+ "_img.nii.gz /MA_Prior/atlas/atlas_" + str(k)+ "/warp_image" + str(l)+ "_to_" + str(k)+ ".nii.gz -dofin /MA_Prior/FFD_malibo/ffd_" + str(l)+ "_to_atlas_" + str(k)+ ".dof.gz -interp Linear -target /malibo_C19/" + str(l)+ "/vol.nii.gz", shell=True)
        subprocess.call("MIRTK_bin/bin/mirtk transform-image /MA_Prior/atlas/atlas_" + str(k)+ "/" + str(k)+ "_label.nii.gz /MA_Prior/atlas/atlas_" + str(k)+ "/warp_label" + str(l)+ "_to_" + str(k)+ ".nii.gz -dofin /MA_Prior/FFD_malibo/ffd_" + str(l)+ "_to_atlas_" + str(k)+ ".dof.gz -interp NN -target /malibo_C19/" + str(l)+ "/label.nii.gz", shell=True)
        

if __name__ == "__main__":
    main()
