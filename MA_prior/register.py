#!/usr/bin/env python
# 
# Script for segmenting cardiac images using multi-atlas label fusion.
# Note that these parameters may need to be adjusted for specific applications.

import sys
import imp
import os
#import nibabel as nib
import subprocess

train1 =  ['11', '14', '18', '20', '25','32', '35', '41', '45', '49', '54', '58', '12', '15', '19E', '23', '28','33', '38', '42', '47', '51', '55', '59','13', '16', '19', '24', '31','34', '40', '43', '48', '52', '56']

path = "/malibo_C19/"
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

    k = sys.argv[1] #testing subject
      
    # Atlas images and label maps that are already warped or registered to the target image space.
    # Image registration can be performed using either IRTK or other registration tools (SyN, ITK etc).
#    n_target = 23
    #n_atlas = 1

#    for i in range(0,n_target):
    for i in range (0, 48, 2):
        l = listing[i]   #training
#        l = train1[i]
        print l, k

        subprocess.call("MIRTK_bin/bin/mirtk register malibo_C19/"+ str(l) + "/vol.nii.gz malibo_C19/" + str(k) +"/vol.nii.gz -parin ffd.cfg -dofout ffd_" + str(l) + "_to_atlas_" + str(k)+ ".dof.gz", shell=True)


if __name__ == "__main__":
    main()
