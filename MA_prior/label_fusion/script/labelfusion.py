#!/usr/bin/env python
# 
# Script for segmenting cardiac images using multi-atlas label fusion.
# Note that these parameters may need to be adjusted for specific applications.


#import sys
#import imp
import os
#import subprocess
import SimpleITK as sitk
import numpy as np

#sys.path.append('/usr/lib/python2.7/dist-packages:/vol/medic02/users/vvv214/sw')
#path = "/vol/bitbucket/vvv214/malibo_C19_4mm/"
#listing = sorted(os.listdir(path))
train = ['11', '12', '13', '14', '15', '16',  '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '28', '31', '32', '33', '34', '35', '37', '38', '39', '40', '41', '42', '43','44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62']

small = [5, 6, 10, 11, 13, 14, 17, 18] #classes of small organs
n_atlas = 24

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Try patch-based segmentation with augmented features
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Parameter file
# You may need to tune the parameters here.

par_file = 'pbaf_par.txt'
f = open(par_file, 'w')
f.write('Search volume size in X (pixel)   = 5\n')
f.write('Search volume size in Y (pixel)   = 5\n')
f.write('Search volume size in Z (pixel)   = 1\n')
f.write('Patch size in X (pixel)           = 7\n')
f.write('Patch size in Y (pixel)           = 7\n')
f.write('Patch size in Z (pixel)           = 1\n')
f.write('Use spatial feature               = False\n')
f.write('Use gradient feature              = True\n')
f.write('Use context feature               = True\n')
f.write('Use binary context feature        = True\n')
f.write('K                                 = 100\n') # Number of pre-selected patches
f.write('h                                 = 16\n')  # Gaussian kernel parameter
f.write('Starting radius                   = 4\n')
f.close()


for i in range(0, 48, 2):
    n = train[i] 
# Target image
    target_image = '/malibo_C19/' + str(n) + '/vol.nii.gz'
# Atlas images and label maps that are already warped or registered to the target image space.
# Image registration can be performed using either IRTK or other registration tools (SyN, ITK etc).
    images_warped = ''
    labels_warped = ''
    
    for k in range(1, 48, 2):
            
        m = train[k]
        print(m)
        images_warped += '/MA_prior/atlas/atlas_' + str(n) + '/warp_image' + str(m) + '_to_' + str(n) + '.nii.gz '
 

        labels_warped += '/MA_prior/atlas/atlas_' + str(n) + '/warp_label' + str(m) + '_to_' + str(n) + '.nii.gz '

        seg = 'pbaf_seg_' + str(n) + '.nii.gz'
        print(seg)
        output_prob_file = 'probmap_pbaf_' + str(n) + '.nii.gz'
        os.system('time ../label_fusion {0} {1} {2} {3} {4} -method PBAF -par {5} -output_prob {6}'.format(target_image, n_atlas, images_warped, labels_warped, seg, par_file, output_prob_file))

