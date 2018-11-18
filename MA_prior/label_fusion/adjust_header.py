# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:34:38 2016

@author: vvv214
"""

import nibabel as nib

nim1 = nib.load('set2/multi_atlas/atlas_image_0.nii.gz')
nim2 = nib.load('set2/multi_atlas/atlas_label_0.nii.gz')
nim3 = nib.Nifti1Image(nim2.get_data(), nim1.get_affine())
nib.save(nim3, 'set2/multi_atlas/wrap_label_0.nii.gz')