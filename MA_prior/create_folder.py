# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 12:01:34 2016

@author: vvv214
"""
import os
path = "/vol/bitbucket/vvv214/malibo_C19/"
listing = sorted(os.listdir(path))


for i in range (0, 48):
    k = listing[i]

#newpath = r'C:\Program Files\arbitrary' atlas_" + str(k)
    newpath = r"MA_prior/atlas/atlas_" + str(k)
    if not os.path.exists(newpath):
        os.makedirs(newpath)
