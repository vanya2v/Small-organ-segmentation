##Small organ segmentation
Small organ segmentation in whole-body scans is challenging. A coarse-to-fine, hierarchical strategy is a common approach to alleviate this problem, with weighting schemes based on auto-context and spatial atlas priors to deal with anatomical variation and class imbalance. This repository is build based on DLTK for FCN baseline and multi-atlas segmentation (http://www.doc.ic.ac.uk/~wbai/software(label_fusion_v2.0.tar.gz) for spatial atlas priors.


### Referencing and citing 
If you use this multi-modal learning repository in your work please refer to this citation:

Vanya V. Valindria, Ioannis Lavdas, Juan Cerrolaza, Eric O. Aboagye, Andrea G. Rockall, Daniel Rueckert, Ben Glocker. "Small Organ Segmentation in Whole-body MRI using a Two-stage FCN and Weighting Schemes," in 9th International Conference on Machine Learning in Medical Imaging (MLMI), 2018.


## Deep Learning Toolkit (DLTK) for Medical Imaging
[![Gitter](https://badges.gitter.im/DLTK/DLTK.svg)](https://gitter.im/DLTK/DLTK?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![Coverage Status](https://coveralls.io/repos/github/DLTK/DLTK/badge.svg?branch=master)](https://coveralls.io/github/DLTK/DLTK?branch=dev)
[![Build Status](https://travis-ci.org/DLTK/DLTK.svg?branch=master)](https://travis-ci.org/DLTK/DLTK)


### Installation
1. Setup a virtual environment and activate it. If you intend to run this on 
machines with different system versions, use the --always-copy flag:

   ```shell
   virtualenv -p python3 --always-copy venv_tf
   source venv_tf/bin/activate
   ```
   
2. Install TensorFlow (>=1.4.0) (preferred: with GPU support) for your system
 as described [here](https://www.tensorflow.org/install/):
   
   ```shell
   pip install tensorflow-gpu>=1.4.0
   ```
   
3. Install DLTK:
   There are two installation options available: You can simply install dltk as is from pypi via
   
   ```shell
   pip install dltk
   ```

### Start playing
1. Coarse-scale segmentation for multi-organ segmentation (```coarse_segm```)
   First stage segmentation for multi-class segmentation with weighting schemes on FCN-based segmentation.
   ```shell
    reader_sampling.py    : set the weights per class for imbalance sampling, set each weight relative to the size/proportion of the structure in the volumes
    train_fcn_weighted.py : train the weighted FCN for coarse segmentation
    deploy_save_pmap.py   : get the coarse multi-organ segmentation results and their probability maps for each class
    ```
   
2. Multi-atlas Prior (```MA_prior```)
   Create spatial prior based on multi-organ multi-atlas segmentation.
      ```shell
    create_folder.py      : create folder for atlas target and warped images & labels
    register.py           : please refer to MRITK software (```https://github.com/BioMedIA/MIRTK```) to register images
    transform_label.py    : create warp images and labels
    script/labelfusion.py: multi-organ segmentation using PBAF label fusion strategy, and get the spatial prior (save the probabilty maps per class)
    ```

3. Fine-scale segmentation of small organ segmentation (```fine_segm```)
   Second stage segmentation for binary small organ segmentation. The cropped ROI input are the images from multiplication of auto-context (probability maps from coarse-scale) and spatial prior (probability maps from multi-atlas segmentation), combined with auto-context (probability maps of the organ from coarse-scale segmentation.

   ```shell
    reader_bladder.py     : reader for binary small organ (here the organ is bladder) segmentation. Set the weight high for foreground class.
    train_fcn_bladder.py : train the weighted FCN for bladder segmentation. Input: cropped ROI (multiplication with spatial prior) of images and probabilty maps of bladder segmentation (from coarse-scale). Set the example size small (8x8x8) as the input size now is cropped (smaller).
    deploy_bladder.py   : get the binary small-organ segmentation result
    ```

### License
See [LICENSE](https://github.com/DLTK/DLTK/blob/master/LICENSE)

