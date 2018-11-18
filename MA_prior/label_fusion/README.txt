The label_fusion tool implements the multi-atlas segmentation method. It
currently includes four algorithms:
(1) majority-vote (MV)
(2) patch-based segmentation (PB)
(3) patch-based segmentation with augmented features (PBAF)
(4) SVM segmentation with augmented features (SVMAF)

A statically linked executable has already been provided. It was compiled in
Ubuntu 12.04. To test the algorithms, please go to the demo_script directory
and run the demo script on the testing data.

You can also compile the code on your own. To that end, Prof. Daniel Rueckert's
IRTK needs to be downloaded and compiled at the first place. The label_fusion
tool will be included into the next release of IRTK.

Related publication:

[1] W Bai et al. A probabilistic patch-based label fusion model for multi-atlas
segmentation with registration refinement: Application to cardiac MR images.
IEEE Transactions on Medical Imaging, 32(7):1302-1315, 2013.

Wenjia Bai (w.bai@imperial.ac.uk)
Imperial College London
9 Sep 2014
