# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function

import argparse
import os
import time

import numpy as np
import pandas as pd
import tensorflow as tf
import SimpleITK as sitk

from tensorflow.contrib import predictor

from dltk.core import metrics as metrics

from dltk.utils import sliding_window_segmentation_inference

#from reader import read_fn
#from reader_pmap import read_fn
#from reader_sampling import read_fn
from reader_sampling_deploy import read_fn

READER_PARAMS = {'extract_examples': False}
N_VALIDATION_SUBJECTS = 1
DSC_all = []

def predict(args):
    # Read in the csv with the file names you would want to predict on
    file_names = pd.read_csv(
        args.csv,
        dtype=object,
        keep_default_na=False,
        na_values=[]).as_matrix()

    # We trained on the first 4 subjects, so we predict on the rest
#    file_names = file_names[-N_VALIDATION_SUBJECTS:]
#
#    print('filenames', file_names)
    file_names = file_names[0:48]
#    file_names = file_names[1:48:2]
#    print('filenames', file_names)
    #3fold, test A
#    file_names = file_names[30:47]
    # From the model_path, parse the latest saved model and restore a
    # predictor from it
    export_dir = [os.path.join(args.model_path, o) for o in os.listdir(args.model_path)
                  if os.path.isdir(os.path.join(args.model_path, o)) and
                  o.isdigit()][-1]
    print('Loading from {}'.format(export_dir))
    my_predictor = predictor.from_saved_model(export_dir)
#    for o in os.listdir(args.model_path):
#        if os.path.isdir(os.path.join(args.model_path, o)) and o.isdigit():
#            print(o)             
#    print('Loading from {}'.format(export_dir))
#    my_predictor = predictor.from_saved_model(export_dir)

    # Fetch the output probability op of the trained network
    y_prob = my_predictor._fetch_tensors['y_prob']
    num_classes = y_prob.get_shape().as_list()[-1]

    # Iterate through the files, predict on the full volumes and compute a Dice
    # coefficient
    for output in read_fn(file_references=file_names,
                          mode=tf.estimator.ModeKeys.EVAL,
                          params=READER_PARAMS):
        t0 = time.time()

        # Parse the read function output and add a dummy batch dimension as
        # required
        img = np.expand_dims(output['features']['x'], axis=0)
        lbl = np.expand_dims(output['labels']['y'], axis=0)
        
        print('Id={}'.format(output['subject_id']))

        # Do a sliding window inference with our DLTK wrapper
        prob = sliding_window_segmentation_inference(
            session=my_predictor.session,
            ops_list=[y_prob],
            sample_dict={my_predictor._feed_tensors['x']: img},
            batch_size=64)[0]
        
        newpath = r"data_fcn_weighted/" + output['subject_id']
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        for c in range(0,19): 
            output_pm = "data_fcn_weighted/" + output['subject_id'] + "/WB_prob_"+ str(c) + ".nii.gz"
#        output_pm = os.path.join(args.model_path, '{}_prob.nii.gz'.format(output['subject_id']))

            probmap = sitk.GetImageFromArray(prob[0, :,:,:, c].astype(np.float32))
            probmap.CopyInformation(output['sitk'])
#        print('size prob', np.size(prob))
#            print('size pmap', probmap.GetSize())
#            print('size pmap'+str(c), probmap.GetSize())
            sitk.WriteImage(probmap, output_pm)
                 
        # Calculate the prediction from the probabilities
        pred = np.argmax(prob, -1)

        # Calculate the Dice coefficient
        dsc = metrics.dice(pred, lbl, num_classes)[1:].mean()
        dsc1 = metrics.dice(pred, lbl, num_classes)[1]
        dsc2 = metrics.dice(pred, lbl, num_classes)[2]
        dsc3 = metrics.dice(pred, lbl, num_classes)[3]
        dsc4 = metrics.dice(pred, lbl, num_classes)[4]
        dsc5 = metrics.dice(pred, lbl, num_classes)[5]
        dsc6 = metrics.dice(pred, lbl, num_classes)[6]
        dsc7 = metrics.dice(pred, lbl, num_classes)[7]
        dsc8 = metrics.dice(pred, lbl, num_classes)[8]
        dsc9 = metrics.dice(pred, lbl, num_classes)[9]
        dsc10 = metrics.dice(pred, lbl, num_classes)[10]
        dsc11 = metrics.dice(pred, lbl, num_classes)[11]
        dsc12 = metrics.dice(pred, lbl, num_classes)[12]
        dsc13 = metrics.dice(pred, lbl, num_classes)[13]
        dsc14 = metrics.dice(pred, lbl, num_classes)[14]
        dsc15 = metrics.dice(pred, lbl, num_classes)[15]
        dsc16 = metrics.dice(pred, lbl, num_classes)[16]
        dsc17 = metrics.dice(pred, lbl, num_classes)[17]
        dsc18 = metrics.dice(pred, lbl, num_classes)[18]
        ds = [dsc1, dsc2, dsc3, dsc4, dsc5, dsc6, dsc7, dsc8, dsc9, dsc10, dsc11, dsc12, dsc13, dsc14, dsc15, dsc16, dsc17, dsc18]
         
        DSC_all.append(ds)
        print(np.shape(DSC_all))
        # Save the file as .nii.gz using the header information from the
        # original sitk image
        output_fn = os.path.join(args.model_path, '{}_seg.nii.gz'.format(output['subject_id']))

        new_sitk = sitk.GetImageFromArray(pred[0].astype(np.int32))
        new_sitk.CopyInformation(output['sitk'])

        sitk.WriteImage(new_sitk, output_fn)

        # Print outputs
#        print('Id={}; Dice={:0.4f}; time={:0.2} secs; output_path={};'.format(
#            output['subject_id'], dsc, time.time() - t0, output_fn))
        print('Id={}; Dice adrnl= {:0.4f}; glbd= {:0.4f}; panc = {:0.4f}; rectum = {:0.4f}; output_path={};'.format(
            output['subject_id'], dsc5, dsc6, dsc10, dsc17, output_fn))
    np.save(args.save_npy, DSC_all) 

if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description='MRBrainS13 example segmentation deploy script')
    parser.add_argument('--verbose', default=False, action='store_true')
    parser.add_argument('--cuda_devices', '-c', default='0')

    parser.add_argument('--model_path', '-p', default='train_fcn_weighted')
    parser.add_argument('--csv', default='malibo_C19.csv')
#    parser.add_argument('--csv', default='malibo_2mm_pmap_WB.csv')
    parser.add_argument('--save_npy', '-s', default='DSC_fcn_weighted.npy')
    
    args = parser.parse_args()

    # Set verbosity
    if args.verbose:
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
        tf.logging.set_verbosity(tf.logging.INFO)
    else:
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
        tf.logging.set_verbosity(tf.logging.ERROR)

    # GPU allocation options
    os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_devices

    # Call training
    predict(args)
