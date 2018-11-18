from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import SimpleITK as sitk
import tensorflow as tf
import os
import numpy as np

from dltk.io.augmentation import add_gaussian_noise, flip, extract_class_balanced_example_array
from dltk.io.preprocessing import whitening, normalise_one_one

def read_fn(file_references, mode, params=None):
    """A custom python read function for interfacing with nii image files.

    Args:
        file_references (list): A list of lists containing file references, such
            as [['id_0', 'image_filename_0', target_value_0], ...,
            ['id_N', 'image_filename_N', target_value_N]].
        mode (str): One of the tf.estimator.ModeKeys strings: TRAIN, EVAL or
            PREDICT.
        params (dict, optional): A dictionary to parameterise read_fn ouputs
            (e.g. reader_params = {'n_examples': 10, 'example_size':
            [64, 64, 64], 'extract_examples': True}, etc.).

    Yields:
        dict: A dictionary of reader outputs for dltk.io.abstract_reader.
    """

    def _augment(img, lbl):
        """An image augmentation function"""
        img = add_gaussian_noise(img, sigma=0.1)
        [img, lbl] = flip([img, lbl], axis=1)

        return img, lbl

    for f in file_references:
        subject_id = f[0]
        img_fn = f[1]

        # Read the image nii with sitk and keep the pointer to the sitk.Image
        # of an input

        t1_sitk = sitk.ReadImage(str(os.path.join(img_fn, 'vol.nii.gz')))
        t1 = sitk.GetArrayFromImage(t1_sitk)
        t2 = sitk.GetArrayFromImage(
            sitk.ReadImage(str(os.path.join(img_fn, 'WB_prob_11.nii.gz'))))
        
        # Normalise volume images        
        t1 = whitening(t1)
        t2 = whitening(t2)
        
        t1 = normalise_one_one(t1)
        t2 = normalise_one_one(t2)

        # Create a 4D multi-sequence image (i.e. [channels, x, y, z])
        images = np.stack([t1,  t2], axis=-1).astype(np.float32)
        

        if mode == tf.estimator.ModeKeys.PREDICT:
            yield {'features': {'x': images},
                   'labels': None,
                   'sitk': t1_sitk,
                   'subject_id': subject_id}

        lbl = sitk.GetArrayFromImage(sitk.ReadImage(str(os.path.join(
            img_fn,
            'GT11.nii.gz')))).astype(np.int32)
       
        
        # Augment if used in training mode
        if mode == tf.estimator.ModeKeys.TRAIN:
            images, lbl = _augment(images, lbl)
       
        # Check if the reader is supposed to return training examples or full
        #  images
        if len(np.argwhere(lbl == 1) != 0):
            if params['extract_examples']:
                n_examples = params['n_examples']
                example_size = params['example_size']
                images, lbl = extract_class_balanced_example_array(
                    image=images,
                    label=lbl,
                    example_size=example_size,
                    n_examples=n_examples,
                    classes=2,
                    class_weights=[1, 14]) #set foreground weight high

                for e in range(n_examples):

                    yield {'features': {'x': images[e].astype(np.float32)},
                           'labels': {'y': lbl[e].astype(np.int32)},
                           'subject_id': subject_id}
            else:
                yield {'features': {'x': images},
                       'labels': {'y': lbl},
                       'sitk': t1_sitk,
                       'subject_id': subject_id}

    return
