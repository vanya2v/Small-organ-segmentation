/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKLABELFUSIONBASE_H
#define _IRTKLABELFUSIONBASE_H

#include <irtkImage.h>
#include <set>
#include <map>

/*
 * Base class for multi-atlas label fusion.
 *
 * VoxelType: data type of atlas images, which can be byte, short or float.
 * The data type of the label maps is always short.
 */

template<typename VoxelType> class irtkLabelFusionBase : public irtkObject
{
public:
  // Constructor
  irtkLabelFusionBase();

  // Destructor
  virtual ~irtkLabelFusionBase();

  // Name of class
  virtual const char *NameOfClass();

  // Set image
  void SetInput(irtkGenericImage<VoxelType> *image);

  // Set atlas images and label maps
  // Label map is always of data type short.
  void SetAtlas(int n_atlas, irtkGenericImage<VoxelType> **atlas, irtkGreyImage **label_map);

  // Initialise
  virtual void Initialize();

  // Run the algorithm, connecting all the parts
  virtual void Run();

  // Return hard segmentation
  irtkGreyImage *GetOutput();

  // Return soft segmentation
  irtkRealImage *GetProbabilityMap();

  // Initial soft segmentation estimate
  void SetInitialProbabilityMap(irtkRealImage *init_prob_map);

  // Read parameter file
  void Read(char *filename);

  // Parse a parameter line
  int ReadLine(istream &in, char *buffer1, char *&buffer2);

  // Parse parameters
  virtual bool Read(char *buffer1, char *buffer2);

protected:
  // Input image
  irtkGenericImage<VoxelType> *_image;
  VoxelType *_image_ptr;
  int _X, _Y, _Z;
  int _number_of_voxels;
  double _xsize, _ysize, _zsize;

  // Multiple atlases, in the same image space and of the same dimension as the input image
  irtkGenericImage<VoxelType> **_atlas;
  irtkGreyImage **_label_map;
  VoxelType **_atlas_ptr;
  irtkGreyPixel **_label_ptr;
  set<int> _label_set;
  map<int, int> _label_idx;

  // Hard segmentation
  irtkGreyImage *_segmentation;
  irtkGreyPixel *_seg_ptr;

  // Soft segmentation: probabilistic maps for each class
  irtkRealImage *_prob_map;
  irtkRealPixel **_prob_ptr;
  irtkRealImage *_init_prob_map;

  // Mask image to record whether a voxel has already been labelled or not
  irtkByteImage *_already;
  irtkBytePixel *_already_ptr;

  // Number of atlases
  int _number_of_atlases;

  // Number of classes
  int _number_of_classes;

  // Debug mode, which outputs intermediate results
  int _debug_mode;
};

#include "../src/irtkLabelFusionBase.txx"

#endif
