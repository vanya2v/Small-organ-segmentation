/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONTEXTPATCHBASEDLABELFUSION_H
#define _IRTKCONTEXTPATCHBASEDLABELFUSION_H

#include "irtkPatchBasedLabelFusionBase.h"
#include <irtkGradientImageFilter.h>
#include "irtkImageMedianFilter.h"
#include <math.h>
#include <vector>

#define PI 3.1415926535

/*
 * Class for context patch-based label fusion.
 */
template<typename VoxelType> class irtkContextPatchBasedLabelFusion : public irtkPatchBasedLabelFusionBase<VoxelType>
{
public:
  // Constructor
  irtkContextPatchBasedLabelFusion();

  // Destructor
  virtual ~irtkContextPatchBasedLabelFusion();

  // Name of class
  virtual const char *NameOfClass();

  // Initialise
  void Initialize();

  // Run the algorithm, connecting all the parts
  void Run();

  // Compute the features
  // (x, y, z): the origin of the current coordinate system, i.e. target voxel index
  // (x2, y2, z2): the centre of a moving patch
  void ComputeFeature(irtkGenericImage<VoxelType> *image, irtkRealImage *prob, irtkRealImage *grad, irtkGenericImage<VoxelType> *median, int x, int y, int z, int x2, int y2, int z2, double *feature);

  // Compute intensity feature
  inline void ComputeIntensityFeature(irtkGenericImage<VoxelType> *image, int x, int y, int z, double *feature);

  // Compute gradient feature
  inline void ComputeGradientFeature(irtkRealImage *grad, int x, int y, int z, double *feature);

  // Compute contextual feature
  inline void ComputeContextualFeature(irtkGenericImage<VoxelType> *image, irtkGenericImage<VoxelType> *median, int x, int y, int z, double *feature);

  // Compute binary contextual feature
  inline void ComputeBinaryContextualFeature(irtkGenericImage<VoxelType> *image, irtkGenericImage<VoxelType> *median, int x, int y, int z, double *feature);

  // Compute probability contextual feature
  inline void ComputeProbabilityContextualFeature(irtkRealImage *prob, int x, int y, int z, double *feature);

  // Scale each element of the feature vector to the range [-1, 1]
  inline void ScaleFeature(double *feature, double *feature_max, double *feature_min);

  // Parse parameters
  bool Read(char *buffer1, char *buffer2);

protected:
  // Gradient images
  irtkRealImage *_image_grad;
  irtkRealImage **_atlas_grad;

  // Median images
  irtkGenericImage<VoxelType> *_image_median;
  irtkGenericImage<VoxelType> **_atlas_median;

  // Probability maps
  irtkRealImage **_atlas_prob;

  // Features for label fusion
  bool _use_spatial_feature;
  bool _use_intensity_feature;
  bool _use_gradient_feature;
  bool _3D_context;
  bool _use_context_feature;
  bool _use_binary_context_feature;
  bool _use_prob_context_feature;

  int _n_ray;
  int _n_sample_per_ray;
  int _n_sample;
  int _start_radius;
  int *_offset_context_x;
  int *_offset_context_y;
  int *_offset_context_z;

  int _n_feature;
  int _n_spatial_feature;
  int _n_intensity_feature;
  int _n_gradient_feature;
  int _n_context_feature;
  int _n_binary_context_feature;
  int _n_prob_context_feature;
  bool *_require_scaling;

  // Parameters
  double _h;
  int _K;
};

#include "../src/irtkContextPatchBasedLabelFusion.txx"

#endif
