/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONTEXTSVMLABELFUSION_H
#define _IRTKCONTEXTSVMLABELFUSION_H

#include "irtkContextPatchBasedLabelFusion.h"
#include "svm.h"

struct SVM {
  svm_model *model;
  svm_node *SV_memory;
  int label;
  double *feature_max;
  double *feature_min;
};

/*
 * Class for context SVM label fusion.
 */
template<typename VoxelType> class irtkContextSVMLabelFusion : public irtkContextPatchBasedLabelFusion<VoxelType>
{
public:
  // Constructor
  irtkContextSVMLabelFusion();

  // Destructor
  virtual ~irtkContextSVMLabelFusion();

  // Name of class
  virtual const char *NameOfClass();

  // Initialise
  void Initialize();

  // Run the algorithm, connecting all the parts
  void Run();

  // Train a SVM classifier at voxel (x, y, z)
  inline void TrainSVMClassifier(SVM &node, int x, int y, int z, int i);

  // Apply the SVM at control point (si, sj, sk) to the patch at voxel (x2, y2, z2)
  inline int ApplySVMClassifier(const SVM &node, int si, int sj, int sk, int x2, int y2, int z2);

  // Parse parameters
  bool Read(char *buffer1, char *buffer2);

protected:
  // Offset and weight tables in the SVM support region
  int _n_support;
  int *_offset_support;
  double *_weight_support;

  // SVM control point spacing
  int _svm_dx;
  int _svm_dy;
  int _svm_dz;

  // The number of SVM control points along each dimension
  int _svm_X;
  int _svm_Y;
  int _svm_Z;

  // Parameters
  double _gamma;
  double _C;
};

#include "../src/irtkContextSVMLabelFusion.txx"

#endif
