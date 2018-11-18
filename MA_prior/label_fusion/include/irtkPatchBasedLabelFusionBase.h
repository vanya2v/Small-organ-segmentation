/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPATCHBASEDLABELFUSIONBASE_H
#define _IRTKPATCHBASEDLABELFUSIONBASE_H

#include "irtkLabelFusionBase.h"

/*
 * Base class for patch-based label fusion.
 */

template<typename VoxelType> class irtkPatchBasedLabelFusionBase : public irtkLabelFusionBase<VoxelType>
{
public:
  // Constructor
  irtkPatchBasedLabelFusionBase();

  // Destructor
  virtual ~irtkPatchBasedLabelFusionBase();

  // Name of class
  virtual const char *NameOfClass();

  // Initialise
  virtual void Initialize();

  // Run the algorithm, connecting all the parts
  virtual void Run();

  // Check the voxel index
  // If it is outside the image boundary, duplicate the boundary voxel.
  inline void DuplicateBoundaryVoxel(int &x, int &y, int &z);

  // Parse parameters
  virtual bool Read(char *buffer1, char *buffer2);

protected:
  // Patch size
  int _patch_size_x;
  int _patch_size_y;
  int _patch_size_z;
  int _px;
  int _py;
  int _pz;
  int _n_patch;
  int *_offset_patch;
  int _offset_x;
  int _offset_y;
  int _offset_z;

  // Search volume
  int _search_vol_x;
  int _search_vol_y;
  int _search_vol_z;
  int _sx;
  int _sy;
  int _sz;
  int _n_search;
  int *_offset_search;

  // _best_patch_mode
  // true: use one best matched patch per atlas
  // false: use a lot of patches per atlas
  bool _best_patch_mode;

  // _multipoint_mode
  // true: fuse labels from the whole atlas patches
  // false: only fuse labels from atlas patch centres
  bool _multipoint_mode;

  // Mask of boundary
  irtkByteImage *_bound_check;
  irtkBytePixel *_bound_ptr;
};

#include "../src/irtkPatchBasedLabelFusionBase.txx"

#endif
