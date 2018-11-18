/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPATCHBASEDMSDLABELFUSION_H
#define _IRTKPATCHBASEDMSDLABELFUSION_H

#include "irtkPatchBasedLabelFusionBase.h"

/*
 * Class for patch-based label fusion.
 */

template<typename VoxelType> class irtkPatchBasedMSDLabelFusion : public irtkPatchBasedLabelFusionBase<VoxelType>
{
public:
  // Constructor
  irtkPatchBasedMSDLabelFusion();

  // Destructor
  virtual ~irtkPatchBasedMSDLabelFusion();

  // Name of class
  virtual const char *NameOfClass();

  // Initialise
  virtual void Initialize();

  // Run the algorithm, connecting all the parts
  virtual void Run();

  // Fuse the label from the atlas patch (n, x2, y2, z2, i2) to the target patch (x, y, z, i) with weight w
  inline void FuseAtlasLabelWithBoundCheck(int x, int y, int z, int i, int n, int x2, int y2, int z2, int i2, double w);

  // Fuse the label from the atlas patch (n, i2) to the target patch (i) with weight w
  inline void FuseAtlasLabelWithoutBoundCheck(int i, int n, int i2, double w);

  // Parse parameters
  virtual bool Read(char *buffer1, char *buffer2);

protected:
  // Weight due to distance
  double *_weight_dist;

  // Parameters
  double _h;
  double _sigma_dist_x;
  double _sigma_dist_y;
  double _sigma_dist_z;
  bool _use_dist_wt;
};

#include "../src/irtkPatchBasedMSDLabelFusion.txx"

#endif
