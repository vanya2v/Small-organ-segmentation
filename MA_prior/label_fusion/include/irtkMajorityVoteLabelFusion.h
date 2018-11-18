/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMAJORITYVOTELABELFUSION_H
#define _IRTKMAJORITYVOTELABELFUSION_H

#include "irtkLabelFusionBase.h"

/*
 * Class for multi-atlas majority vote.
 */

template<typename VoxelType> class irtkMajorityVoteLabelFusion : public irtkLabelFusionBase<VoxelType>
{
public:
  // Constructor
  irtkMajorityVoteLabelFusion();

  // Destructor
  virtual ~irtkMajorityVoteLabelFusion();

  // Name of class
  virtual const char *NameOfClass();

  // Initialise
  void Initialize();

  // Finalise
  void Finalize();

  // Run the algorithm, connecting all the parts
  void Run();

protected:
};

#include "../src/irtkMajorityVoteLabelFusion.txx"

#endif
