/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _IRTKMAJORITYVOTELABELFUSION_TXX
#define _IRTKMAJORITYVOTELABELFUSION_TXX

template<typename VoxelType> irtkMajorityVoteLabelFusion<VoxelType>::irtkMajorityVoteLabelFusion()
{
}

template<typename VoxelType> irtkMajorityVoteLabelFusion<VoxelType>::~irtkMajorityVoteLabelFusion()
{
}

template<typename VoxelType> const char *irtkMajorityVoteLabelFusion<VoxelType>::NameOfClass()
{
  return "irtkMajorityVoteLabelFusion";
}

template<typename VoxelType> void irtkMajorityVoteLabelFusion<VoxelType>::Initialize()
{
  this->irtkLabelFusionBase<VoxelType>::Initialize();
}

template<typename VoxelType> void irtkMajorityVoteLabelFusion<VoxelType>::Finalize()
{
  this->irtkLabelFusionBase<VoxelType>::Finalize();
}

template<typename VoxelType> void irtkMajorityVoteLabelFusion<VoxelType>::Run()
{
  double *vote = new double[this->_number_of_classes];

  for(int i=0; i<this->_number_of_voxels; i++){
    if(this->_already_ptr[i] == 1){
      continue;
    }

    // Initialise the vote vector
    for(int k=0; k<this->_number_of_classes; k++){
      vote[k] = 0;
    }

    // For each atlas
    for(int n=0; n<this->_number_of_atlases; n++){
      int l = this->_label_ptr[n][i];
      vote[this->_label_idx[l]]++;
    }
    
    // Label with the majority vote
    double max_vote = -1E10;
    int max_k = 0;
    for(int k=0; k<this->_number_of_classes; k++){
      if(vote[k] > max_vote){
	max_vote = vote[k];
	max_k = k;
      }
    }

    set<int>::iterator it = this->_label_set.begin();
    advance(it, max_k);
    this->_seg_ptr[i] = *it;

    // Probability
    for(int k=0; k<this->_number_of_classes; k++){
      this->_prob_ptr[k][i] = vote[k] / (double)this->_number_of_atlases;
    }    
  }

  delete []vote;
  vote = NULL;
}

#endif
