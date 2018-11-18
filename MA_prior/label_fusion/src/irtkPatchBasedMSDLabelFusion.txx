/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _IRTKPATCHBASEDMSDLABELFUSION_TXX
#define _IRTKPATCHBASEDMSDLABELFUSION_TXX

template<typename VoxelType> irtkPatchBasedMSDLabelFusion<VoxelType>::irtkPatchBasedMSDLabelFusion()
{
  _h = 1;
  _use_dist_wt = false;

  _weight_dist = NULL;
}

template<typename VoxelType> irtkPatchBasedMSDLabelFusion<VoxelType>::~irtkPatchBasedMSDLabelFusion()
{
  if(_weight_dist != NULL){
    delete []_weight_dist;
    _weight_dist = NULL;
  }
}

template<typename VoxelType> const char *irtkPatchBasedMSDLabelFusion<VoxelType>::NameOfClass()
{
  return "irtkPatchBasedMSDLabelFusion";
}

template<typename VoxelType> void irtkPatchBasedMSDLabelFusion<VoxelType>::Initialize()
{
  this->irtkPatchBasedLabelFusionBase<VoxelType>::Initialize();

  // Determine the sigma of distance weighting by dividing the search volume by 4,
  // i.e. search volume is equal to +/- 2 \times sigma_dist.
  _sigma_dist_x = (2 * this->_sx + 1) * this->_xsize / 4;
  _sigma_dist_y = (2 * this->_sy + 1) * this->_ysize / 4;
  _sigma_dist_z = (2 * this->_sz + 1) * this->_zsize / 4;

  cout << "Sigma distance in X = " << _sigma_dist_x << " mm." << endl;
  cout << "Sigma distance in Y = " << _sigma_dist_y << " mm." << endl;
  cout << "Sigma distance in Z = " << _sigma_dist_z << " mm." << endl;

  // Pre-compute the weight due to the distance between the target patch and the atlas patch
  if(_use_dist_wt){
    _weight_dist = new double[this->_n_search];

    int j = 0;
    for(int dz=-this->_sz; dz<=this->_sz; dz++){
      for(int dy=-this->_sy; dy<=this->_sy; dy++){
	for(int dx=-this->_sx; dx<=this->_sx; dx++, j++){
	  double dist_x = dx * this->_xsize;
	  double dist_y = dy * this->_ysize;
	  double dist_z = dz * this->_zsize;

	  _weight_dist[j] = exp(- (dist_x*dist_x)/(2*_sigma_dist_x*_sigma_dist_x) - (dist_y*dist_y)/(2*_sigma_dist_y*_sigma_dist_y) - (dist_z*dist_z)/(2*_sigma_dist_z*_sigma_dist_z));
	}
      }
    }
  }
}

template<typename VoxelType> void irtkPatchBasedMSDLabelFusion<VoxelType>::Run()
{
  // Best matched patch for each atlas
  struct PATCH *best_patch = new struct PATCH[this->_number_of_atlases];

  int i = 0;
  for(int z=0; z<this->_Z; z++){
    for(int y=0; y<this->_Y; y++){
      for(int x=0; x<this->_X; x++, i++){
	if(this->_already_ptr[i] == 1){
	  continue;
	}
	
	// For each atlas
	for(int n=0; n<this->_number_of_atlases; n++){
	  // The best matched patch for this atlas
	  struct PATCH *p = &best_patch[n];
	  p->msd = 1E10;

	  // Search in a small volume for voxels which require bound-check
	  if(this->_bound_ptr[i] == 1){
	    int j = 0;
	    for(int dz=-this->_sz; dz<=this->_sz; dz++){
	      for(int dy=-this->_sy; dy<=this->_sy; dy++){
		for(int dx=-this->_sx; dx<=this->_sx; dx++, j++){
		  // Atlas patch centre
		  int x2 = x + dx;
		  int y2 = y + dy;
		  int z2 = z + dz;
		  int i2 = i + this->_offset_search[j];

		  if((x2 >= 0) && (x2 < this->_X)
		     && (y2 >= 0) && (y2 < this->_Y)
		     && (z2 >= 0) && (z2 < this->_Z)){
		    // Compute the mean squared difference between target and patch voxels
		    double msd = 0;
		    int n_voxel_count = 0;

		    int j2 = 0;
		    for(int dz2=-this->_pz; dz2<=this->_pz; dz2++){
		      for(int dy2=-this->_py; dy2<=this->_py; dy2++){
			for(int dx2=-this->_px; dx2<=this->_px; dx2++, j2++){
			  int x_target = x + dx2;
			  int y_target = y + dy2;
			  int z_target = z + dz2;

			  int x_atlas = x2 + dx2;
			  int y_atlas = y2 + dy2;
			  int z_atlas = z2 + dz2;

			  int i_target = i + this->_offset_patch[j2];
			  int i_atlas = i2 + this->_offset_patch[j2];

			  if((x_target >= 0) && (x_target < this->_X)
			     && (y_target >= 0) && (y_target < this->_Y)
			     && (z_target >= 0) && (z_target < this->_Z)
			     && (x_atlas >= 0) && (x_atlas < this->_X)
			     && (y_atlas >= 0) && (y_atlas < this->_Y)
			     && (z_atlas >= 0) && (z_atlas < this->_Z)){
			    double diff = this->_image_ptr[i_target] - this->_atlas_ptr[n][i_atlas];
			    msd += diff * diff;
			    n_voxel_count++;
			  } // if
			} // for dx2
		      } // for dy2
		    } // for dz2

		    msd /= n_voxel_count;

		    if(this->_best_patch_mode){
		      // Record the best matched patch and perform label fusion later
		      if(msd < p->msd){
			p->n = n;
			p->i = i2;
			p->x = x2;
			p->y = y2;
			p->z = z2;
			p->l = this->_label_ptr[n][i2];
			p->msd = msd;
		      }
		    }
		    else{
		      // Perform label fusion now
		      // Intensity weight
		      double w_inten = exp(- msd / _h);

		      // Total weight
		      double w = w_inten;

		      // If distance weight is also used
		      if(_use_dist_wt){
			double w_dist = _weight_dist[j];
			w = w_inten * w_dist;
		      }

		      // Fuse the atlas labels
		      FuseAtlasLabelWithBoundCheck(x, y, z, i, n, x2, y2, z2, i2, w);
		    } // if-else _best_patch_mode
		  } // if x2 >= 0
		} // for dz
	      } // for dy
	    } // for dx

	    // Fuse the label from the best matched patch for this atlas
	    if(this->_best_patch_mode){
	      int i2 = p->i;
	      int x2 = p->x;
	      int y2 = p->y;
	      int z2 = p->z;
	      double msd = p->msd;

	      // Intensity weight
	      double w_inten = exp(- msd / _h);

	      // Total weight
	      // Do not use distance weight in this mode
	      double w = w_inten;
	      
	      // Fuse the atlas labels
	      FuseAtlasLabelWithBoundCheck(x, y, z, i, n, x2, y2, z2, i2, w);
	    }
	  }
	  // Search in a small volume for voxels which do not require bound-check
	  else{
	    for(int j=0; j<this->_n_search; j++){
	      // Atlas patch centre
	      int i2 = i + this->_offset_search[j];
	      
	      // Compute the mean squared difference between target and patch voxels
	      double msd = 0;

	      for(int j2=0; j2<this->_n_patch; j2++){
		int i_target = i + this->_offset_patch[j2];
		int i_atlas = i2 + this->_offset_patch[j2];

		double diff = this->_image_ptr[i_target] - this->_atlas_ptr[n][i_atlas];
		msd += diff * diff;
	      }
	      msd /= this->_n_patch;

	      if(this->_best_patch_mode){
		// Record the best matched patch and perform label fusion later
		if(msd < p->msd){
		  p->n = n;
		  p->i = i2;
		  p->x = -1;
		  p->y = -1;
		  p->z = -1;
		  p->l = this->_label_ptr[n][i2];
		  p->msd = msd;
		}
	      }
	      else{
		// Perform label fusion now
		// Intensity weight
		double w_inten = exp(- msd / _h);

		// Total weight
		double w = w_inten;

		// If distance weight is also used
		if(_use_dist_wt){
		  double w_dist = _weight_dist[j];
		  w = w_inten * w_dist;
		}

		// Fuse the atlas labels
		FuseAtlasLabelWithoutBoundCheck(i, n, i2, w);
	      } // if-else _best_patch_mode
	    } // for j

	    // Fuse the label from the best matched patch for this atlas
	    if(this->_best_patch_mode){
	      int i2 = p->i;
	      double msd = p->msd;

	      // Intensity weight
	      double w_inten = exp(- msd / _h);

	      // Total weight
	      // Do not use distance weight in this mode
	      double w = w_inten;
	      
	      // Fuse the atlas labels
	      FuseAtlasLabelWithoutBoundCheck(i, n, i2, w);
	    }
	  } // if-else bound
	} // for each atlas
      } // for x
    } // for y
  } // for z

  // Label with the maximum vote
  i = 0;
  for(int z=0; z<this->_Z; z++){
    for(int y=0; y<this->_Y; y++){
      for(int x=0; x<this->_X; x++, i++){
	if(this->_already_ptr[i] == 1){
	  continue;
	}

	double max_vote = -1E10;
	double sum_vote = 0;
	int max_k = 0;
	for(int k=0; k<this->_number_of_classes; k++){
	  double vote = this->_prob_ptr[k][i];
	  if(vote > max_vote){
	    max_vote = vote;
	    max_k = k;
	  }
	  sum_vote += vote;
	}

	set<int>::iterator it = this->_label_set.begin();
	advance(it, max_k);
	this->_seg_ptr[i] = *it;

	// Probability
	for(int k=0; k<this->_number_of_classes; k++){
	  this->_prob_ptr[k][i] /= sum_vote;
	}
      } // for x
    } // for y
  } // for z

  // Free memory
  if(best_patch != NULL){
    delete []best_patch;
  }
}

// Fuse the label from the atlas patch (n, x2, y2, z2, i2) to the target patch (x, y, z, i) with weight w
template<typename VoxelType> inline void irtkPatchBasedMSDLabelFusion<VoxelType>::FuseAtlasLabelWithBoundCheck(int x, int y, int z, int i, int n, int x2, int y2, int z2, int i2, double w)
{
  // Multi-point mode
  if(this->_multipoint_mode){
    // Fuse labels from the whole atlas patch
    int j2 = 0;
    for(int dz2=-this->_pz; dz2<=this->_pz; dz2++){
      for(int dy2=-this->_py; dy2<=this->_py; dy2++){
	for(int dx2=-this->_px; dx2<=this->_px; dx2++, j2++){
	  int x_target = x + dx2;
	  int y_target = y + dy2;
	  int z_target = z + dz2;
	  
	  int x_atlas = x2 + dx2;
	  int y_atlas = y2 + dy2;
	  int z_atlas = z2 + dz2;
	  
	  int i_target = i + this->_offset_patch[j2];
	  int i_atlas = i2 + this->_offset_patch[j2];

	  if((x_target >= 0) && (x_target < this->_X)
	     && (y_target >= 0) && (y_target < this->_Y)
	     && (z_target >= 0) && (z_target < this->_Z)
	     && (x_atlas >= 0) && (x_atlas < this->_X)
	     && (y_atlas >= 0) && (y_atlas < this->_Y)
	     && (z_atlas >= 0) && (z_atlas < this->_Z)){
	    int l = this->_label_ptr[n][i_atlas];
	    this->_prob_ptr[this->_label_idx[l]][i_target] += w;
	  } // if
	} // for dx2
      } // for dy2
    } // for dz2
  }
  // Single-point mode
  else{
    // Only fuse the label from the atlas patch centre
    int l = this->_label_ptr[n][i2];
    this->_prob_ptr[this->_label_idx[l]][i] += w;
  }
}

// Fuse the label from the atlas patch (n, i2) to the target patch (i) with weight w
template<typename VoxelType> inline void irtkPatchBasedMSDLabelFusion<VoxelType>::FuseAtlasLabelWithoutBoundCheck(int i, int n, int i2, double w)
{
  // Multi-point mode
  if(this->_multipoint_mode){
    // Fuse labels from the whole atlas patch
    for(int j2=0; j2<this->_n_patch; j2++){
      int i_target = i + this->_offset_patch[j2];
      int i_atlas = i2 + this->_offset_patch[j2];
      
      int l = this->_label_ptr[n][i_atlas];
      this->_prob_ptr[this->_label_idx[l]][i_target] += w;
    }
  }
  // Single-point mode
  else{
    // Only fuse the label from the atlas patch centre
    int l = this->_label_ptr[n][i2];
    this->_prob_ptr[this->_label_idx[l]][i] += w;
  }
}

template<typename VoxelType> bool irtkPatchBasedMSDLabelFusion<VoxelType>::Read(char *buffer1, char *buffer2)
{
  bool ok = false;

  if(strcmp(buffer1, "h") == 0){
    this->_h = atof(buffer2);
    cout << "h is ... " << this->_h << endl;
    ok = true;
  }

  if(strcmp(buffer1, "Use distance weight") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_dist_wt = false;
      cout << "Use distance weight is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_dist_wt = true;
      cout << "Use distance weight is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(ok == false){
    ok = this->irtkPatchBasedLabelFusionBase<VoxelType>::Read(buffer1, buffer2);
  }

  return ok;
}

#endif
