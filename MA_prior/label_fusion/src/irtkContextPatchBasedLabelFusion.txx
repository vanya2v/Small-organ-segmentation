/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONTEXTPATCHBASEDLABELFUSION_TXX
#define _IRTKCONTEXTPATCHBASEDLABELFUSION_TXX

template<typename VoxelType> irtkContextPatchBasedLabelFusion<VoxelType>::irtkContextPatchBasedLabelFusion()
{
  _image_grad = NULL;
  _atlas_grad = NULL;
  _image_median = NULL;
  _atlas_median = NULL;
  _atlas_prob = NULL;

  _use_spatial_feature = true;
  _use_intensity_feature = true;
  _use_gradient_feature = true;
  _3D_context = false;
  _use_context_feature = true;
  _use_binary_context_feature = false;
  _use_prob_context_feature = false;

  _n_ray = 8;
  _n_sample_per_ray = 4;
  _n_sample = 0;
  _start_radius = 4;
  _offset_context_x = NULL;
  _offset_context_y = NULL;
  _offset_context_z = NULL;

  _require_scaling = NULL;

  _h = 1;
  _K = 0;
}

template<typename VoxelType> irtkContextPatchBasedLabelFusion<VoxelType>::~irtkContextPatchBasedLabelFusion()
{
  if(_image_grad != NULL){
    delete _image_grad;
    _image_grad = NULL;
  }

  if(_atlas_grad != NULL){
    for(int n=0; n<this->_number_of_atlases; n++){
      if(_atlas_grad[n] != NULL){
	delete _atlas_grad[n];
      }
    }
    delete []_atlas_grad;
    _atlas_grad = NULL;
  }

  if(_image_median != NULL){
    delete _image_median;
    _image_median = NULL;
  }

  if(_atlas_median != NULL){
    for(int n=0; n<this->_number_of_atlases; n++){
      if(_atlas_median[n] != NULL){
	delete _atlas_median[n];
      }
    }
    delete []_atlas_median;
    _atlas_median = NULL;
  }

  if(_atlas_prob != NULL){
    for(int n=0; n<this->_number_of_atlases; n++){
      if(_atlas_prob[n] != NULL){
	delete _atlas_prob[n];
      }
    }
    delete []_atlas_prob;
    _atlas_prob = NULL;
  }

  if(_offset_context_x != NULL){
    delete []_offset_context_x;
    _offset_context_x = NULL;
  }

  if(_offset_context_y != NULL){
    delete []_offset_context_y;
    _offset_context_y = NULL;
  }

  if(_offset_context_z != NULL){
    delete []_offset_context_z;
    _offset_context_z = NULL;
  }

  if(_require_scaling != NULL){
    delete []_require_scaling;
    _require_scaling = NULL;
  }
}

template<typename VoxelType> const char *irtkContextPatchBasedLabelFusion<VoxelType>::NameOfClass()
{
  return "irtkContextPatchBasedLabelFusion";
}

template<typename VoxelType> void irtkContextPatchBasedLabelFusion<VoxelType>::Initialize()
{
  this->irtkPatchBasedLabelFusionBase<VoxelType>::Initialize();

  // Features for classification include:
  // (1) Spatial feature: shift from the target voxel
  _n_spatial_feature = _use_spatial_feature? 3 : 0;
	      
  // (2) Intensity feature
  _n_intensity_feature = _use_intensity_feature? this->_n_patch : 0;

  // (3) Gradient feature
  _n_gradient_feature = _use_gradient_feature? (this->_n_patch * 3) : 0;

  // Number of sample regions in the surroundings
  _n_sample = _3D_context? (_n_sample_per_ray * _n_ray * 3) : (_n_sample_per_ray * _n_ray);

  // (4) Contextual feature
  _n_context_feature = _use_context_feature? _n_sample : 0;

  // (5) Binary contextual feature
  _n_binary_context_feature = _use_binary_context_feature? _n_sample : 0;

  // (6) Probability contextual feature
  _n_prob_context_feature = _use_prob_context_feature? (_n_sample * this->_number_of_classes) : 0;

  // Total number of features
  _n_feature = _n_spatial_feature + _n_intensity_feature + _n_gradient_feature + _n_context_feature + _n_binary_context_feature + _n_prob_context_feature;

  // Determine which bits of the feature vector require scaling
  _require_scaling = new bool[_n_feature];
  int index = 0;
  
  for(int m=0; m<_n_spatial_feature; m++, index++){
    _require_scaling[index] = true;
  }

  for(int m=0; m<_n_intensity_feature; m++, index++){
    _require_scaling[index] = true;
  }

  for(int m=0; m<_n_gradient_feature; m++, index++){
    _require_scaling[index] = true;
  }

  for(int m=0; m<_n_context_feature; m++, index++){
    _require_scaling[index] = true;
  }

  for(int m=0; m<_n_binary_context_feature; m++, index++){
    _require_scaling[index] = false;
  }

  for(int m=0; m<_n_prob_context_feature; m++, index++){
    _require_scaling[index] = false;
  }

  // Pre-compute the gradient images
  _atlas_grad = new irtkRealImage *[this->_number_of_atlases];
  for(int n=0; n<this->_number_of_atlases; n++){
    _atlas_grad[n] = NULL;
  }

  if(_use_gradient_feature){
    // Create the gradient image filter (data type: float) which calculates the normalised gradient vector.
    // TODO: may need to be changed after updating irtk.
    irtkGradientImageFilter<irtkRealPixel> gradient(irtkGradientImageFilter<irtkRealPixel>::NORMALISED_GRADIENT_VECTOR);

    // Convert image data type to float
    irtkRealImage image = *this->_image; 

    // Compute the gradient image
    _image_grad = new irtkRealImage;
    gradient.SetInput(&image);
    gradient.SetOutput(_image_grad);
    gradient.Run();

    // Compute the atlas gradient images
    for(int n=0; n<this->_number_of_atlases; n++){
      // Convert image data type to float
      image = *this->_atlas[n];

      _atlas_grad[n] = new irtkRealImage;
      gradient.SetInput(&image);
      gradient.SetOutput(_atlas_grad[n]);
      gradient.Run();
    }
  }

  // Pre-compute the median images
  _atlas_median = new irtkGenericImage<VoxelType> *[this->_number_of_atlases];
  for(int n=0; n<this->_number_of_atlases; n++){
    _atlas_median[n] = NULL;
  }

  if(_use_context_feature || _use_binary_context_feature){
    // Create the median image filter
    irtkImageMedianFilter<VoxelType> median;

    if(_3D_context){
      median.SetWindowSize(1, 1, 1);
    }
    else{
      median.SetWindowSize(1, 1, 0);
    }

    // Compute the median image
    _image_median = new irtkGenericImage<VoxelType>;
    median.SetInput(this->_image);
    median.SetOutput(_image_median);
    median.Run();

    // Compute the atlas median images
    for(int n=0; n<this->_number_of_atlases; n++){
      _atlas_median[n] = new irtkGenericImage<VoxelType>;
      median.SetInput(this->_atlas[n]);
      median.SetOutput(_atlas_median[n]);
      median.Run();
    }
  }

  // Pre-compute the atlas probability maps
  _atlas_prob = new irtkRealImage *[this->_number_of_atlases];
  for(int n=0; n<this->_number_of_atlases; n++){
    _atlas_prob[n] = NULL;
  }

  if(_use_prob_context_feature){
    if(this->_init_prob_map == NULL){
      cerr << "Error: no initial probability map!" << endl;
      exit(1);
    }
    
    irtkImageAttributes attr = this->_prob_map->GetImageAttributes();

    for(int n=0; n<this->_number_of_atlases; n++){
      _atlas_prob[n] = new irtkRealImage;
      _atlas_prob[n]->Initialize(attr);
      
      for(int z=0; z<this->_Z; z++){
	for(int y=0; y<this->_Y; y++){
	  for(int x=0; x<this->_X; x++){
	    int l = this->_label_map[n]->Get(x, y, z);
	    int k = this->_label_idx[l];
	    _atlas_prob[n]->Put(x, y, z, k, 1);
	  }
	}
      }
    }
  }

  // Pre-compute the context sample positions
  if(_use_context_feature || _use_binary_context_feature || _use_prob_context_feature){
    _offset_context_x = new int[_n_sample];
    _offset_context_y = new int[_n_sample];
    _offset_context_z = new int[_n_sample];

    int index = 0;
    int n_plane = (_3D_context)? 3 : 1;

    // For each plane
    for(int i=0; i<n_plane; i++){
      // For each direction
      for(int j=0; j<_n_ray; j++){
	double degree = j * (2 * PI) / (double)_n_ray;
	double du = cos(degree);
	double dv = sin(degree);

	// For each distance
	int radius = _start_radius;
	for(int k=0; k<_n_sample_per_ray; k++){
	  // Offset
	  switch(i){
	  case 0:
	    // XY-plane
	    _offset_context_x[index] = (int)round(radius * du);
	    _offset_context_y[index] = (int)round(radius * dv);
	    _offset_context_z[index] = 0;
	    break;

	  case 1:
	    // YZ-plane
	    _offset_context_x[index] = 0;
	    _offset_context_y[index] = (int)round(radius * du);
	    _offset_context_z[index] = (int)round(radius * dv);
	    break;

	  case 2:
	    // ZX-plane
	    _offset_context_x[index] = (int)round(radius * dv);
	    _offset_context_y[index] = 0;
	    _offset_context_z[index] = (int)round(radius * du);
	  break;

	  default:
	    cerr << "Error: unknown plane index:" << i << "!" << endl;
	    exit(1);
	  }

	  // Increase the index
	  index++;

	  // Increase the distance from the patch centre
	  radius *= 2;
	} // for k
      } // for j
    } // for i
  }
}

template<typename VoxelType> void irtkContextPatchBasedLabelFusion<VoxelType>::Run()
{
  // Record the atlas patches and its similarity metric
  vector<PATCH> atlas_patches;

  int i = 0;
  for(int z=0; z<this->_Z; z++){
    for(int y=0; y<this->_Y; y++){
      for(int x=0; x<this->_X; x++, i++){
	if(this->_already_ptr[i] == 1){
	  continue;
	}

	// Initialise the atlas patch set
	atlas_patches.clear();

	// For each atlas
	for(int n=0; n<this->_number_of_atlases; n++){
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

		    // Add this patch
		    struct PATCH p;
		    p.n = n;
		    p.i = i2;
		    p.x = x2;
		    p.y = y2;
		    p.z = z2;
		    p.l = this->_label_ptr[n][i2];
		    p.msd = msd;

		    atlas_patches.push_back(p);
		  } // if
		} // for dz
	      } // for dy
	    } // for dx
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

	      // Add this patch
	      struct PATCH p;
	      p.n = n;
	      p.i = i2;
	      p.x = -1; // (x,y,z) is unknown here. We will compute it later.
	      p.y = -1;
	      p.z = -1;
	      p.l = this->_label_ptr[n][i2];
	      p.msd = msd;

	      atlas_patches.push_back(p);
	    } // for j
	  } // if-else bound
	} // for each atlas n

	// Sort the patches in the ascending order of MSD
	sort(atlas_patches.begin(), atlas_patches.end(), patch_comp);

	// Use the K most-similar patches
	int K;
	if(_K == 0){
	  K = atlas_patches.size();
	}
	else{
	  K = min((int)atlas_patches.size(), (int)_K);
	}

	// The labels and features of training samples
	int *training_labels = new int[K];
	double *training_features = new double[K * _n_feature];

	// For each sample in the training set
	for(int k=0; k<K; k++){
	  struct PATCH *p = &atlas_patches[k];
	  int n = p->n;
	  int i2 = p->i;
	  int x2 = p->x;
	  int y2 = p->y;
	  int z2 = p->z;
	  int l = p->l;

	  if(p->x == -1){
	    // Compute (x,y,z) from i
	    int tmp = i2;

	    // Hopefully, the compiler will optimise the code so that only one DIV instruction will be used.
	    z2 = tmp / this->_offset_z;
	    tmp = tmp % this->_offset_z;

	    y2 = tmp / this->_offset_y; 
	    x2 = tmp % this->_offset_y;
	  }

	  // Sample label
	  training_labels[k] = l;

	  // Sample feature
	  double *p_feature = training_features + k * _n_feature;
	  ComputeFeature(this->_atlas[n], _atlas_prob[n], _atlas_grad[n], _atlas_median[n], x, y, z, x2, y2, z2, p_feature);
	}
	
	// Scale the feature range to be [-1, 1]
	double *feature_max = new double[_n_feature];
	double *feature_min = new double[_n_feature];

	for(int m=0; m<_n_feature; m++){
	  feature_max[m] = -1E10;
	  feature_min[m] = 1E10;
	}

	double *p_feature = training_features;
	for(int k=0; k<K; k++){
	  for(int m=0; m<_n_feature; m++, p_feature++){
	    double val = *p_feature;
	    feature_max[m] = max(val, feature_max[m]);
	    feature_min[m] = min(val, feature_min[m]);
	  }
	}

	for(int k=0; k<K; k++){
	  double *p_feature = training_features + k * _n_feature;
	  ScaleFeature(p_feature, feature_max, feature_min);
	}

	// The testing feature vector
	double *testing_feature = new double[_n_feature];
	ComputeFeature(this->_image, this->_init_prob_map, _image_grad, _image_median, x, y, z, x, y, z, testing_feature);

	// Scale the testing feature
	ScaleFeature(testing_feature, feature_max, feature_min);

	// Weighted label fusion
	for(int k=0; k<K; k++){
	  struct PATCH *p = &atlas_patches[k];
	  int n = p->n;
	  int i2 = p->i;
	  int x2 = p->x;
	  int y2 = p->y;
	  int z2 = p->z;

	  double ssd = 0;
	  double *ptr1 = training_features + k * _n_feature;
	  double *ptr2 = testing_feature;

	  for(int m=0; m<_n_feature; m++, ptr1++, ptr2++){
	    double diff = *ptr1 - *ptr2;
	    ssd += diff * diff;
	  }
	  double w = exp(- ssd / _h);

	  // Multi-point mode
	  if(this->_multipoint_mode){
	    // Fuse labels from the whole atlas patch
	    if(this->_bound_ptr[i] == 1){
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
	    else{
	      for(int j2=0; j2<this->_n_patch; j2++){
		int i_target = i + this->_offset_patch[j2];
		int i_atlas = i2 + this->_offset_patch[j2];
      
		int l = this->_label_ptr[n][i_atlas];
		this->_prob_ptr[this->_label_idx[l]][i_target] += w;
	      }
	    }
	  }
	  // Single-point mode
	  else{
	    // Only fuse the label from the atlas patch centre
	    int l = this->_label_ptr[n][i2];
	    this->_prob_ptr[this->_label_idx[l]][i] += w;
	  } // if-else bound
	} // for k

	// Free memory
	delete []training_labels;
	training_labels = NULL;
	delete []training_features;
	training_features = NULL;
	delete []feature_max;
	feature_max = NULL;
	delete []feature_min;
	feature_min = NULL;
	delete []testing_feature;
	testing_feature = NULL;
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
}

// Compute the features
// (x, y, z): the origin of the current coordinate system, i.e. target voxel index
// (x2, y2, z2): the centre of a moving patch
template<typename VoxelType> void irtkContextPatchBasedLabelFusion<VoxelType>::ComputeFeature(irtkGenericImage<VoxelType> *image, irtkRealImage *prob, irtkRealImage *grad, irtkGenericImage<VoxelType> *median, int x, int y, int z, int x2, int y2, int z2, double *feature)
{
  // Pointer
  double *ptr = feature;

  // Spatial feature
  if(_use_spatial_feature){
    *ptr = (x2 - x) * this->_xsize;
    ptr++;
    *ptr = (y2 - y) * this->_ysize;
    ptr++;
    *ptr = (z2 - z) * this->_zsize;
    ptr++;
  }

  // Intensity feature
  if(_use_intensity_feature){
    ComputeIntensityFeature(image, x2, y2, z2, ptr);
    ptr += _n_intensity_feature;
  }

  // Gradient feature
  if(_use_gradient_feature){
    ComputeGradientFeature(grad, x2, y2, z2, ptr);
    ptr += _n_gradient_feature;
  }

  // Contextual feature
  if(this->_use_context_feature){
    ComputeContextualFeature(image, median, x2, y2, z2, ptr);
    ptr += _n_context_feature;
  }

  // Binary contextual feature
  if(this->_use_binary_context_feature){
    ComputeBinaryContextualFeature(image, median, x2, y2, z2, ptr);
    ptr += _n_binary_context_feature;
  }

  // Probability contextual feature
  if(this->_use_prob_context_feature){
    ComputeProbabilityContextualFeature(prob, x2, y2, z2, ptr);
    ptr += _n_prob_context_feature;
  }
}

// Compute intensity feature
template<typename VoxelType> inline void irtkContextPatchBasedLabelFusion<VoxelType>::ComputeIntensityFeature(irtkGenericImage<VoxelType> *image, int x, int y, int z, double *feature)
{
  double *ptr = feature;
  for(int dz2=-this->_pz; dz2<=this->_pz; dz2++){
    for(int dy2=-this->_py; dy2<=this->_py; dy2++){
      for(int dx2=-this->_px; dx2<=this->_px; dx2++){
	// Image voxel index
	int x3 = x + dx2;
	int y3 = y + dy2;
	int z3 = z + dz2;

	// Duplicate at the boundary
	this->DuplicateBoundaryVoxel(x3, y3, z3);

	// Voxel intensity
	double val = image->Get(x3, y3, z3);

	*ptr = val;
	ptr++;
      }
    }
  }
}

// Compute gradient feature
template<typename VoxelType> inline void irtkContextPatchBasedLabelFusion<VoxelType>::ComputeGradientFeature(irtkRealImage *grad, int x, int y, int z, double *feature)
{
  double *ptr = feature;
  for(int dz2=-this->_pz; dz2<=this->_pz; dz2++){
    for(int dy2=-this->_py; dy2<=this->_py; dy2++){
      for(int dx2=-this->_px; dx2<=this->_px; dx2++){
	// Image voxel index
	int x3 = x + dx2;
	int y3 = y + dy2;
	int z3 = z + dz2;

	// Duplicate at the boundary
	this->DuplicateBoundaryVoxel(x3, y3, z3);
		    
	// Normalised gradient
	*ptr = grad->Get(x3, y3, z3, 0);
	ptr++;
	*ptr = grad->Get(x3, y3, z3, 1);
	ptr++;
	*ptr = grad->Get(x3, y3, z3, 2);
	ptr++;
      }
    }
  }
}

// Compute contextual feature
template<typename VoxelType> inline void irtkContextPatchBasedLabelFusion<VoxelType>::ComputeContextualFeature(irtkGenericImage<VoxelType> *image, irtkGenericImage<VoxelType> *median, int x, int y, int z, double *feature)
{
  double *ptr = feature;
  double val_central = image->Get(x, y, z);

  for(int index=0; index<_n_sample; index++){
    int x3 = x + _offset_context_x[index];
    int y3 = y + _offset_context_y[index];
    int z3 = z + _offset_context_z[index];

    this->DuplicateBoundaryVoxel(x3, y3, z3);

    // The median intensity in the sample patch
    double val_median = median->Get(x3, y3, z3);

    // Record the difference between the sample intensity and the central intensity as a feature
    *ptr = val_median - val_central;
    ptr++;
  }
}

// Compute binary contextual feature
template<typename VoxelType> inline void irtkContextPatchBasedLabelFusion<VoxelType>::ComputeBinaryContextualFeature(irtkGenericImage<VoxelType> *image, irtkGenericImage<VoxelType> *median, int x, int y, int z, double *feature)
{
  double *ptr = feature;
  double val_central = image->Get(x, y, z);

  for(int index=0; index<_n_sample; index++){
    int x3 = x + _offset_context_x[index];
    int y3 = y + _offset_context_y[index];
    int z3 = z + _offset_context_z[index];

    this->DuplicateBoundaryVoxel(x3, y3, z3);

    // The median intensity in the sample patch
    double val_median = median->Get(x3, y3, z3);

    // Record the difference between the sample intensity and the central intensity as a feature
    *ptr = (val_median > val_central);
    ptr++;
  }
}

// Compute probability contextual feature
template<typename VoxelType> inline void irtkContextPatchBasedLabelFusion<VoxelType>::ComputeProbabilityContextualFeature(irtkRealImage *prob, int x, int y, int z, double *feature)
{
  double *ptr = feature;

  for(int index=0; index<_n_sample; index++){
    int x3 = x + _offset_context_x[index];
    int y3 = y + _offset_context_y[index];
    int z3 = z + _offset_context_z[index];

    this->DuplicateBoundaryVoxel(x3, y3, z3);

    // The probability in the sample voxel for each class
    // Questions: do we need to compute the average probability in the sample region?
    for(int k=0; k<this->_number_of_classes; k++){
      double val = prob->Get(x3, y3, z3, k);
      *ptr = val;
      ptr++;
    }
  }
}

// Scale each element of the feature vector to the range [-1, 1]
template<typename VoxelType> inline void irtkContextPatchBasedLabelFusion<VoxelType>::ScaleFeature(double *feature, double *feature_max, double *feature_min)
{
  double upper = 1;
  double lower = -1;

  double *ptr = feature;
  double *ptr_max = feature_max;
  double *ptr_min = feature_min;

  for(int m=0; m<_n_feature; m++, ptr++, ptr_max++, ptr_min++){
    if(_require_scaling[m]){
      double val = *ptr;
      double val_max = *ptr_max;
      double val_min = *ptr_min;

      if(val_max == val_min){
	continue;
      }
      
      val = lower + (upper-lower) * (val-val_min) / (val_max-val_min);
      *ptr = val;
    }
  }
}

template<typename VoxelType> bool irtkContextPatchBasedLabelFusion<VoxelType>::Read(char *buffer1, char *buffer2)
{
  bool ok = false;

  if(strcmp(buffer1, "Use spatial feature") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_spatial_feature = false;
      cout << "Use spatial feature is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_spatial_feature = true;
      cout << "Use spatial feature is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Use intensity feature") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_intensity_feature = false;
      cout << "Use intensity feature is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_intensity_feature = true;
      cout << "Use intensity feature is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Use gradient feature") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_gradient_feature = false;
      cout << "Use gradient feature is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_gradient_feature = true;
      cout << "Use gradient feature is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "3D context mode") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_3D_context = false;
      cout << "3D context mode is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_3D_context = true;
      cout << "3D context mode is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Use context feature") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_context_feature = false;
      cout << "Use context feature is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_context_feature = true;
      cout << "Use context feature is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Use binary context feature") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_binary_context_feature = false;
      cout << "Use binary context feature is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_binary_context_feature = true;
      cout << "Use binary context feature is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Use probability context feature") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_use_prob_context_feature = false;
      cout << "Use probability context feature is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_use_prob_context_feature = true;
      cout << "Use probability context feature is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Number of rays") == 0){
    this->_n_ray = atoi(buffer2);
    cout << "Number of rays is ... " << this->_n_ray << endl;
    ok = true;
  }

  if(strcmp(buffer1, "Number of samples per ray") == 0){
    this->_n_sample_per_ray = atoi(buffer2);
    cout << "Number of samples per ray is ... " << this->_n_sample_per_ray << endl;
    ok = true;
  }

  if(strcmp(buffer1, "Starting radius") == 0){
    this->_start_radius = atoi(buffer2);
    cout << "Starting radius is ... " << this->_start_radius << endl;
    ok = true;
  }

  if(strcmp(buffer1, "h") == 0){
    this->_h = atof(buffer2);
    cout << "h is ... " << this->_h << endl;
    ok = true;
  }

  if(strcmp(buffer1, "K") == 0){
    this->_K = atoi(buffer2);
    cout << "K is ... " << this->_K << endl;
    ok = true;
  }

  if(ok == false){
    ok = this->irtkPatchBasedLabelFusionBase<VoxelType>::Read(buffer1, buffer2);
  }

  return ok;
}

#endif

