/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONTEXTSVMLABELFUSION_TXX
#define _IRTKCONTEXTSVMLABELFUSION_TXX

template<typename VoxelType> irtkContextSVMLabelFusion<VoxelType>::irtkContextSVMLabelFusion()
{
  _n_support = 0;
  _offset_support = NULL;
  _weight_support = NULL;

  _svm_dx = 1;
  _svm_dy = 1;
  _svm_dz = 1;
  _svm_X = 0;
  _svm_Y = 0;
  _svm_Z = 0;
  _gamma = 0.01;
  _C = 1;
}

template<typename VoxelType> irtkContextSVMLabelFusion<VoxelType>::~irtkContextSVMLabelFusion()
{
  if(_offset_support != NULL){
    delete []_offset_support;
    _offset_support = NULL;
  }

  if(_weight_support != NULL){
    delete []_weight_support;
    _weight_support = NULL;
  }
}

template<typename VoxelType> const char *irtkContextSVMLabelFusion<VoxelType>::NameOfClass()
{
  return "irtkContextSVMLabelFusion";
}

template<typename VoxelType> void irtkContextSVMLabelFusion<VoxelType>::Initialize()
{
  this->irtkContextPatchBasedLabelFusion<VoxelType>::Initialize();

  // The number of SVM control points along each dimension
  _svm_X = floor((this->_X-1) / (double)_svm_dx) + 1;
  _svm_Y = floor((this->_Y-1) / (double)_svm_dy) + 1;
  _svm_Z = floor((this->_Z-1) / (double)_svm_dz) + 1;

  // Support region of each SVM control point
  _n_support = (2*_svm_dx - 1) * (2*_svm_dy - 1) * (2*_svm_dz - 1);

  _offset_support = new int[_n_support];
  _weight_support = new double[_n_support];
  
  int j = 0;
  for(int dz=-(_svm_dz-1); dz<=(_svm_dz-1); dz++){
    for(int dy=-(_svm_dy-1); dy<=(_svm_dy-1); dy++){
      for(int dx=-(_svm_dx-1); dx<=(_svm_dx-1); dx++, j++){
	_offset_support[j] = dz * this->_offset_z + dy * this->_offset_y + dx * this->_offset_x;
	_weight_support[j] = (1 - fabs(dx) / (double)_svm_dx) * (1 - fabs(dy) / (double)_svm_dy) * (1 - fabs(dz) / (double)_svm_dz);
      }
    }
  }
}

template<typename VoxelType> void irtkContextSVMLabelFusion<VoxelType>::Run()
{
  // Train each SVM, apply it, and delete the SVM.
  // Do not store each SVM due to memory saving purpose. If we store each SVM in a lattice, 
  // we will use huge memory (about 1G for cropped cardiac image and 60G for brain image).
  for(int sk=0; sk<_svm_Z; sk++){
    for(int sj=0; sj<_svm_Y; sj++){
      for(int si=0; si<_svm_X; si++){
	// Control point
	int x = si * _svm_dx;
	int y = sj * _svm_dy;
	int z = sk * _svm_dz;
	int i = z * this->_offset_z + y * this->_offset_y + x * this->_offset_x;

	// Train the classifier at this place
	SVM classifier;
	classifier.model = NULL;
	classifier.SV_memory = NULL;
	classifier.label = 0;
	classifier.feature_max = NULL;
	classifier.feature_min = NULL;

	TrainSVMClassifier(classifier, x, y, z, i);

	// Apply it to the neighouring voxels within its support (-/+ spacing)
	int xa = x - _svm_dx + 1;
	int xb = x + _svm_dx - 1;
	int ya = y - _svm_dy + 1;
	int yb = y + _svm_dy - 1;
	int za = z - _svm_dz + 1;
	int zb = z + _svm_dz - 1;

	// Iterate over the voxels in the support
	int j = 0;
	for(int z2=za; z2<=zb; z2++){
	  for(int y2=ya; y2<=yb; y2++){
	    for(int x2=xa; x2<=xb; x2++, j++){
	      int i2 = i + _offset_support[j];

	      if((x2 >= 0) && (x2 < this->_X)
		 && (y2 >= 0) && (y2 < this->_Y)
		 && (z2 >= 0) && (z2 < this->_Z)){
		if(this->_already_ptr[i2]){
		  continue;
		}

		// Output label from the SVM classifier
		int l = this->ApplySVMClassifier(classifier, si, sj, sk, x2, y2, z2);

		// The fraction weight for linear interpolation
		double wt = _weight_support[j];
		
		// Accumulate the probability
		int k = this->_label_idx[l];
		this->_prob_ptr[k][i2] += wt;
	      } // if
	    } // for x2
	  } // for y2
	} // for z2

	// Release memory
	if(classifier.model != NULL){
	  svm_free_and_destroy_model(&classifier.model);
	  classifier.model = NULL;
	}
	if(classifier.SV_memory != NULL){
	  delete []classifier.SV_memory;
	  classifier.SV_memory = NULL;
	}
	if(classifier.feature_max != NULL){
	  delete []classifier.feature_max;
	  classifier.feature_max = NULL;
	}
	if(classifier.feature_min != NULL){
	  delete []classifier.feature_min;
	  classifier.feature_min = NULL;
	}
      } // for si
    } // for sj
  } // for sk

  // Segmentation for each voxel
  double *vote = new double[this->_number_of_classes];

  int i = 0;
  for(int z=0; z<this->_Z; z++){
    for(int y=0; y<this->_Y; y++){
      for(int x=0; x<this->_X; x++, i++){
	if(this->_already_ptr[i] == 1){
	  continue;
	}
	  
	// Vote at this voxel
	for(int k=0; k<this->_number_of_classes; k++){
	  vote[k] = this->_prob_ptr[k][i];
	}

	// Label with the maximum vote
	double max_vote = -1E10;
	double sum_vote = 0;
	int max_k = 0;
	for(int k=0; k<this->_number_of_classes; k++){
	  if(vote[k] > max_vote){
	    max_vote = vote[k];
	    max_k = k;
	  }
	  sum_vote += vote[k];
	}
	
	set<int>::iterator it = this->_label_set.begin();
	advance(it, max_k);
	this->_seg_ptr[i] = *it;

	// Probability
	for(int k=0; k<this->_number_of_classes; k++){
	  this->_prob_ptr[k][i] = vote[k] / sum_vote;
	}
      } // for x
    } // for y
  } // for z

  delete []vote;
  vote = NULL;
}

// Train a SVM classifier at voxel (x, y, z)
template<typename VoxelType> inline void irtkContextSVMLabelFusion<VoxelType>::TrainSVMClassifier(SVM &node, int x, int y, int z, int i)
{
  // Atlas patch set
  vector<PATCH> atlas_patches;

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
  if(this->_K == 0){
    K = atlas_patches.size();
  }
  else{
    K = min((int)atlas_patches.size(), (int)this->_K);
  }

  // Check whether the K patches come from the same label
  bool unanimous = true;
  int l = atlas_patches[0].l;

  for(int i=1; i<K; i++){
    if(atlas_patches[i].l != l){
      unanimous = false;
      break;
    }
  }

  if(unanimous){
    // This node does not contain a SVM model. It is a unanimous label node.
    node.model = NULL;
    node.SV_memory = NULL;
    node.label = l;
    node.feature_max = NULL;
    node.feature_min = NULL;
  }
  else{
    // Otherwise, this node contains a SVM model.
    // We need to train the SVM here.

    // The labels and features of training samples
    int *training_labels = new int[K];
    double *training_features = new double[K * this->_n_feature];

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
      double *p_feature = training_features + k * this->_n_feature;
      this->ComputeFeature(this->_atlas[n], this->_atlas_prob[n], this->_atlas_grad[n], this->_atlas_median[n], x, y, z, x2, y2, z2, p_feature);
    }

    // Scale the feature range to be [-1, 1]
    node.feature_max = new double[this->_n_feature];
    node.feature_min = new double[this->_n_feature];

    for(int m=0; m<this->_n_feature; m++){
      node.feature_max[m] = -1E10;
      node.feature_min[m] = 1E10;
    }

    double *p_feature = training_features;
    for(int k=0; k<K; k++){
      for(int m=0; m<this->_n_feature; m++, p_feature++){
	double val = *p_feature;
	node.feature_max[m] = max(val, node.feature_max[m]);
	node.feature_min[m] = min(val, node.feature_min[m]);
      }
    }

    for(int k=0; k<K; k++){
      double *p_feature = training_features + k * this->_n_feature;
      this->ScaleFeature(p_feature, node.feature_max, node.feature_min);
    }

    // Convert the training features into the format that the SVM library needs
    struct svm_node *svm_training_features = new svm_node[K * (this->_n_feature + 1)];

    for(int k=0; k<K; k++){
      double *ptr1 = training_features + k * this->_n_feature;
      svm_node *ptr2 = svm_training_features + k * (this->_n_feature + 1);

      for(int m=0; m<this->_n_feature; m++, ptr1++, ptr2++){
	(*ptr2).index = m + 1; // Index for the SVM feature vector starts from 1
	(*ptr2).value = *ptr1; // Copy value
      }

      // End of feature vector
      (*ptr2).index = -1;
    }

    // Train the SVM
    struct svm_parameter param;
    struct svm_problem prob;

    // Default values
    param.svm_type = C_SVC;
    param.kernel_type = RBF;
    param.degree = 3;
    param.gamma = this->_gamma;
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 1000;
    param.C = this->_C;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;

    // Problem
    prob.l = K;
    prob.y = new double[K];
    prob.x = new svm_node *[K];

    for(int k=0; k<K; k++){
      prob.y[k] = training_labels[k];
      prob.x[k] = svm_training_features + k * (this->_n_feature + 1);
    }

    const char *error_msg = svm_check_parameter(&prob, &param);
    if(error_msg){
      cerr << "ERROR: " << error_msg << endl;
      exit(1);
    }

    // Training
    node.model = svm_train(&prob, &param);

    // model->SV stores the pointers to the support vectors in the training sample set. It
    // does not store the support vectors themselves.
    // As a result, we can not free the memory of svm_training_features (size O(K)) at the moment.
    // However, if we copy the support vectors to a separate memory block (size O(nSV)), we can
    // then free the memory of training_features. The new memory block will be freed at the end
    // of the whole programme.
    int nSV = node.model->l;
    node.SV_memory = new svm_node[nSV * (this->_n_feature + 1)];

    for(int k=0; k<nSV; k++){
      svm_node *dest = node.SV_memory + k * (this->_n_feature + 1); // This is the pointer to the new memory block
      svm_node *source = node.model->SV[k]; // This is a pointer to a sample in svm_training_features

      // Copy the support vector
      for(int m=0; m<(this->_n_feature + 1); m++){
	dest[m].index = source[m].index;
	dest[m].value = source[m].value;
      }

      // Modify the pointer
      node.model->SV[k] = dest;
    }

    // Free memory
    svm_destroy_param(&param);
    delete []prob.x;
    prob.x = NULL;
    delete []prob.y;
    prob.y = NULL;
    delete []training_labels;
    training_labels = NULL;
    delete []training_features;
    training_features = NULL;
    delete []svm_training_features;
    svm_training_features = NULL;
  }
}

// Apply the SVM at control point (si, sj, sk) to the patch at voxel (x2, y2, z2)
template<typename VoxelType> inline int irtkContextSVMLabelFusion<VoxelType>::ApplySVMClassifier(const SVM &node, int si, int sj, int sk, int x2, int y2, int z2)
{
  if(node.model == NULL){
    // This is a one-label node
    return node.label;
  }

  // The node position
  int x = si * _svm_dx;
  int y = sj * _svm_dy;
  int z = sk * _svm_dz;

  // The testing feature vector
  double *testing_feature = new double[this->_n_feature];
  this->ComputeFeature(this->_image, this->_init_prob_map, this->_image_grad, this->_image_median, x, y, z, x2, y2, z2, testing_feature);

  // Scale the testing feature
  this->ScaleFeature(testing_feature, node.feature_max, node.feature_min);

  // Convert the testing feature into the format that the SVM library needs
  struct svm_node *svm_testing_feature = new svm_node[this->_n_feature + 1];

  double *ptr1 = testing_feature;
  svm_node *ptr2 = svm_testing_feature;

  for(int m=0; m<this->_n_feature; m++, ptr1++, ptr2++){
    (*ptr2).index = m + 1; // Index for the SVM feature vector starts from 1
    (*ptr2).value = *ptr1; // Copy value
  }

  // End of feature vector
  (*ptr2).index = -1;

  // Predict the target label
  int predict_label = (int)svm_predict(node.model, svm_testing_feature);

  // Release memory
  delete []testing_feature;
  testing_feature = NULL;
  delete []svm_testing_feature;
  svm_testing_feature = NULL;

  return predict_label;
}

template<typename VoxelType> bool irtkContextSVMLabelFusion<VoxelType>::Read(char *buffer1, char *buffer2)
{
  bool ok = false;
  
  if(strcmp(buffer1, "SVM spacing in X (pixel)") == 0){
    _svm_dx = atoi(buffer2);
    cout << "SVM spacing in X is ... " << _svm_dx << endl;
    ok = true;
  }
  if(strcmp(buffer1, "SVM spacing in Y (pixel)") == 0){
    _svm_dy = atoi(buffer2);
    cout << "SVM spacing in Y is ... " << _svm_dy << endl;
    ok = true;
  }
  if(strcmp(buffer1, "SVM spacing in Z (pixel)") == 0){
    _svm_dz = atoi(buffer2);
    cout << "SVM spacing in Z is ... " << _svm_dz << endl;
    ok = true;
  }

  if(strcmp(buffer1, "gamma") == 0){
    _gamma = atof(buffer2);
    cout << "gamma is ... " << _gamma << endl;
    ok = true;
  }

  if(strcmp(buffer1, "C") == 0){
    _C = atof(buffer2);
    cout << "C is ... " << _C << endl;
    ok = true;
  }

  if(ok == false){
    ok = this->irtkContextPatchBasedLabelFusion<VoxelType>::Read(buffer1, buffer2);
  }

  return ok;
}

#endif
