/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _IRTKPATCHBASEDLABELFUSIONBASE_TXX
#define _IRTKPATCHBASEDLABELFUSIONBASE_TXX

struct PATCH {
  int n;  // Atlas index
  int i;  // Voxel index in 1D
  int x;  // Voxel index in 3D
  int y;
  int z;
  int l;  // Atlas label
  double msd; // Mean squared difference
};

bool patch_comp(PATCH a, PATCH b)
{
  return (a.msd < b.msd);
}

template<typename VoxelType> irtkPatchBasedLabelFusionBase<VoxelType>::irtkPatchBasedLabelFusionBase()
{
  _patch_size_x = 1;
  _patch_size_y = 1;
  _patch_size_z = 1;
  _px = 0;
  _py = 0;
  _pz = 0;
  _search_vol_x = 1;
  _search_vol_y = 1;
  _search_vol_z = 1;
  _sx = 0;
  _sy = 0;
  _sz = 0;
  _best_patch_mode = false;
  _multipoint_mode = false;

  _offset_patch = NULL;
  _offset_search = NULL;
  _bound_check = NULL;
}

template<typename VoxelType> irtkPatchBasedLabelFusionBase<VoxelType>::~irtkPatchBasedLabelFusionBase()
{
  if(_offset_patch != NULL){
    delete []_offset_patch;
    _offset_patch = NULL;
  }

  if(_offset_search != NULL){
    delete []_offset_search;
    _offset_search = NULL;
  }

  if(_bound_check != NULL){
    delete _bound_check;
    _bound_check = NULL;
  }
}

template<typename VoxelType> const char *irtkPatchBasedLabelFusionBase<VoxelType>::NameOfClass()
{
  return "irtkPatchBasedLabelFusionBase";
}

template<typename VoxelType> void irtkPatchBasedLabelFusionBase<VoxelType>::Initialize()
{
  this->irtkLabelFusionBase<VoxelType>::Initialize();

  // Half size of the patch in pixels
  _px = (int)floor(_patch_size_x * 0.5);
  _py = (int)floor(_patch_size_y * 0.5);
  _pz = (int)floor(_patch_size_z * 0.5);

  cout << "Patch size in X = " << (2 * _px + 1) << " pixel." << endl;
  cout << "Patch size in Y = " << (2 * _py + 1) << " pixel." << endl;
  cout << "Patch size in Z = " << (2 * _pz + 1) << " pixel." << endl;

  // Half size of the search volume in pixels
  _sx = (int)floor(_search_vol_x * 0.5);
  _sy = (int)floor(_search_vol_y * 0.5);
  _sz = (int)floor(_search_vol_z * 0.5);

  cout << "Search volume size in X = " << (2 * _sx + 1) << " pixel." << endl;
  cout << "Search volume size in Y = " << (2 * _sy + 1) << " pixel." << endl;
  cout << "Search volume size in Z = " << (2 * _sz + 1) << " pixel." << endl;

  // Offset table for evaluating the patch voxels
  _n_patch = (2*_px + 1) * (2*_py + 1) * (2*_pz + 1);
  _offset_patch = new int[_n_patch];
  
  _offset_x = 1;
  _offset_y = this->_X;
  _offset_z = this->_X * this->_Y;

  int j = 0;
  for(int dz=-_pz; dz<=_pz; dz++){
    for(int dy=-_py; dy<=_py; dy++){
      for(int dx=-_px; dx<=_px; dx++, j++){
	_offset_patch[j] = dz * _offset_z + dy * _offset_y + dx * _offset_x;
      }
    }
  }

  // Offset table for evaluating the search volume
  _n_search = (2*_sx + 1) * (2*_sy + 1) * (2*_sz + 1);
  _offset_search = new int[_n_search];
  
  j = 0;
  for(int dz=-_sz; dz<=_sz; dz++){
    for(int dy=-_sy; dy<=_sy; dy++){
      for(int dx=-_sx; dx<=_sx; dx++, j++){
	_offset_search[j] = dz * _offset_z + dy * _offset_y + dx * _offset_x;
      }
    }
  }

  // --------------------
  // |        1         |
  // |  --------------  |
  // |  |            |  |
  // |  |            |  |
  // |  |     0      |  |
  // |  |            |  |
  // |  |            |  |
  // |  --------------  |
  // |        1         |
  // --------------------
  //
  // Create the mask of bound-check region, as shown above
  // 
  // If the target voxel is in this region (labelled 1), we need to perform boundary check for the atlas voxel index when applying an offset (search volume + patch voxel).
  // Outside this region (labelled 0), we can simply create the atlas index by applying the offset without any check.
  // 
  if(_bound_check != NULL){
    delete _bound_check;
  }
  _bound_check = new irtkByteImage;
  _bound_check->Initialize(this->_image->GetImageAttributes());
  _bound_ptr = (irtkBytePixel *)_bound_check->GetScalarPointer();

  int bound_x = _sx + _px;
  int bound_y = _sy + _py;
  int bound_z = _sz + _pz;
  int n_bound_check_voxels = 0;

  int i = 0;
  for(int z=0; z<this->_Z; z++){
    for(int y=0; y<this->_Y; y++){
      for(int x=0; x<this->_X; x++, i++){
	if((x >= bound_x) && (x < (this->_X - bound_x))
	   && (y >= bound_y) && (y < (this->_Y - bound_y))
	   && (z >= bound_z) && (z < (this->_Z - bound_z))){
	  _bound_ptr[i] = 0;
	}
	else{
	  _bound_ptr[i] = 1;

	  if(this->_already_ptr[i] == 0){
	    n_bound_check_voxels++;
	  }
	}
      }
    }
  }  

  double percent_bound = (double)n_bound_check_voxels * 100.0 / this->_number_of_voxels;
  cout << n_bound_check_voxels << " voxels (" << percent_bound << "%) will be applied image bound check." << endl;
}

template<typename VoxelType> void irtkPatchBasedLabelFusionBase<VoxelType>::Run()
{
}

// Check the voxel index
// If it is outside the image boundary, duplicate the boundary voxel.
template<typename VoxelType> inline void irtkPatchBasedLabelFusionBase<VoxelType>::DuplicateBoundaryVoxel(int &x, int &y, int &z)
{
  x = (x >= 0)? x : 0;
  x = (x < this->_X)? x : (this->_X-1);

  y = (y >= 0)? y : 0;
  y = (y < this->_Y)? y : (this->_Y-1);
   
  z = (z >= 0)? z : 0;
  z = (z < this->_Z)? z : (this->_Z-1);
}

template<typename VoxelType> bool irtkPatchBasedLabelFusionBase<VoxelType>::Read(char *buffer1, char *buffer2)
{
  bool ok = false;

  if(strcmp(buffer1, "Patch size in X (pixel)") == 0){
    this->_patch_size_x = atoi(buffer2);
    cout << "Patch size in X is ... " << this->_patch_size_x << endl;
    ok = true;
  }
  if(strcmp(buffer1, "Patch size in Y (pixel)") == 0){
    this->_patch_size_y = atoi(buffer2);
    cout << "Patch size in Y is ... " << this->_patch_size_y << endl;
    ok = true;
  }
  if(strcmp(buffer1, "Patch size in Z (pixel)") == 0){
    this->_patch_size_z = atoi(buffer2);
    cout << "Patch size in Z is ... " << this->_patch_size_z << endl;
    ok = true;
  }

  if(strcmp(buffer1, "Search volume size in X (pixel)") == 0){
    this->_search_vol_x = atoi(buffer2);
    cout << "Search volume size in X is ... " << this->_search_vol_x << endl;
    ok = true;
  }
  if(strcmp(buffer1, "Search volume size in Y (pixel)") == 0){
    this->_search_vol_y = atoi(buffer2);
    cout << "Search volume size in Y is ... " << this->_search_vol_y << endl;
    ok = true;
  }
  if(strcmp(buffer1, "Search volume size in Z (pixel)") == 0){
    this->_search_vol_z = atoi(buffer2);
    cout << "Search volume size in Z is ... " << this->_search_vol_z << endl;
    ok = true;
  }

  if(strcmp(buffer1, "Best patch mode") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_best_patch_mode = false;
      cout << "Best patch mode is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_best_patch_mode = true;
      cout << "Best patch mode is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(strcmp(buffer1, "Multipoint mode") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_multipoint_mode = false;
      cout << "Multipoint mode is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_multipoint_mode = true;
      cout << "Multipoint mode is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  if(ok == false){
    ok = this->irtkLabelFusionBase<VoxelType>::Read(buffer1, buffer2);
  }

  return ok;
}

#endif
