/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _IRTKLABELFUSIONBASE_TXX
#define _IRTKLABELFUSIONBASE_TXX

template<typename VoxelType> irtkLabelFusionBase<VoxelType>::irtkLabelFusionBase()
{
  _atlas = NULL;
  _label_map = NULL;
  _atlas_ptr = NULL;
  _label_ptr = NULL;
  _segmentation = NULL;
  _prob_map = NULL;
  _prob_ptr = NULL; 
  _init_prob_map = NULL;
  _already = NULL;
}

template<typename VoxelType> irtkLabelFusionBase<VoxelType>::~irtkLabelFusionBase()
{
  if(_atlas != NULL){
    delete []_atlas;
    _atlas = NULL;
  }

  if(_label_map != NULL){
    delete []_label_map;
    _label_map = NULL;
  }

  if(_atlas_ptr != NULL){
    delete []_atlas_ptr;
    _atlas_ptr = NULL;
  }

  if(_label_ptr != NULL){
    delete []_label_ptr;
    _label_ptr = NULL;
  }

  if(_segmentation != NULL){
    delete _segmentation;
    _segmentation = NULL;
  }

  if(_prob_map != NULL){
    delete _prob_map;
    _prob_map = NULL;
  }

  if(_prob_ptr != NULL){
    delete []_prob_ptr;
    _prob_ptr = NULL;
  }

  if(_init_prob_map != NULL){
    delete _init_prob_map;
    _init_prob_map = NULL;
  }

  if(_already != NULL){
    delete _already;
    _already = NULL;
  }
}

template<typename VoxelType> const char *irtkLabelFusionBase<VoxelType>::NameOfClass()
{
  return "irtkLabelFusionBase";
}

template<typename VoxelType> void irtkLabelFusionBase<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image)
{
  _image = image;
  _image_ptr = (VoxelType *)_image->GetScalarPointer();
  _X = _image->GetX();
  _Y = _image->GetY();
  _Z = _image->GetZ();
  _number_of_voxels = _X * _Y * _Z;
  _xsize = _image->GetXSize();
  _ysize = _image->GetYSize();
  _zsize = _image->GetZSize();
}

template<typename VoxelType> void irtkLabelFusionBase<VoxelType>::SetAtlas(int n_atlas, irtkGenericImage<VoxelType> **atlas, irtkGreyImage **label_map)
{
  // Check the dimension of the atlases
  _number_of_atlases = n_atlas;

  _atlas = new irtkGenericImage<VoxelType> *[_number_of_atlases];
  _label_map = new irtkGreyImage *[_number_of_atlases];
  _atlas_ptr = new VoxelType *[_number_of_atlases];
  _label_ptr = new irtkGreyPixel *[_number_of_atlases];

  for(int n=0; n<_number_of_atlases; n++){
    _atlas[n] = atlas[n]; // operator= for assignment of pointer
    _label_map[n] = label_map[n];
    _atlas_ptr[n] = (VoxelType *)_atlas[n]->GetScalarPointer();
    _label_ptr[n] = (irtkGreyPixel *)_label_map[n]->GetScalarPointer();

    if((_atlas[n]->GetX() != _X)
       || (_atlas[n]->GetY() != _Y)
       || (_atlas[n]->GetZ() != _Z)){
      cerr << "Dimension mismatch between target image and atlas image [" << n << "] !" << endl;
      exit(1);
    }

    if((_label_map[n]->GetX() != _X)
       || (_label_map[n]->GetY() != _Y)
       || (_label_map[n]->GetZ() != _Z)){
      cerr << "Dimension mismatch between target image and atlas label map [" << n << "] !" << endl;
      exit(1);
    }
  }

  // Count the label set
  for(int n=0; n<_number_of_atlases; n++){
    irtkGreyPixel *ptr = _label_ptr[n];
    
    for(int i=0; i<_number_of_voxels; i++, ptr++){
      _label_set.insert(*ptr);
    }
  }

  _number_of_classes = _label_set.size();

  cout << _number_of_classes << " unique labels found. They are namely:" << endl;
  for(set<int>::iterator it=_label_set.begin(); it!=_label_set.end(); ++it){
    cout << ' ' << *it;
  }
  cout << endl << endl;

  // The mapping from a label to its index in the labelset
  int index = 0;
  for(set<int>::iterator it=_label_set.begin(); it!=_label_set.end(); ++it, ++index){
    _label_idx.insert(pair<int, int>(*it, index));
  }

  cout << "Mapping from label to index:" << endl;
  for(map<int, int>::iterator it=_label_idx.begin(); it!=_label_idx.end(); ++it){
    cout << it->first << " => " << it->second << endl;
  }
  cout << endl;
}

template<typename VoxelType> void irtkLabelFusionBase<VoxelType>::Initialize()
{
  // Initialise segmentation
  if(_segmentation != NULL){
    delete _segmentation;
  }
  _segmentation = new irtkGreyImage;
  _segmentation->Initialize(_image->GetImageAttributes());
  _seg_ptr = (irtkGreyPixel *)_segmentation->GetScalarPointer();

  // Initialise the probabilistic label map
  // Use the t-dimension to store the probability vector for each class
  irtkImageAttributes attr = _image->GetImageAttributes();
  attr._t = _number_of_classes;

  if(_prob_map != NULL){
    delete _prob_map;
  }
  _prob_map = new irtkRealImage;
  _prob_map->Initialize(attr);

  _prob_ptr = new irtkRealPixel *[_number_of_classes];
  for(int k=0; k<_number_of_classes; k++){
    _prob_ptr[k] = (irtkRealPixel *)_prob_map->GetScalarPointer(0, 0, 0, k);
  }
  
  for(int k=0; k<_number_of_classes; k++){
    irtkRealPixel *ptr = _prob_ptr[k];
    for(int i=0; i<_number_of_voxels; i++, ptr++){
      *ptr = 0;
    }
  }

  // Create the mask of operation
  // If a target voxel is unanimously labelled by all the atlases, we do not perform sophisticated label fusion to save time.
  if(_already != NULL){
    delete _already;
  }
  _already = new irtkByteImage;
  _already->Initialize(_image->GetImageAttributes());
  _already_ptr = (irtkBytePixel *)_already->GetScalarPointer();
  int n_already_voxels = 0;

  for(int i=0; i<_number_of_voxels; i++){
    // Label at this voxel from atlas 0
    int l = _label_ptr[0][i];

    // Check whether the label is the same for the other atlases
    bool unanimous;

    if(_number_of_atlases == 1){
      // If only ONE atlas is available, we need to perform label fusion for all the voxels.
      unanimous = false;
    }
    else{
      // If there is more than one atlas, we check the unanimity.
      unanimous = true;

      for(int n=1; n<_number_of_atlases; n++){
	int l2 = _label_ptr[n][i];
	if(l2 != l){
	  unanimous = false;
	  break;
	}
      }
    }

    if(unanimous){
      _already_ptr[i] = 1;
      _seg_ptr[i] = l;
      _prob_ptr[_label_idx[l]][i] = 1.0;
      n_already_voxels++;
    }
    else{
      _already_ptr[i] = 0;
    }
  }

  double percent_already = (double)n_already_voxels * 100.0 / _number_of_voxels;
  cout << "The number of voxels is " << _number_of_voxels << "." << endl;
  cout << n_already_voxels << " voxels (" << percent_already << "%) are unanimously labelled." << endl;
  cout << (_number_of_voxels - n_already_voxels) << " voxels (" << (100 - percent_already) << "%) will be applied multi-atlas label fusion." << endl;
}

template<typename VoxelType> void irtkLabelFusionBase<VoxelType>::Run()
{
}

template<typename VoxelType> irtkGreyImage *irtkLabelFusionBase<VoxelType>::GetOutput()
{
  return _segmentation;
}

template<typename VoxelType> irtkRealImage *irtkLabelFusionBase<VoxelType>::GetProbabilityMap()
{
  return _prob_map;
}

template<typename VoxelType> void irtkLabelFusionBase<VoxelType>::SetInitialProbabilityMap(irtkRealImage *init_prob_map)
{
  _init_prob_map = init_prob_map;
}

template<typename VoxelType> void irtkLabelFusionBase<VoxelType>::Read(char *filename)
{
  char buffer1[255], *buffer2;
  ifstream from(filename);

  // Open file
  if(!from){
    cerr << "irtkLabelFusionBase::Read(): Can't open file " << filename << endl;
    exit(1);
  }
  cout << "--------------------------------------------------" << endl;
  cout << "Reading parameters ..." << endl;

  // Parse each uncommented line
  while(from.eof() != true){
    if(this->ReadLine(from, buffer1, buffer2) != 0){
      if(this->Read(buffer1, buffer2) == false){
        cerr << "Couldn't parse line: " << buffer1 << " = " << buffer2 << endl;
	exit(1);
      }
    }
  }
  cout << endl;
}


// Parse a parameter line like this,
// A       =  B
// ^          ^
// buffer1    buffer2
//
// This function was rewritten from read_line() in registration/src/irtkUtil.cc,
// in which buffer1 is wrongly parsed as the whole line, instead of A.
template<typename VoxelType> int irtkLabelFusionBase<VoxelType>::ReadLine(istream &in, char *buffer1, char *&buffer2)
{
  char c;

  do{
    if(in.eof() == true) return 0;
    in.getline(buffer1, 255);
    c = buffer1[0];
  } while((strlen(buffer1) == 0) || (c == '#') || (c == 13));

  // Separator
  if((buffer2 = strchr(buffer1, '=')) == NULL){
    cerr << "Not valid line format" << endl;
    exit(1);
  }
  
  // Parse the first part
  char *p = buffer2;
  do{
    p--;
    if(p < buffer1){
      cerr << "No content in buffer1" << endl;
      exit(1);
    }
  } while((*p == ' ') || (*p == '\t'));
  *(p+1) = '\0';

  // Parse the second part
  do{
    buffer2++;
  } while((*buffer2 == ' ') || (*buffer2 == '\t'));

  return strlen(buffer1);
}

template<typename VoxelType> bool irtkLabelFusionBase<VoxelType>::Read(char *buffer1, char *buffer2)
{
  bool ok = false;

  if(strcmp(buffer1, "Debug mode") == 0){
    if((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)){
      this->_debug_mode = false;
      cout << "Debug mode is ... false" << endl;
    }
    else if((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)){
      this->_debug_mode = true;
      cout << "Debug mode is ... true" << endl;
    }
    else{
      cerr << "Can't read boolean value = " << buffer2 << endl;
      exit(1);
    }
    ok = true;
  }

  return ok;
}

#endif
