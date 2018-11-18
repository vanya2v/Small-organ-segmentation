/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

/*
 * Multi-atlas label fusion.
 */

#include "irtkVoxel.h"
#include "irtkLabelFusionBase.h"
#include "irtkMajorityVoteLabelFusion.h"
#include "irtkPatchBasedMSDLabelFusion.h"
#include "irtkContextPatchBasedLabelFusion.h"
#include "irtkContextSVMLabelFusion.h"
#include <time.h>

void usage(char *command)
{
  cerr << "Usage: " << command << " image n_atlas atlas_1 ... atlas_n label_1 ... label_n output <options>" << endl;
  cerr << "Label fusion from multiple atlases." << endl;
  cerr << endl;

  cerr << "image: target image to be segmented" << endl;
  cerr << "n_atlas:                 number of atlases" << endl;
  cerr << "atlas_1 ... atlas_n:     atlas images (warped to the target image space)" << endl;
  cerr << "label_1 ... label_n:     atlas label maps" << endl;
  cerr << "output:                  segmentation given by label fusion" << endl;
  cerr << endl;

  cerr << "<options> can be one or more of the following:" << endl;
  cerr << "<-method name>           MV = majority vote" << endl;
  cerr << "                         PB = patch-based segmentation" << endl;
  cerr << "                         PBAF = patch-based segmentation with augmented features" << endl;
  cerr << "                         SVMAF = SVM segmentation with augmented features" << endl;
  cerr << "<-par file>              Read parameters (such as patch size, search volume etc.) from a file." << endl;
  cerr << "<-input_prob file>       Initial probability map estimate." << endl;
  cerr << "<-output_prob file>      Output the probability map for each class. The probability map can be used for EM or graphcuts segmentation refinement. The output is a single 4D Nifit image. Its temporal dimension represents each tissue class." << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  // Check command line
  char *command = argv[0];
  if(argc < 6){
    usage(command);
  }

  // Input and output
  char *image_name = argv[1];
  argc--;
  argv++;
  int n_atlas = atoi(argv[1]);
  argc--;
  argv++;
  char **atlas_name = new char *[n_atlas];
  for(int i=0; i<n_atlas; i++){
    atlas_name[i] = argv[1];
    argc--;
    argv++;
    if((argc == 1) && (i != (n_atlas-1))){
      cerr << "The number of atlas images is not enough!" << endl;
      cerr << endl;
      usage(command);
    }
  }
  char **label_name = new char *[n_atlas];
  for(int i=0; i<n_atlas; i++){
    label_name[i] = argv[1];
    argc--;
    argv++;
    if((argc == 1) && (i != (n_atlas-1))){
      cerr << "The number of atlas labels is not enough!" << endl;
      cerr << endl;
      usage(command);
    }
  }
  char *seg_name = argv[1];
  argc--;
  argv++;

  // Label fusion filter
  irtkLabelFusionBase<irtkGreyPixel> *fusion = NULL;

  // Parameters
  bool ok;
  char *method_name = "";
  char *par_name = NULL;
  char *input_prob_name = NULL;
  char *output_prob_name = NULL;

  // Parse the parameters
  while(argc > 1){
    ok = false;
    if((ok == false) && (strcmp(argv[1], "-method") == 0)){
      argc--;
      argv++;
      method_name = argv[1];
      argc--;
      argv++;
      if(strcmp(method_name, "MV") == 0){
	fusion = new irtkMajorityVoteLabelFusion<irtkGreyPixel>;
      }
      else if(strcmp(method_name, "PB") == 0){
	fusion = new irtkPatchBasedMSDLabelFusion<irtkGreyPixel>;
      }
      else if(strcmp(method_name, "PBAF") == 0){
	fusion = new irtkContextPatchBasedLabelFusion<irtkGreyPixel>;
      }
      else if(strcmp(method_name, "SVMAF") == 0){
	fusion = new irtkContextSVMLabelFusion<irtkGreyPixel>;
      }
      else{
	cerr << "No valid method name specified!" << endl;
	usage(command);
      }
      ok = true;
    }
    if((ok == false) && (strcmp(argv[1], "-par") == 0)){
      argc--;
      argv++;
      par_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if((ok == false) && (strcmp(argv[1], "-input_prob") == 0)){
      argc--;
      argv++;
      input_prob_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if((ok == false) && (strcmp(argv[1], "-output_prob") == 0)){
      argc--;
      argv++;
      output_prob_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if(ok == false){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage(command);
    }
  }

  // Clock for counting computation time
  clock_t start_clock, end_clock;
  start_clock = clock();

  // Read the images
  // irtk will do the type casting from the input data type to irtkGreyPixel (short).
  cout << "Reading image ..." << endl;
  irtkGreyImage *image = new irtkGreyImage(image_name);
  
  cout << "Reading atlas images ..." << endl;
  irtkGreyImage **atlas = new irtkGreyImage *[n_atlas];
  for(int i=0; i<n_atlas; i++){
    cout << " " << i << ": " << atlas_name[i] << endl;
    atlas[i] = new irtkGreyImage(atlas_name[i]);
  }

  cout << "Reading atlas label maps ..." << endl;
  irtkGreyImage **label_map = new irtkGreyImage *[n_atlas];
  for(int i=0; i<n_atlas; i++){
    cout << " " << i << ": " << label_name[i] << endl;
    label_map[i] = new irtkGreyImage(label_name[i]);
  }

  // Construct the label fusion classifier
  fusion->SetInput(image);
  fusion->SetAtlas(n_atlas, atlas, label_map);
  if(par_name){
    fusion->Read(par_name);
  }
  if(input_prob_name){
    irtkRealImage *input_prob = new irtkRealImage(input_prob_name);
    fusion->SetInitialProbabilityMap(input_prob);
  }

  // Run the algorithm
  fusion->Initialize();
  fusion->Run();

  // Output results
  irtkGreyImage *seg = fusion->GetOutput();
  seg->Write(seg_name);

  if(output_prob_name){
    irtkRealImage *output_prob = fusion->GetProbabilityMap();
    output_prob->Write(output_prob_name);
  }

  // Free memory
  delete []atlas_name;
  delete []label_name;
  for(int i=0; i<n_atlas; i++){
    delete atlas[i];
  }
  delete []atlas;
  for(int i=0; i<n_atlas; i++){
    delete label_map[i];
  }
  delete []label_map;
  delete fusion;
  delete image;

  end_clock = clock();
  cout << "It takes " << (end_clock - start_clock) / (double)CLOCKS_PER_SEC << " seconds." << endl;

  return 0;
}
