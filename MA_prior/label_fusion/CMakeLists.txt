cmake_minimum_required(VERSION 2.6)
PROJECT(label_fusion)

# Please download and compile IRTK first
# This programme will be in the next release of IRTK.
# If not, you can still compile the programme using the following lines.
SET(IRTK_DIR /vol/medic02/users/wbai/git/irtk_bin)
INCLUDE(${IRTK_DIR}/lib/IRTKInclude.cmake)

# Include file locations and linking libraries
INCLUDE_DIRECTORIES(./include ./external/libsvm-3.17)
LINK_LIBRARIES(vtkCommon vtkIO vtkGraphics /usr/lib/libgsl.a /usr/lib/libgslcblas.a /usr/lib/x86_64-linux-gnu/libpng.a  /usr/lib/x86_64-linux-gnu/libz.a)

# Target
ADD_EXECUTABLE(label_fusion label_fusion.cc external/libsvm-3.17/svm.cpp)
