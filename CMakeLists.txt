#Change PROJECT_NAME to the name of your project
PROJECT(connectedComponent)

#include some macros from another file...
INCLUDE(IJMacros.txt)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
SET(Required_Packages
  ITK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
  ITKCommon
  ITKBasicFilters
  ITKIO
)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project
SET(CurrentExe "testEquiv")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "testLabelling")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})


#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.

ADD_TEST(lab2d_a testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_a.png 100 2 0 1 0)
ADD_TEST(lab2d_b testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_b.png 100 2 1 1 0)
ADD_TEST(lab2d_c testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_c.png 100 2 2 1 0)
ADD_TEST(lab2d_d testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_d.png 100 2 3 1 0)

ADD_TEST(comp2d_1 testEquiv lab2d_a.png lab2d_b.png 2)
ADD_TEST(comp2d_2 testEquiv lab2d_a.png lab2d_c.png 2)
ADD_TEST(comp2d_3 testEquiv lab2d_a.png lab2d_d.png 2)

ADD_TEST(lab3d_a testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_a.img 100 3 0 1 0)
ADD_TEST(lab3d_b testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_b.img 100 3 1 1 0)
ADD_TEST(lab3d_c testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_c.img 100 3 2 1 0)
ADD_TEST(lab3d_d testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_d.img 100 3 3 1 0)

ADD_TEST(comp3d_1 testEquiv lab3d_a.img lab3d_b.img 3)
ADD_TEST(comp3d_2 testEquiv lab3d_a.img lab3d_c.img 3)
ADD_TEST(comp3d_3 testEquiv lab3d_a.img lab3d_d.img 3)

## fully connected tests
ADD_TEST(lab2d_af testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_af.png 100 2 0 1 1)
ADD_TEST(lab2d_bf testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_bf.png 100 2 1 1 1)
ADD_TEST(lab2d_cf testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_cf.png 100 2 2 1 1)
ADD_TEST(lab2d_df testLabelling ${CMAKE_SOURCE_DIR}/images/BrainMidSagittalSlice.png lab2d_df.png 100 2 3 1 1)

ADD_TEST(comp2d_1f testEquiv lab2d_af.png lab2d_bf.png 2)
ADD_TEST(comp2d_2f testEquiv lab2d_af.png lab2d_cf.png 2)
ADD_TEST(comp2d_3f testEquiv lab2d_af.png lab2d_df.png 2)

ADD_TEST(lab3d_af testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_af.img 100 3 0 1 1)
ADD_TEST(lab3d_bf testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_bf.img 100 3 1 1 1)
ADD_TEST(lab3d_cf testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_cf.img 100 3 2 1 1)
ADD_TEST(lab3d_df testLabelling ${CMAKE_SOURCE_DIR}/images/ESCells.img lab3d_df.img 100 3 3 1 1)

ADD_TEST(comp3d_1f testEquiv lab3d_af.img lab3d_bf.img 3)
ADD_TEST(comp3d_2f testEquiv lab3d_af.img lab3d_cf.img 3)
ADD_TEST(comp3d_2f testEquiv lab3d_af.img lab3d_df.img 3)

## timing test
# ADD_TEST(timing2d sh ${CMAKE_SOURCE_DIR}/runTimer2d.sh)
# ADD_TEST(timing3d sh ${CMAKE_SOURCE_DIR}/runTimer3d.sh)
#ADD_TEST(Dilate2DCompare ImageCompare bdilate2D.tif hdilate2D.tif)

