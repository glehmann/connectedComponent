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

ADD_TEST(lab2d_a testLabelling BrainMidSagittalSlice.png orig.tif 100 2 0 1 0)
ADD_TEST(lab2d_b testLabelling BrainMidSagittalSlice.png meth1.tif 100 2 1 1 0)
ADD_TEST(lab2d_c testLabelling BrainMidSagittalSlice.png meth2.tif 100 2 2 1 0)
ADD_TEST(lab2d_d testLabelling BrainMidSagittalSlice.png meth3.tif 100 2 3 1 0)

ADD_TEST(comp2d_1 testEquiv orig.tif meth1.tif 2)
ADD_TEST(comp2d_2 testEquiv orig.tif meth2.tif 2)
ADD_TEST(comp2d_3 testEquiv orig.tif meth3.tif 2)

ADD_TEST(lab3d_a testLabelling ESCells.img orig.img 100 3 0 1 0)
ADD_TEST(lab3d_b testLabelling ESCells.img meth1.img 100 3 1 1 0)
ADD_TEST(lab3d_c testLabelling ESCells.img meth2.img 100 3 2 1 0)
ADD_TEST(lab3d_d testLabelling ESCells.img meth3.img 100 3 3 1 0)

ADD_TEST(comp3d_1 testEquiv orig.img meth1.img 3)
ADD_TEST(comp3d_2 testEquiv orig.img meth2.img 3)
ADD_TEST(comp3d_3 testEquiv orig.img meth3.img 3)

ADD_TEST(lab2d_af testLabelling BrainMidSagittalSlice.png orig.tif 100 2 0 1 1)
ADD_TEST(lab2d_bf testLabelling BrainMidSagittalSlice.png meth1.tif 100 2 1 1 1)
ADD_TEST(lab2d_cf testLabelling BrainMidSagittalSlice.png meth2.tif 100 2 2 1 1)
ADD_TEST(lab2d_df testLabelling BrainMidSagittalSlice.png meth3.tif 100 2 3 1 1)

## fully connected tests
ADD_TEST(comp2d_1f testEquiv orig.tif meth1.tif 2)
ADD_TEST(comp2d_2f testEquiv orig.tif meth2.tif 2)
ADD_TEST(comp2d_3f testEquiv orig.tif meth3.tif 2)

ADD_TEST(lab3d_af testLabelling ESCells.img orig.img 100 3 0 1 1)
ADD_TEST(lab3d_bf testLabelling ESCells.img meth1.img 100 3 1 1 1)
ADD_TEST(lab3d_cf testLabelling ESCells.img meth2.img 100 3 2 1 1)
ADD_TEST(lab3d_df testLabelling ESCells.img meth3.img 100 3 3 1 1)

ADD_TEST(comp2d_1f testEquiv orig.img meth1.img 3)
ADD_TEST(comp2d_2f testEquiv orig.img meth2.img 3)
ADD_TEST(comp2d_2f testEquiv orig.img meth3.img 3)

## timing test
ADD_TEST(timing2d runTimer2d.sh)
ADD_TEST(timing3d runTimer3d.sh)
#ADD_TEST(Dilate2DCompare ImageCompare bdilate2D.tif hdilate2D.tif)

