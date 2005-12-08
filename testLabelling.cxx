// test and time connected component labelling in itk

#include <stdlib.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkConnectedComponentImageFilter.h"
#include "itkConnectedComponentImageFilter1.h"
#include "itkConnectedComponentImageFilter2.h"
#include "itkConnectedComponentImageFilter3.h"

#include "itkBinaryThresholdImageFilter.h"

//const int dimensions = 2;

// 2d case
// ./testLabelling /home/richardb/Build/itk-mima2/Insight/Examples/Data/BrainMidSagittalSlice.png crap.tif 100

// 3d case
// ./testLabelling brain.tif crap.tif 100 3


template <int dimensions, class InputPixelType, class OutputPixelType>
void ThreshAndLabel(std::string infile, std::string outfile, int thresh,
		    int repeats, bool fullyConnected, int whichMethod,
		    int OutVal=255, int InVal=0)
{
  typedef itk::Image< InputPixelType,  dimensions >    InputImageType;
  typedef itk::Image< OutputPixelType,  dimensions >   OutputImageType;
  
  typedef itk::BinaryThresholdImageFilter< InputImageType, InputImageType >  ThresholdFilterType;

  typedef itk::ConnectedComponentImageFilter<InputImageType, OutputImageType> LabelFilt0;
  typedef itk::ConnectedComponentImageFilter1<InputImageType, OutputImageType> LabelFilt1;
  typedef itk::ConnectedComponentImageFilter2<InputImageType, OutputImageType> LabelFilt2;
  typedef itk::ConnectedComponentImageFilter3<InputImageType, OutputImageType> LabelFilt3;

  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  // the reader, writer and thresholding objects
  typename ReaderType::Pointer reader;
  typename WriterType::Pointer writer;
  typename ThresholdFilterType::Pointer threshfilter;

  reader = ReaderType::New();
  writer = WriterType::New();
  threshfilter = ThresholdFilterType::New();

  // Labelling objects
  typename LabelFilt0::Pointer labeller0;
  typename LabelFilt1::Pointer labeller1;
  typename LabelFilt2::Pointer labeller2;
  typename LabelFilt3::Pointer labeller3;

  labeller0 = LabelFilt0::New();
  labeller1 = LabelFilt1::New();
  labeller2 = LabelFilt2::New();
  labeller3 = LabelFilt3::New();

  // set up reader
  reader->SetFileName(infile.c_str());
  threshfilter->SetInput(reader->GetOutput());

  // set up the thresholder
  threshfilter->SetOutsideValue(OutVal);
  threshfilter->SetInsideValue(InVal);
  threshfilter->SetUpperThreshold( thresh );
  threshfilter->SetLowerThreshold( 0 );
  threshfilter->Update();

  std::cout << "Read file and thresholded. Now labelling " << repeats << " times" << std::endl;



  switch (whichMethod) {
  case 0:
    labeller0->SetFullyConnected(fullyConnected);
    for (int i=0;i<repeats;i++) 
      {
      // set first labeller
      labeller0->SetInput(NULL);
      labeller0->SetInput(threshfilter->GetOutput());
      labeller0->Update();
      }
    writer->SetInput(labeller0->GetOutput());
    break;
  case 1:
    labeller1->SetFullyConnected(fullyConnected);
    for (int i=0;i<repeats;i++) 
      {
      // set first labeller
      labeller1->SetInput(NULL);
      labeller1->SetInput(threshfilter->GetOutput());
      labeller1->Update();
      }
    writer->SetInput(labeller1->GetOutput());
    break;
  case 2:
    labeller2->SetFullyConnected(fullyConnected);
    for (int i=0;i<repeats;i++) 
      {
      // set first labeller
      labeller2->SetInput(NULL);
      labeller2->SetInput(threshfilter->GetOutput());
      labeller2->Update();
      }
    writer->SetInput(labeller2->GetOutput());
    break;
  case 3:
    labeller3->SetFullyConnected(fullyConnected);
    for (int i=0;i<repeats;i++) 
      {
      // set first labeller
      labeller3->SetInput(NULL);
      labeller3->SetInput(threshfilter->GetOutput());
      labeller3->Update();
      }
    writer->SetInput(labeller3->GetOutput());
    break;
  default:
    std::cerr << "Unknown labelling option" << std::endl;
  }

  writer->SetFileName(outfile.c_str());
  writer->Update();

}


int main( int argc, char * argv[] )
{
  if( argc < 8 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImageFile ";
    std::cerr << " outputImageFile threshold dimensions labellernumber[0-3] repeats connect[0,1]" << std::endl;
    return EXIT_FAILURE;
    }

  // Assume the input type is char
  typedef  unsigned char InputPixelType;
  typedef  unsigned short OutputPixelType;
  //typedef int OutputPixelType;

  std::string Infile = argv[1];
  std::string Outfile = argv[2];
  int thresh = atoi(argv[3]);
  int dimensions=atoi(argv[4]);
  int whichMethod=atoi(argv[5]);
  int repeats = atoi(argv[6]);
  bool connect = (bool)atoi(argv[7]);

  std::cout << "Labelling summary:" <<  std::endl;
  std::cout << "         source      : " << Infile << std::endl;
  std::cout << "         threshold   : " << thresh << std::endl;
  std::cout << "         dimensions  : " << dimensions << std::endl;
  std::cout << "         method      : " << whichMethod << std::endl;
  std::cout << "         repeats     : " << repeats << std::endl;
  std::cout << " fully connected     : " << connect << std::endl;

  switch(dimensions) {
  case 2:
    ThreshAndLabel<2, 
      InputPixelType, OutputPixelType>(Infile, 
				       Outfile, thresh,
				       repeats, connect, whichMethod);
    break;
  case 3:
    ThreshAndLabel<3, 
      InputPixelType, OutputPixelType>(Infile, 
				       Outfile, thresh,
				       repeats, connect, whichMethod);
    break;
  default:
    std::cerr << "Unsupported dimensions" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
