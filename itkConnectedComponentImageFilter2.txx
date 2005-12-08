/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConnectedComponentImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/26 16:17:55 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkConnectedComponentImageFilter2_txx
#define _itkConnectedComponentImageFilter2_txx

#include "itkConnectedComponentImageFilter2.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkEquivalencyTable.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkNeighborhoodAlgorithm.h"


namespace itk
{
template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter2< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // We need all the input.
  InputImagePointer input = const_cast<InputImageType *>(this->GetInput());
  
  input->SetRequestedRegion( input->GetLargestPossibleRegion() );
}



template <class TInputImage, class TOutputImage>
void 
ConnectedComponentImageFilter2<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()
    ->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter2< TInputImage, TOutputImage >
::GenerateData()
{


  // create an equivalency table
  EquivalencyTable::Pointer eqTable = EquivalencyTable::New();

//  OutputPixelType    label, originalLabel, neighborLabel;
  OutputPixelType    maxLabel = NumericTraits<OutputPixelType>::Zero;
  const OutputPixelType maxPossibleLabel=NumericTraits<OutputPixelType>::max();

  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();

  // Allocate the output and initialize to zeros
  this->AllocateOutputs();
  output->FillBuffer( NumericTraits<OutputPixelType>::Zero );
  


  // along with a neighborhood iterator on the input, use a standard
  // iterator on the input and output
  ImageRegionConstIterator<InputImageType> it;
  ImageRegionIterator<OutputImageType> oit;
  it = ImageRegionConstIterator<InputImageType>(input,
                                                output->GetRequestedRegion());
  oit = OutputIteratorType(output,
       	                   output->GetRequestedRegion());
  

  // progress reporter overhead is similar to the other parts of the 
  // filter. We don't want this to be a significant part of the filter
  // execution time, so pretend we have 1 pixel per pass of the
  // algorithm.
  ProgressReporter progress(this, 0, 3);

  // Mark the output image as either background or unlabeled
  it.GoToBegin();
  oit.GoToBegin();
  while (!it.IsAtEnd())
    {
    if (it.Get() != NumericTraits<InputPixelType>::Zero)
      {
      // mark pixel as unlabeled
      oit.Set(maxPossibleLabel);
      }
    
    ++it;
    ++oit;
    }

  progress.CompletedPixel();
  // iterate over the image, labeling the objects and defining
  // equivalence classes.  Use the neighborhood iterator to access the
  // "previous" neighbor pixels and an output iterator to access the
  // current pixel

  SizeType kernelRadius;
  kernelRadius.Fill(1);
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<TInputImage> FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  typename FaceCalculatorType::FaceListType faceList;
  faceList = faceCalculator(input, output->GetRequestedRegion(), kernelRadius);
  typename FaceCalculatorType::FaceListType::iterator fit;
  
  
  oit = OutputIteratorType(output, *(faceList.begin()));

  // Set up the boundary condition to be zero padded (used on output image)
  ConstantBoundaryCondition<TOutputImage> BC;
  BC.SetConstant(NumericTraits<OutputPixelType>::Zero);

  // Neighborhood iterator.  Let's use a shaped neighborhood so we can
  // restrict the access to face connected neighbors. This iterator
  // will be applied to the output image
  

  NeighborhoodIteratorType nit(kernelRadius, output,
                               *(faceList.begin()));
  //nit.OverrideBoundaryCondition(&BC); // assign the boundary condition

  setNeighborhood(nit);
  // do the region that doesn't need bounds checking
  // on the first region - bounds checking should be automatically
  // disabled, but lets be explicit
  //nit.NeedToUseBoundaryConditionOff();
  doRegion(eqTable, oit, nit, maxLabel);

  // now do the faces. We need to create a new neighborhood iterator that
  // looks forwards and backwards, because we are visiting the
  // faces last and we don't know whether a given face is in front
  // of or behind the body
  for (fit = (++faceList.begin()); fit != faceList.end(); ++fit)
    {
       oit = OutputIteratorType(output, *fit);
       NeighborhoodIteratorType nit2(kernelRadius,output, *fit);	
       nit2.OverrideBoundaryCondition(&BC); // assign the boundary condition
       setFaceNeighborhood(nit2);
       doRegion(eqTable, oit, nit2, maxLabel);
    }


  progress.CompletedPixel();

  // Flatten the equavalency table
  eqTable->Flatten();

  // remap the labels
  // reset the output iterator to reference the complete volume
  oit = OutputIteratorType(output,
       	                   output->GetRequestedRegion());

  oit.GoToBegin();
  while ( !oit.IsAtEnd() )
    {
    OutputPixelType label = oit.Get();
    // if pixel has a label, write out the final equivalence
    if (label != NumericTraits<OutputPixelType>::Zero)
      {
      oit.Set( eqTable->Lookup( label ) );
      }
    ++oit;
    }
  progress.CompletedPixel();

}

template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter2< TInputImage, TOutputImage >
::doRegion(EquivalencyTable::Pointer &eqTable, 
OutputIteratorType &oit, NeighborhoodIteratorType &nit, 
OutputPixelType &maxLabel)
{
  OutputPixelType    label, originalLabel, neighborLabel;
  const OutputPixelType maxPossibleLabel=NumericTraits<OutputPixelType>::max();
  nit.GoToBegin();
  oit.GoToBegin();
  while ( !oit.IsAtEnd() )
    {
    // Get the current pixel label
    label = oit.Get();
    originalLabel = label;

    // If the pixel is not background
    if (label != NumericTraits<OutputPixelType>::Zero)
      {
      nit += oit.GetIndex() - nit.GetIndex();
      // loop over the "previous" neighbors to find labels.  this loop
      // may establish one or more new equivalence classes
      typename NeighborhoodIteratorType::ConstIterator sIt;
      for (sIt = nit.Begin(); !sIt.IsAtEnd(); ++sIt)
        {
        // get the label of the pixel previous to this one along a
        // particular dimension (neighbors activated in neighborhood iterator)
        neighborLabel = sIt.Get();

        // if the previous pixel has a label, verify equivalence or
        // establish a new equivalence
        if (neighborLabel != NumericTraits<OutputPixelType>::Zero)
          {
          // if current pixel is unlabeled, copy the label from neighbor
          if (label == maxPossibleLabel)
            {
            // copy the label from the previous pixel
            label = neighborLabel;
            }
          // else if current pixel has a label that is not already
          // equivalent to the label of the previous pixel, then setup
          // a new equivalence.  
          else if ((label != neighborLabel)
                   && (eqTable->RecursiveLookup(label)
                       != eqTable->RecursiveLookup(neighborLabel))) 
            {
            eqTable->Add(label, neighborLabel);
            }
          }
        }

      // if none of the "previous" neighbors were set, then make a new label
      if (originalLabel == label)
        {
        // create a new entry label
        if (maxLabel == maxPossibleLabel)
          {
          itkWarningMacro(<< "ConnectedComponentImageFilter2::GenerateData: Number of labels exceeds number of available labels for the output type." );
          }
        else
          {
          ++maxLabel;
          }

        // assign the new label
        label = maxLabel;
        }

      // Finally, set the output pixel to whatever label we have
      if (label != originalLabel)
        {
        oit.Set( label );
        }
      }

    // move the iterators
    ++oit;
    }
}

template< class TInputImage, class TOutputImage >
void ConnectedComponentImageFilter2< TInputImage, TOutputImage >::
setFaceNeighborhood(NeighborhoodIteratorType &nit)
{
  typename NeighborhoodIteratorType::OffsetType offset;
  unsigned int d;
  if (!m_FullyConnected)
    {
    // only activate the neighbors that are face connected
    // to the current pixel. do not include the center pixel
    offset.Fill(0);
    for (d=0; d < InputImageType::ImageDimension; ++d)
      {
      offset[d] = -1;
      nit.ActivateOffset(offset);
      offset[d] = 1;
      nit.ActivateOffset(offset);
      offset[d] = 0;
      }
    }
  else
    {
    // activate all neighbors that are face+edge+vertex
    // connected to the current pixel. do not include the center pixel
    unsigned int centerIndex = nit.GetCenterNeighborhoodIndex();
    for (d=0; d < 2*centerIndex + 1; d++)
      {
      offset = nit.GetOffset(d);
      nit.ActivateOffset(offset);
      }
      nit.DeactivateOffset(nit.GetOffset(centerIndex));
    }

}

template< class TInputImage, class TOutputImage >
void ConnectedComponentImageFilter2< TInputImage, TOutputImage >::
setNeighborhood(NeighborhoodIteratorType &nit)
{
  // only activate the indices that are "previous" to the current
  // pixel and face connected (exclude the center pixel from the
  // neighborhood)
  //
  unsigned int d;
  typename NeighborhoodIteratorType::OffsetType offset;

  if (!m_FullyConnected)
    {
    // only activate the "previous" neighbors that are face connected
    // to the current pixel. do not include the center pixel
    offset.Fill(0);
    for (d=0; d < InputImageType::ImageDimension; ++d)
      {
      offset[d] = -1;
      nit.ActivateOffset(offset);
      offset[d] = 0;
      }
    }
  else
    {
    // activate all "previous" neighbors that are face+edge+vertex
    // connected to the current pixel. do not include the center pixel
    unsigned int centerIndex = nit.GetCenterNeighborhoodIndex();
    for (d=0; d < centerIndex; d++)
      {
      offset = nit.GetOffset(d);
      nit.ActivateOffset(offset);
      }
    }
}

template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter2< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
}

} // end namespace itk

#endif
