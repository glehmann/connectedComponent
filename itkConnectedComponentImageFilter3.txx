/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConnectedComponentImageFilter2.txx,v $
  Language:  C++
  Date:      $Date: 2005/05/26 16:17:55 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkConnectedComponentImageFilter3_txx
#define _itkConnectedComponentImageFilter3_txx

#include "itkConnectedComponentImageFilter3.h"
#include "itkNumericTraits.h"
#include "itkEquivalencyTable.h"
// don't think we need the indexed version as we only compute the
// index at the start of each run, but there isn't a choice
#include "itkImageLinearConstIteratorWithIndex.h"  
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter3< TInputImage, TOutputImage >
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
ConnectedComponentImageFilter3<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()
    ->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter3< TInputImage, TOutputImage >
::GenerateData()
{
  // create an equivalency table
  EquivalencyTable::Pointer eqTable = EquivalencyTable::New();
  // create a line iterator
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputLineIteratorType;

  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();
  long LineIdx = 0;
  InputLineIteratorType inLineIt(input, output->GetRequestedRegion());
  inLineIt.SetDirection(0);
  LineMapType LineMap;
  // Allocate the output
  this->AllocateOutputs();

  OutputPixelType lab = NumericTraits<OutputPixelType>::Zero;
  OffsetVec LineOffsets;

  // set the progress reporter to deal with the number of lines
  long pixelcount = this->GetOutput()->GetRequestedRegion().GetNumberOfPixels();
  long xsize = this->GetOutput()->GetRequestedRegion().GetSize()[0];
  long linecount = pixelcount/xsize;
  ProgressReporter progress(this, 0, linecount/2);

  SetupLineOffsets(LineOffsets);

  for (inLineIt.GoToBegin(); ! inLineIt.IsAtEnd(); inLineIt.NextLine(), ++LineIdx)
    {
    inLineIt.GoToBeginOfLine();
    lineEncoding ThisLine;
    while (! inLineIt.IsAtEndOfLine())
      {
      InputPixelType PVal = inLineIt.Get();
      //std::cout << inLineIt.GetIndex() << std::endl;
      if (PVal != NumericTraits<InputPixelType>::Zero)
	{
	// We've hit the start of a run
	runLength thisRun;
	long length=0;
	typename InputImageType::IndexType thisIndex;
	++lab;
	thisIndex = inLineIt.GetIndex();
	//std::cout << thisIndex << std::endl;
	while ((PVal != NumericTraits<InputPixelType>::Zero) &&
	       (! inLineIt.IsAtEndOfLine()))
	  {
	  ++length;
	  ++inLineIt;
	  PVal = inLineIt.Get();
	  }
	// create the run length object to go in the vector
	thisRun.length=length;
	thisRun.label=lab;
	thisRun.where = thisIndex;
	thisRun.relabelled = false;
	ThisLine.push_back(thisRun);
	//std::cout << thisIndex[0] << " " << thisIndex[1] << " " << length << std::endl;
	}
      else 
	{
	++inLineIt;
	}
      }
    if (ThisLine.size() != 0)
      {
      // There are some runs on this line, so insert it into the map
      LineMap[LineIdx] = ThisLine;
      }
    progress.CompletedPixel();
    }


  // now process the map and make appropriate entries in an equivalence
  // table


  typename LineMapType::iterator MapBegin, MapEnd, LineIt;

  MapBegin = LineMap.begin();
  MapEnd = LineMap.end(); 
  LineIt = MapBegin;

  //while( LineIt != MapEnd)
  for (LineIt = MapBegin; LineIt != MapEnd; ++LineIt)
    {
    //lineEncoding L = LineIt->second;
    long ThisIdx = LineIt->first;
    //std::cout << "Line number = " << LineIt->first << std::endl;
    for (OffsetVec::const_iterator I = LineOffsets.begin();
	 I != LineOffsets.end(); ++I)
      {
      long NeighIdx = ThisIdx + (*I);
      // check if the neighbor is in the map
      typename LineMapType::const_iterator NN = LineMap.find(NeighIdx);
      if (NN != MapEnd) 
	{
	// Compare the two lines
	CompareLines(LineIt->second, NN->second, eqTable);
	}
      }
    }
  
  // Flatten the equavalency table
  eqTable->Flatten();
  // create the output

#if 0
  output->FillBuffer( NumericTraits<OutputPixelType>::Zero );
  ImageRegionIterator<OutputImageType> oit(output,
					   output->GetRequestedRegion());

  for (LineIt = MapBegin; LineIt != MapEnd; ++LineIt)
    {
    //lineEncoding *L = &(LineIt->second);
    typename lineEncoding::const_iterator cIt;
    for (cIt = LineIt->second.begin();cIt != LineIt->second.end();++cIt)
      {
      //runLength cL = *cIt;
      OutputPixelType lab = eqTable->Lookup( cIt->label);
      oit += (cIt->where) - oit.GetIndex();
      for (long i = 0; i < cIt->length; ++i, ++oit)
	{
	oit.Set(lab);
	}
      }
    }
#else
  // A more complex version that is intended to minimize the number of
  // visits to the output image which should improve cache
  // performance on large images. We also want to optimize the
  // performance of the map by being able to iterate through it,
  // rather than do lots of look ups. Don't know whether that will
  // make much difference in practice.
  // Note - this is unnecessary if AllocateOutputs initalizes to zero

  FillOutput(LineMap, eqTable, progress);

#endif
}


template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter3< TInputImage, TOutputImage >
::SetupLineOffsets(OffsetVec &LineOffsets)
{
  // Create a neighborhood so that we can generate a table of offsets
  // to "previous" line indexes
  // We are going to mis-use the neighborhood iterators to compute the
  // offset for us. All this messing around produces an array of
  // offsets that will be used to index the map
  typename TOutputImage::Pointer output = this->GetOutput();
  typedef Image<long, TOutputImage::ImageDimension - 1> PretendImageType;
  typedef typename PretendImageType::RegionType::SizeType PretendSizeType;
  typedef ConstShapedNeighborhoodIterator<PretendImageType> LineNeighborhoodType;

  typename PretendImageType::Pointer fakeImage;
  fakeImage = PretendImageType::New();

  typename PretendImageType::RegionType LineRegion;
  //LineRegion = PretendImageType::RegionType::New();

  OutSizeType OutSize;
  OutSize = output->GetRequestedRegion().GetSize();

  PretendSizeType PretendSize;
  PretendSize = fakeImage->GetRequestedRegion().GetSize();
  // The first dimension has been collapsed
  for (unsigned int i = 0; i<PretendSize.GetSizeDimension(); i++)
    {
    PretendSize[i] = OutSize[i+1];
    }

  LineRegion.SetSize(PretendSize);
  PretendSizeType kernelRadius;
  kernelRadius.Fill(1);
  LineNeighborhoodType lnit(kernelRadius, fakeImage, LineRegion);

  // only activate the indices that are "previous" to the current
  // pixel and face connected (exclude the center pixel from the
  // neighborhood)
  //
  unsigned int d;
  typename LineNeighborhoodType::OffsetType offset;

  if (!m_FullyConnected)
    {
    // only activate the "previous" neighbors that are face connected
    // to the current pixel. do not include the center pixel
    offset.Fill(0);
    for (d=0; d < PretendImageType::ImageDimension; ++d)
      {
      offset[d] = -1;
      lnit.ActivateOffset(offset);
      offset[d] = 0;
      }
    }
  else
    {
    // activate all "previous" neighbors that are face+edge+vertex
    // connected to the current pixel. do not include the center pixel
    unsigned int centerIndex = lnit.GetCenterNeighborhoodIndex();
    for (d=0; d < centerIndex; d++)
      {
      offset = lnit.GetOffset(d);
      lnit.ActivateOffset(offset);
      }
    }

  typename LineNeighborhoodType::IndexListType ActiveIndexes;
  ActiveIndexes = lnit.GetActiveIndexList();

  typedef std::vector<long> OffsetVec;

  typename LineNeighborhoodType::IndexListType::const_iterator LI;
  
  PretendSizeType SizeBuff;
  SizeBuff[0] = 1;
  // should do this with a copy
  for (int i =0; i < PretendImageType::ImageDimension - 1; i++)
    {
    SizeBuff[i+1] = PretendSize[i];
    }
  unsigned int pos = 0;
  for (LI=ActiveIndexes.begin(); LI != ActiveIndexes.end(); LI++, pos++)
    {
    unsigned int idx = *LI;
    offset = lnit.GetOffset(idx);
    //std::cout << offset << std::endl;
    int vv = 0;
    for (int J = 0; J < PretendImageType::ImageDimension;J++)
      {
      vv += offset[J] * SizeBuff[J];
      }
    //std::cout << vv << std::endl;
    LineOffsets.push_back( vv);
    }
  
  // LineOffsets is the thing we wanted.
}
template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter3< TInputImage, TOutputImage >
::CompareLines(lineEncoding &current, const lineEncoding &Neighbour,
	       EquivalencyTable::Pointer &eqTable)
{
  long offset = 0;
  if (m_FullyConnected)
    offset = 1;

  typename lineEncoding::const_iterator nIt, mIt;
  typename lineEncoding::iterator cIt;

  mIt = Neighbour.begin(); // out marker iterator

  for (cIt = current.begin();cIt != current.end();++cIt)
    {
    //runLength cL = *cIt;
    long cStart = cIt->where[0];  // the start x position
    long cLast = cStart + cIt->length - 1;

    for (nIt=mIt; nIt != Neighbour.end(); ++nIt)
      {
      //runLength nL = *nIt;
      long nStart = nIt->where[0];
      long nLast = nStart + nIt->length - 1;
      // there are a few ways that neighbouring lines might overlap
      //   neighbor      S                  E
      //   current    S                        E
      //------------------------------------------
      //   neighbor      S                  E
      //   current    S                E
      //------------------------------------------
      //   neighbor      S                  E
      //   current             S                  E
      //------------------------------------------
      //   neighbor      S                  E
      //   current             S       E
      //------------------------------------------
      long ss1 = nStart - offset;
      long ss2 = nStart + offset;
      long ee1 = nLast - offset;
      long ee2 = nLast + offset;
      bool eq = false;
      // the logic here can probably be improved a lot
      if ((ss1 >= cStart) && (ee2 <= cLast))
	{
	// case 1
	eq = true;
	} 
      else 
	{
	if ((ss1 <= cLast) && (ee2 >= cLast))
	  {
	  // case 2
	  eq = true;
	  }
	else 
	  {
	  if ((ss1 <= cStart) && (ee2 >= cStart))
	    {
	    // case 3 
	    eq = true;
	    }					
	  else 
	    {
	    if ((ss1 <= cStart) && (ee2 >= cLast))
	      {
	      // case 4
	      eq = true;
	      }
	    }
	  }
	}
      if (eq) 
	{
        // this one has already been relabelled, so bung it straight
        // into the table
	if (cIt->relabelled) 
	  {
	  eqTable->Add(cIt->label, nIt->label);
	  }
	else
	  {
	  // Handle some equivalency here to help reduce the load
	  // on the equivalency table
	  cIt->label = nIt->label;
	  cIt->relabelled = true;
	  }
	} 

      if (ee1 >= cLast)
	{
	// No point looking for more overlaps with the current run
	// because the neighbor run is either case 2 or 4
	mIt = nIt;
	break;
	}
      }
    }

}

template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter3< TInputImage, TOutputImage >
::FillOutput(const LineMapType &LineMap,
	     const EquivalencyTable::Pointer &eqTable,
	     ProgressReporter &progress)
{

  typename LineMapType::const_iterator MapBegin, MapEnd, LineIt;
  typename TOutputImage::Pointer output = this->GetOutput();
  MapBegin = LineMap.begin();
  MapEnd = LineMap.end(); 
  LineIt = MapBegin;

  ImageRegionIterator<OutputImageType> oit(output,
					   output->GetRequestedRegion());

  ImageRegionIterator<OutputImageType> fstart=oit, fend=oit;
  fstart.GoToBegin();
  fend.GoToEnd();

  for (LineIt = MapBegin; LineIt != MapEnd; ++LineIt)
    {
    // now fill the labelled sections
    typename lineEncoding::const_iterator cIt;

    //std::cout << LineIt->first << std::endl;

    for (cIt = LineIt->second.begin();cIt != LineIt->second.end();++cIt)
      {
      //runLength cL = *cIt;
      OutputPixelType lab = eqTable->Lookup( cIt->label);
      oit.SetIndex(cIt->where);
      // initialize the non labelled pixels
      for (; fstart != oit; ++fstart)
	{
	fstart.Set(NumericTraits<OutputPixelType>::Zero );
	}
      for (long i = 0; i < cIt->length; ++i, ++oit)
	{
	oit.Set(lab);
	}
      fstart = oit;
      //++fstart;
      }
    progress.CompletedPixel();
    }
  // fill the rest of the image with zeros
  for (; fstart != fend; ++fstart)
    {
    fstart.Set(NumericTraits<OutputPixelType>::Zero );
    }
}

template< class TInputImage, class TOutputImage >
void
ConnectedComponentImageFilter3< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
}

} // end namespace itk

#endif
