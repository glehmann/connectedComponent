/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConnectedComponentImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2004/11/04 20:40:37 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkConnectedComponentImageFilter3_h
#define __itkConnectedComponentImageFilter3_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include <vector>
#include <map>

namespace itk
{

/**
 * \class ConnectedComponentImageFilter
 * \brief Label the objects in a binary image
 *
 * ConnectedComponentImageFilter labels the objects in a binary image.
 * Each distinct object is assigned a unique label. The filter experiments
 * with some improvements to the existing implementation, and is based on
 * run length encoding along lines
 * The final object labels are in no particular order (and some object
 * labels may not be used on the final objects).  You can reorder the
 * labels such that object labels are consecutive and sorted based on
 * object size by passing the output of this filter to a
 * RelabelComponentImageFilter. 
 *
 * \sa ImageToImageFilter
 */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ConnectedComponentImageFilter3 : 
    public ImageToImageFilter< TInputImage, TOutputImage > 
{
public:
  /**
   * Standard "Self" & Superclass typedef.
   */
  typedef ConnectedComponentImageFilter3 Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;

  /**
   * Types from the Superclass
   */
  typedef typename Superclass::InputImagePointer InputImagePointer;

  /**
   * Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same.
   */
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputInternalPixelType;
  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename TInputImage::InternalPixelType InputInternalPixelType;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  
  /**
   * Image typedef support
   */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  typedef   typename TInputImage::IndexType       IndexType;
  typedef   typename TInputImage::SizeType        SizeType;
  typedef   typename TOutputImage::RegionType     RegionType;
  typedef   std::list<IndexType>                  ListType;

  /** 
   * Smart pointer typedef support 
   */
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /**
   * Run-time type information (and related methods)
   */
  itkTypeMacro(ConnectedComponentImageFilter3, ImageToImageFilter);
  
  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  /**
   * Set/Get whether the connected components are defined strictly by
   * face connectivity or by face+edge+vertex connectivity.  Default is
   * FullyConnectedOff.  For objects that are 1 pixel wide, use
   * FullyConnectedOn.
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);

protected:
  ConnectedComponentImageFilter3() 
    {
    m_FullyConnected = false;
    }
  virtual ~ConnectedComponentImageFilter3() {}
  ConnectedComponentImageFilter3(const Self&) {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /**
   * Standard pipeline method. 
   */
  void GenerateData();

  /** ConnectedComponentImageFilter needs the entire input. Therefore
   * it must provide an implementation GenerateInputRequestedRegion().
   * \sa ProcessObject::GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** ConnectedComponentImageFilter will produce all of the output.
   * Therefore it must provide an implementation of
   * EnlargeOutputRequestedRegion().
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));
  
private:

  bool m_FullyConnected;

  // some additional types
  typedef typename TOutputImage::RegionType::SizeType OutSizeType;

  // types to support the run length encoding of lines
  typedef class runLength
    {
    public:
      long length;  // run length information - may be a more type safe way of doing this
      typename InputImageType::IndexType where;  // Index of the start of the run
      OutputPixelType label; // the initial label of the run
      // a flag to indicate whether we have already relablled this entry
      bool relabelled;   
    };

  typedef std::vector<runLength> lineEncoding;

  // the map storing lines
  typedef std::map<long, lineEncoding> LineMapType;
  
  typedef std::vector<long> OffsetVec;

  void CompareLines(lineEncoding &current, const lineEncoding &Neighbour,
		    EquivalencyTable::Pointer &eqTable );

  void FillOutput(const LineMapType &LineMap,
		  const EquivalencyTable::Pointer &eqTable);

  void SetupLineOffsets(OffsetVec &LineOffsets);
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConnectedComponentImageFilter3.txx"
#endif

#endif
