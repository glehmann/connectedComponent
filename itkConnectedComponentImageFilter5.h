
#ifndef __itkConnectedComponentImageFilter5_h
#define __itkConnectedComponentImageFilter5_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include <vector>
#include <map>
#include "itkProgressReporter.h"

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
class ITK_EXPORT ConnectedComponentImageFilter5 : 
    public ImageToImageFilter< TInputImage, TOutputImage > 
{
public:
  /**
   * Standard "Self" & Superclass typedef.
   */
  typedef ConnectedComponentImageFilter5 Self;
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
  itkTypeMacro(ConnectedComponentImageFilter5, ImageToImageFilter);
  
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
  
  // only set after completion
  itkGetConstReferenceMacro(ObjectCount, long);

  // experimental
  itkSetMacro(AreaLambda, long);
  itkGetMacro(AreaLambda, long);

  itkSetMacro(OutValue, OutputPixelType);
  itkGetMacro(OutValue, OutputPixelType);

  itkSetMacro(SelectMax, bool);
  itkGetMacro(SelectMax, bool);
  

protected:
  ConnectedComponentImageFilter5() 
    {
    m_FullyConnected = false;
    m_ObjectCount = -1;
    m_AreaLambda = -1;
    m_SelectMax = false;
    m_OutValue = 1;
    }
  virtual ~ConnectedComponentImageFilter5() {}
  ConnectedComponentImageFilter5(const Self&) {}
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
  long m_ObjectCount;
  bool m_SelectMax;
  long m_AreaLambda;
  OutputPixelType m_OutValue, m_BiggestObjectLab;

  // some additional types
  typedef typename TOutputImage::RegionType::SizeType OutSizeType;

  // types to support the run length encoding of lines
  typedef class runLength
    {
    public:
      long int length;  // run length information - may be a more type safe way of doing this
      typename InputImageType::IndexType where;  // Index of the start of the run
      long int label; // the initial label of the run
    };

  typedef std::vector<runLength> lineEncoding;

  // the map storing lines
  typedef std::map<long, lineEncoding> LineMapType;
  
  typedef std::vector<long> OffsetVec;

  // the types to support union-find operations
  typedef std::vector<long int> UnionFindType;
  UnionFindType m_UnionFind;
  UnionFindType m_Consecutive;
  // functions to support union-find operations
  void InitUnion(const long int size) 
  {
    m_UnionFind = UnionFindType(size + 1);
  }
  void InsertSet(const long int label);
  long int LookupSet(const long int label);
  void LinkLabels(const long int lab1, const long int lab2);
  long int CreateConsecutive();
  //////////////////

  typedef std::vector<long> AttributeType;
  AttributeType m_Attributes;

  void CompareLines(lineEncoding &current, const lineEncoding &Neighbour);

  void FillOutput(const LineMapType &LineMap,
		  ProgressReporter &progress);

  void SetupLineOffsets(OffsetVec &LineOffsets);

  void computeAttributes(const long int ObjCount, const LineMapType &LineMap);
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConnectedComponentImageFilter5.txx"
#endif

#endif