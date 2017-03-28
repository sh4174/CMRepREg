#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

int main( int argc, char* argv[] )
{
  if( argc != 4 )
    {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputFileName> <Label>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;

  typedef unsigned char                      PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  const char * inputFileName = argv[1];
  const char * outputFileName = argv[2];
  const PixelType label = static_cast< PixelType >( atoi( argv[3] ) );

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName );

  typedef itk::LabelObject< PixelType, Dimension >  LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >          LabelMapType;

  typedef itk::LabelImageToLabelMapFilter< ImageType, LabelMapType > LabelImageToLabelMapFilterType;
  LabelImageToLabelMapFilterType::Pointer labelMapConverter = LabelImageToLabelMapFilterType::New();
  labelMapConverter->SetInput( reader->GetOutput() );
  labelMapConverter->SetBackgroundValue( itk::NumericTraits< PixelType >::Zero );

  typedef itk::LabelSelectionLabelMapFilter< LabelMapType > SelectorType;
  SelectorType::Pointer selector = SelectorType::New();
  selector->SetInput( labelMapConverter->GetOutput() );
  selector->SetLabel( label );

  typedef itk::LabelMapToLabelImageFilter< LabelMapType, ImageType > LabelMapToLabelImageFilterType;
  LabelMapToLabelImageFilterType::Pointer labelImageConverter = LabelMapToLabelImageFilterType::New();
  labelImageConverter->SetInput( selector->GetOutput( 0 ) );

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputFileName );
    writer->SetInput( labelImageConverter->GetOutput() );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      std::cerr << "Error: " << error << std::endl;
      return EXIT_FAILURE;
      }
    

  return EXIT_SUCCESS;
}
