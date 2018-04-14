#ifndef __TUBNAV_GENIO_H
#define __TUBNAV_GENIO_H

#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include "../Core/tubNavContainer.h"
#include "tubNavVmtkBranchedReader.h"
#include "../../definitions.h"

#include <itkVTKImageIO.h>

namespace tubNav {

   static vtkSmartPointer<vtkPolyData> ReadPolyData(const char* filename){
       // Read all the data from the file
       vtkSmartPointer<vtkXMLPolyDataReader> reader =
         vtkSmartPointer<vtkXMLPolyDataReader>::New();
       reader->SetFileName(filename);
       reader->Update();
       vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
       return data;
   }



   static vtkSmartPointer<vtkPolyData> ReadSTL(const char* filename){
       // Read all the data from the file
       vtkSmartPointer<vtkSTLReader> reader =
         vtkSmartPointer<vtkSTLReader>::New();
       reader->SetFileName(filename);
       reader->Update();
       vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
       return data;
   }

   static void SavePolyData(const char* filename,   vtkSmartPointer<vtkPolyData> polydata){

       vtkSmartPointer<vtkXMLPolyDataWriter> writer =
               vtkSmartPointer<vtkXMLPolyDataWriter>::New();

       writer->SetInputData(polydata);
       writer->SetFileName(filename);
       writer->Update();
   }


   template<typename T,
            typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
   std::shared_ptr<Container<T>>  ReadData(const char* filename){

       // Check suffix, and call the appropiate reader..
        vtkSmartPointer<vtkPolyData> polydata = ReadPolyData(filename);
        BranchedCenterlinesReader reader;
        if ( reader.IsBranchingData(polydata)){
            auto container = reader.CreateContainer<T>(polydata);
            return container;
        }


        return nullptr;
   }


   static ImageType::Pointer ReadVolume( const char* path){
     typedef itk::ImageFileReader< ImageType > VtkReaderType;
     VtkReaderType::Pointer vtkReader = VtkReaderType::New();

     //vtkReader->SetImageIO(itk::MetaImageIO::New());
     vtkReader->SetImageIO(itk::VTKImageIO::New() );
     vtkReader->SetFileName(path);

     try{
         vtkReader->Update();
     }
      catch ( itk::ExceptionObject&e){
             std::cerr << "exception in file reader" << std::endl;
             std::cerr << path << std::endl;
             std::cerr << e << std::endl;
      }

     return vtkReader->GetOutput();
   }


   static WallShearStressMap::Pointer ReadWSSVolume( const char* path){
     typedef itk::ImageFileReader<WallShearStressMap > VtkReaderType;
     VtkReaderType::Pointer vtkReader = VtkReaderType::New();

     //vtkReader->SetImageIO(itk::MetaImageIO::New());
     vtkReader->SetImageIO(itk::VTKImageIO::New() );
     vtkReader->SetFileName(path);

     try{
         vtkReader->Update();
     }
      catch ( itk::ExceptionObject&e){
             std::cerr << "exception in file reader" << std::endl;
             std::cerr << path << std::endl;
             std::cerr << e << std::endl;
      }

     return vtkReader->GetOutput();
   }

}

#endif
