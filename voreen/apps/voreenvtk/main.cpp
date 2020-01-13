/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "vtkVoreenActor.h"
#include "vtkActor.h" 

#include "vtkRenderWindow.h" 
#include "vtkRenderWindowInteractor.h" 
#include "vtkInteractorStyleTrackballCamera.h" 
#include "vtkRenderer.h" 
#include "vtkObjectFactory.h"
#include "vtkConeSource.h" 
#include "vtkPolyDataMapper.h" 
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkProperty.h"
#include "vtkImageReader2.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkVolumeTextureMapper2D.h"
#include "vtkVolumeTextureMapper3D.h"
#include "vtkVectorText.h"
#include "vtkLinearExtrusionFilter.h"
#include "vtkTextActor3D.h"
#include "vtkContourFilter.h"

#include "vtkVoreenInteractor.h"

#include "voreen/core/voreenapplication.h"

using namespace voreen;

void PerformVtkRendering(vtkRenderer* ren1, vtkImageReader2* rawReader, int ren_mode) {

    //0 = VolumeRayCastMapper  2 = TextureMapper2D  1 = TextureMapper3D
    int mode = ren_mode;  

    double max = rawReader->GetOutput()->GetScalarTypeMax();
    
    vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0,0.8,0.8,0.8);
    colorTransferFunction->AddRGBPoint(max,0.4,0.4,0.4);

    vtkPiecewiseFunction *opacityTransferFunction = vtkPiecewiseFunction::New();
    opacityTransferFunction->AddPoint(0, 0.0);
    opacityTransferFunction->AddPoint(max,1);

    // The property describes how the data will look
    vtkVolumeProperty *volumeProperty = vtkVolumeProperty::New();
    volumeProperty->SetColor(colorTransferFunction);
    volumeProperty->SetScalarOpacity(opacityTransferFunction);
    volumeProperty->SetAmbient( 0.0 );
    volumeProperty->SetDiffuse( 0.9 );
    volumeProperty->SetSpecular( 0.2 );
    volumeProperty->SetSpecularPower( 10.0 );
    volumeProperty->ShadeOn();

     vtkVolume *volData= vtkVolume::New();
     volData->SetProperty(volumeProperty);

     if (mode == 0) {
        vtkVolumeRayCastCompositeFunction* rayCompositeFunction = vtkVolumeRayCastCompositeFunction::New();
        rayCompositeFunction->SetCompositeMethodToClassifyFirst();
        vtkVolumeRayCastMapper *rayCastMapper = vtkVolumeRayCastMapper::New();
        rayCastMapper->SetInput( rawReader->GetOutput());
        rayCastMapper->SetVolumeRayCastFunction(rayCompositeFunction);      
        rayCastMapper->SetSampleDistance(0.5); 
        rayCastMapper->SetImageSampleDistance(0.5); 
        volData->SetMapper(rayCastMapper);
        
     }

     if (mode == 1) {
        vtkVolumeTextureMapper3D* volumeMapper = vtkVolumeTextureMapper3D::New();
        volumeMapper->SetInput(rawReader->GetOutput());
        volData->SetMapper(volumeMapper);
     }

     if (mode == 2) {
        vtkVolumeTextureMapper2D* volumeMapper = vtkVolumeTextureMapper2D::New();
        volumeMapper->SetInput(rawReader->GetOutput());
        volData->SetMapper(volumeMapper);
     }
    ren1->AddVolume(volData);
}

void CreateCone(vtkRenderer* ren1) {

    //Create a Cone
    vtkConeSource *cone = vtkConeSource::New();
    cone->SetHeight(1);
    cone->SetRadius(1.0);
    cone->SetResolution(10);
    vtkPolyDataMapper *coneMapper = vtkPolyDataMapper::New();
    coneMapper->SetInputConnection(cone->GetOutputPort());
    vtkActor *coneActor = vtkActor::New();
    coneActor->SetMapper(coneMapper);
    ren1->AddActor(coneActor);
}

int main(int argc, char** argv) {

    //Create a render window
    vtkRenderer *ren1 = vtkRenderer::New();
    vtkVoreenActor* voreenActor = 0;
    ren1->SetBackground(0,0,0);
    vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);
    renWin->SetSize(1024, 768);
    renWin->SetPosition(0,0);
    renWin->Render();

    //Create 2D-Text Overlay
    vtkTextActor* text = vtkTextActor::New();
    text->SetInput("vtkVoreen");
    vtkTextProperty* tprop = text->GetTextProperty();
    tprop->SetFontFamilyToArial();
    tprop->BoldOn();
    tprop->ShadowOn();
    tprop->SetLineSpacing(1.0);
    tprop->SetFontSize(30);
    tprop->SetColor(1,1,1);
    tprop->SetShadowOffset(2,2);
    text->SetDisplayPosition(10, 10);
    ren1->AddActor2D(text);

    //Import Raw Data from nucleon.raw
    /*
    vtkImageReader2* rawReader= vtkImageReader2::New();
    rawReader->SetFileName("../../data/nucleon.raw");
    rawReader->setDataExtent(0,40,0,40,0,40);
    rawReader->setDataScalarTypeToUnsignedChar();
    rawReader->setDataByteOrderToLittleEndian();
    rawReader->SetFileDimensionality(3);
    rawReader->setDataSpacing(1,1,1);
    rawReader->setDataOrigin(0,0,0);
    rawReader->SetNumberOfScalarComponents(1);
    rawReader->FileLowerLeftOn();
    rawReader->UpdateWholeExtent();
    */
    
    //Import Raw Data from walnut.raw
    
    /*vtkImageReader2* rawReader= vtkImageReader2::New();
    rawReader->SetFileName("../../data/walnut.raw");
    rawReader->setDataExtent(0,127,0,95,0,113);
    rawReader->setDataScalarTypeToUnsignedShort();
    rawReader->SetFileDimensionality(3);
    rawReader->setDataSpacing(1,1,1);
    rawReader->FileLowerLeftOn();
    rawReader->UpdateWholeExtent(); */

    //Use VTK-Renderer

    //PerformVtkRendering(ren1,rawReader,2);

    // or 

    //Create Voreen Actor
    voreenActor = vtkVoreenActor::New();
    voreenActor->connectToRenderer(ren1);
    if (!VoreenApplication::app()) {
        LERRORC("vtkSample.main", "VoreenApplication not initiated");
        return 1;
    }

    try {
        std::string workspace = "standard.vws";
        if (argc == 2) {
            workspace = std::string(argv[1]);
        }
        voreenActor->loadWorkspace(VoreenApplication::app()->getWorkspacePath(workspace));
    }
    catch (voreen::SerializationException& e) {
        LERRORC("vtkSample.main", e.what());
    }

    //voreenActor->SetInput(rawReader);
    //voreenActor->SetRendererToSingleVolumeRayCaster();
    //voreenActor->ResetTransferFunction("all");

    ren1->ResetCamera();
    renWin->Render();
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);


    if (voreenActor){
        vtkVoreenInteractor* voreen_style = vtkVoreenInteractor::New();
        voreen_style->SetVoreenSource(voreenActor);
        iren->SetInteractorStyle(voreen_style);
    }
    else {
         vtkInteractorStyleTrackballCamera* trackball_style = vtkInteractorStyleTrackballCamera::New();
         iren->SetInteractorStyle(trackball_style);
    }
        
    iren->Initialize();
    iren->Start();

    return 0;
}

