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

#ifndef __vtkVoreenActor_h
#define __vtkVoreenActor_h

#include <vector>

#include "vtkVoreenActorConfigure.h" // Include configuration header.
#include "vtkdummycanvas.h"
#include "vtkRenderer.h" 
#include "vtkProp.h"
#include "vtkImageData.h"
#include "vtkImageReader2.h"
#include "vtkCamera.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

#include "voreen/core/io/serialization/serialization.h"

class VTK_vtkVoreenActor_EXPORT vtkVoreenActor : public vtkProp {

public:

    vtkTypeRevisionMacro(vtkVoreenActor, vtkProp);
    static vtkVoreenActor* New();

    /**
     * Rendering method called by VTK
     */
    int RenderVolumetricGeometry(vtkViewport* view);
  
    /**
     * Initializes voreen application
     */
    void voreen_init();

    /**
     * Loads a workspace
     * @param workspace Path&Filename of the workspace to load
     */
    void loadWorkspace(std::string workspace);

    /**
     * Connects vtkVoreenActor to a renderer
     * @param ren vtkRenderer
     */
    void connectToRenderer(vtkRenderer* ren);

    /**
     * Synchronizes workspace cameras with vtk
     * @param cam vtkCamera
     */
    void syncCamera(vtkCamera* cam);

    /**
     * Getter for tgt canvas
     * @return tgt::DummyCanvas used by vreenVtk
     */
    voreen::VtkDummyCanvas* vtkVoreenActor::getCanvas();

    /**
     * Getter for current workspace name
     * @return Path and name of current workspace
     */
    std::string getWorkspaceName();

    /**
     * Check if voreen is initialized
     * @return voreen_initialized_
     */
    bool isInitialized();

    /**
     * Lists all processors in current network
     */
    void listProcessors();

    /**
     * Lists all Properties in of selected processor
     * @param procName Name of the processor
     */
    void listProperties(std::string procName);

    /**
     * Lists all Properties in of selected processor
     * @param procID ID of the processor
     */
    void listProperties(int procID);

    /**
     * Increases or decreases a custom int property
     * @param procName Processor name
     * @param propName Property name
     * @param increase Increase or decrease value
     */
    void changeIntPropByName(std::string procName, std::string propName, bool increase);
    
    /**
     * Changes a custom int property
     * @param procName Processor name
     * @param propName Property name
     * @param val New value
     * @param index 1: Set Value | 2: Set Stepping | 3: SetMinValue | 4: SetMaxValue
     */
    void changeIntPropByName(std::string procName, std::string propName, int val, int index);
    
    /**
     * Increases or decreases a custom float property
     * @param procName Processor name
     * @param propName Property name
     * @param increase Increase or decrease value
     */
    void changeFloatPropByName(std::string procName, std::string propName, bool increase);
    
    /**
     * Changes a custom float property
     * @param procName Processor name
     * @param propName Property name
     * @param val New value
     */
    void changeFloatPropByName(std::string procName, std::string propName, float val, int index);

    /**
     * Changes a custom string option property
     * @param procName Processor name
     * @param propName Property name
     * @param key New string value
     * @param index 1: Set Value | 2: Set Stepping | 3: SetMinValue | 4: SetMaxValue
     */
    void changeStringOptionPropByName(std::string procName, std::string propName, std::string key);
    
    /**
     * Changes a custom string option property
     * @param procName Processor name
     * @param propName Property name
     * @param key_id ID of the selection to be made from available keys
     */
    void changeStringOptionPropByName(std::string procName, std::string propName, int key_id);

    /**
     * Resets transfer function to two standard points
     * @param id index of the transfer function list (0 for first)
     */
    void resetTransferFunction(std::string procName);
    
    /**
     * Remove all keys from transfer function
     * @param id index of the transfer function list (0 for first)
     */
    void clearTransferFunction(std::string procName);
    
    /**
     * Loads a custom transfer function from .tfi file
     * @param procName Processor that has transfer function object
     * @param file Transfer function to load
     */
    void loadTransferFunction(std::string procName,std::string file);
    
    /**
     * Loads selects a transfer function key
     * @param id Id of the key that should be selected
     */
    void selectTransferfunctionKey(int id);
    
    /**
     * Moves the selected tf-key.
     * left/right will change the current intensity value,
     * up/down    will change the current opacity
     * @param direction Direction for the movement (0 1 2 3 -> left up right down)
     */
    void moveTransferFunctionKey(int direction);

    /**
     * Returns the number of the selected key.
     */
    int getSelectedKey();

    /**
     * Adds a new key to transfer function
     * @param procName Processor that has transfer function object
     * @param intensity Intensity [0..1]
     * @param red Red color (0..255)
     * @param green Green color (0..255)
     * @param blue Blue color (0..255)
     * @param opacity Opacity (0..255)
     */
    void addTransferFunctionKey(std::string procName, float intensity, float red, float green, float blue, float opacity);

    /**
     * Imports raw data from vtkImageReader
     * @param rawReader vtkImageReader object to get data from
     */
    void setInput(vtkImageReader2* rawReader);

    //Setters and Getters
    bool isSingleVolumeRaycaster();
    bool isSimpleRaycaster();
    bool isGlslRaycaster();

    void setRendererToSingleVolumeRayCaster();
    void setRendererToSimpleRayCaster();
    void setRendererToGlslRayCaster();

    void setCompositingDVR();
    void setCompositingMIP();
    void setCompositingISO();
    void setCompositingFHP();
    void setCompositingFHN();
    std::string getCompositeStyle();

    void setShading(int i);
    int getShading();

    std::string getRaycasterName();

    void statusToggle();

    void increaseSamplingRate();
    void decreaseSamplingRate();
    void setSamplingRate(float val);
    
    /**
     * Allows to select a clipping plane 
     * @param plane ("left", "right", "top", "bottom", "front", "back")
     */
    void selectClippingPlane(std::string plane);

    /**
     * Moves the selected clipping plane 
     * @param increase Defines if the value should be increased or decreased
     */
    void moveClippingPlane(bool increase);

    std::string getSelectedClippingPlane();

    void setClippingStepping(int value);
    void setClippingStepping(bool increase);
    int getClippingStepping();

    void enableBoundingBox();

    /**
     * Add message to output stream and display it
     * @param s Message string
     */
    void PrintMessage(std::string s);

protected:
    vtkVoreenActor();
    ~vtkVoreenActor();

    void printSelf(ostream& os, vtkIndent indent);

    /**
     * Creates the status output vtkRenderWindow and vtkTextActor
     * @param textWin vtkRenderWindow
     * @param status vtkTextActor
     */
    void createStatusOutput();

    /**
     * Sets the status output vtkRenderWindow and vtkTextActor
     * @param textWin vtkRenderWindow
     * @param status vtkTextActor
     */
    void setStatusOutput(vtkRenderer *textRen, vtkTextActor* status);

    /**
     * Print the current outputstream_ in status window
     */
    void StatusOutput();

private:
      bool voreen_initialized_;
      vtkVoreenActor(const vtkVoreenActor&);  // Not implemented.
      void operator=(const vtkVoreenActor&);  // Not implemented.
};

#endif
