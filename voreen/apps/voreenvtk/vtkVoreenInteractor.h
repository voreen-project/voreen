/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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
#include "vtkObjectFactory.h"
#include "vtkInteractorStyleTrackballCamera.h" 
#include "vtkRenderWindow.h" 
#include "vtkRenderWindowInteractor.h" 

#include "tgt/timer.h"
#include "tgt/glcanvas.h"
#include "tgt/event/keyevent.h"
#include "voreen/core/network/networkevaluator.h"


// Define interaction style
class vtkVoreenInteractor : public vtkInteractorStyleTrackballCamera
{
  public:
    static vtkVoreenInteractor* New();
    vtkTypeRevisionMacro(vtkVoreenInteractor, vtkInteractorStyleTrackballCamera);
	vtkVoreenActor* voreensource;
	tgt::MouseEvent::MouseButtons holdButton_;
	static const bool debug = 0;

	virtual tgt::KeyEvent::KeyCode getKeyCode(const int& key) {
		// Attention!!! This resides on tgt internal (ascii at the moment) settings which could change.
		// Normaly it should be handled in a switch statement.
		// (taken from glutcanvas.cpp)
		if ((key >= 91) && (key <= 127))   /// small letters, some brackets, etc.
			return static_cast<tgt::KeyEvent::KeyCode> (key);
		if ((key >= 65) && (key <= 90))    /// there are no UPPERCASE keys in tgt
			return static_cast<tgt::KeyEvent::KeyCode> (key+32);
		if ((key >= 38) && (key <= 64))    /// numbers, symbols
			return static_cast<tgt::KeyEvent::KeyCode> (key);
		if ((key >= 32) && (key <= 36))    /// symbols
			return static_cast<tgt::KeyEvent::KeyCode> (key);
		// the rest
		switch (key) {
		case   8: return tgt::KeyEvent::K_BACKSPACE;
		case   9: return tgt::KeyEvent::K_TAB;
		case  12: return tgt::KeyEvent::K_CLEAR;
		case  13: return tgt::KeyEvent::K_RETURN;
		case  19: return tgt::KeyEvent::K_PAUSE;
		case  27: return tgt::KeyEvent::K_ESCAPE;
		}
		return tgt::KeyEvent::K_UNKNOWN;
}

	virtual tgt::KeyEvent::KeyCode getSpecialKeyCode(const std::string& special) {
	    
			if (special.compare("F1") == 0) return tgt::KeyEvent::K_F1;
			if (special.compare("F2") == 0) return tgt::KeyEvent::K_F2;
			if (special.compare("F3") == 0) return tgt::KeyEvent::K_F3;
			if (special.compare("F4") == 0) return tgt::KeyEvent::K_F4;
			if (special.compare("F5") == 0) return tgt::KeyEvent::K_F5;
			if (special.compare("F6") == 0) return tgt::KeyEvent::K_F6;
			if (special.compare("F7") == 0) return tgt::KeyEvent::K_F7;
			if (special.compare("F8") == 0) return tgt::KeyEvent::K_F8;
			if (special.compare("F9") == 0) return tgt::KeyEvent::K_F9;
			if (special.compare("F10") == 0) return tgt::KeyEvent::K_F10;
			if (special.compare("F11") == 0) return tgt::KeyEvent::K_F11;
			if (special.compare("F12") == 0) return tgt::KeyEvent::K_F12;

			if (special.compare("Left") == 0) return tgt::KeyEvent::K_LEFT;
			if (special.compare("Up") == 0) return tgt::KeyEvent::K_UP;
			if (special.compare("Right") == 0) return tgt::KeyEvent::K_RIGHT;
			if (special.compare("Down") == 0) return tgt::KeyEvent::K_DOWN;

			if (special.compare("Delete") == 0) return tgt::KeyEvent::K_DELETE;
			if (special.compare("End") == 0) return tgt::KeyEvent::K_END;
			if (special.compare("Next") == 0) return tgt::KeyEvent::K_PAGEDOWN;

			if (special.compare("Insert") == 0) return tgt::KeyEvent::K_INSERT;
			if (special.compare("Home") == 0) return tgt::KeyEvent::K_HOME;
			if (special.compare("Prior") == 0) return tgt::KeyEvent::K_PAGEUP;

			return tgt::KeyEvent::K_UNKNOWN;
			
		
	}

	virtual int getModifier() {
    
		vtkRenderWindowInteractor *rwi = this->Interactor;	

		int result = 0;
		bool shift = rwi->GetShiftKey();
		bool alt = rwi->GetAltKey();
		bool ctrl = rwi->GetControlKey();

		if (debug){
			cout << "SHIFT-Modifier: " << shift << " ALT-Modifier: " << alt << " STRG-Modifier: " << ctrl << endl;
		}

		if (shift)
			result |= tgt::Event::SHIFT;
		if (alt)
			result |= tgt::Event::ALT;
		if (ctrl)
			result |= tgt::Event::CTRL;
		return result;
}

	virtual int* GetMousePosition()
	{
		//the voreen y-axis 0 point starts at the top, in vtk on the bottom
		int* mouse_coords = this->GetInteractor()->GetEventPosition();
		int win_height = voreensource->getCanvas()->getHeight();
		int win_width = voreensource->getCanvas()->getWidth();
		mouse_coords[1] = win_height - mouse_coords[1];
		return mouse_coords;
	}

	virtual void SetVoreenSource(vtkVoreenActor* voreen)
	{
		voreensource = voreen;
        holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_NONE;
		
	}

	virtual void OnMouseMove () 
	{
	  if (voreensource->isInitialized()){
		  int* mouse_coords = GetMousePosition();
		  if (debug) cout << "Moved mouse to." << mouse_coords[0] << " " << mouse_coords[1] << endl;
		  int win_height = voreensource->getCanvas()->getHeight();
		  int win_width = voreensource->getCanvas()->getWidth();
          tgt::MouseEvent::MouseAction action = tgt::MouseEvent::MOTION;
          tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
		  tgt::MouseEvent* moveEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1], action, tgtModifier, holdButton_ ,tgt::ivec2(win_width,win_height));
		  voreensource->getCanvas()->getEventHandler()->broadcast(moveEvent);
		}
	  vtkInteractorStyleTrackballCamera::OnMouseMove();
	}
	
	virtual void OnLeftButtonUp () 
	{
	  int* mouse_coords = GetMousePosition();
	  if (debug) cout << "Released left mouse button." << mouse_coords[0] << " " << mouse_coords[1] << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_BUTTON_LEFT;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::RELEASED;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mouseReleasedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mouseReleasedEvent);
	  holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_NONE;
	  vtkInteractorStyleTrackballCamera::OnLeftButtonUp();

	}

	virtual void OnMiddleButtonUp () {
	
	  int* mouse_coords = GetMousePosition();
	  if (debug) cout << "Released middle mouse button." << mouse_coords[0] << " " << mouse_coords[1] << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_BUTTON_MIDDLE;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::RELEASED;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mouseReleasedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mouseReleasedEvent);
	  holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_NONE;
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
	}
	virtual void OnRightButtonUp () {
		
	  int* mouse_coords = GetMousePosition();
	  if (debug) cout << "Released right mouse button." << mouse_coords[0] << " " << mouse_coords[1] << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_BUTTON_RIGHT;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::RELEASED;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mouseReleasedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mouseReleasedEvent);
	  holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_NONE;
      vtkInteractorStyleTrackballCamera::OnRightButtonUp();
	}
 
    virtual void OnLeftButtonDown() 
    {
      int* mouse_coords = GetMousePosition();
	  if (debug) cout << "Pressed left mouse button." << mouse_coords[0] << " " << mouse_coords[1] << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_BUTTON_LEFT;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::PRESSED;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mousePressedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mousePressedEvent);
	  holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_LEFT;
      vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }
 
    virtual void OnMiddleButtonDown() 
    {
      int* mouse_coords = GetMousePosition();
	  if (debug) cout << "Pressed middle mouse button." << mouse_coords[0] << " " << mouse_coords[1] << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_BUTTON_MIDDLE;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::PRESSED;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mousePressedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mousePressedEvent);
	  holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_MIDDLE;
      vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    }
 
    virtual void OnRightButtonDown() 
    {
      int* mouse_coords = GetMousePosition();
	  if (debug) cout << "Pressed right mouse button." << mouse_coords[0] << " " << mouse_coords[1] << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_BUTTON_RIGHT;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::PRESSED;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mousePressedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mousePressedEvent);
	  holdButton_ = tgt::MouseEvent::MOUSE_BUTTON_RIGHT;
	  vtkInteractorStyleTrackballCamera::OnRightButtonDown();
    }

	virtual void OnMouseWheelForward ()
	{
	  int* mouse_coords = GetMousePosition();
	  if (debug) cout << "MouseWheel up" << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_WHEEL_UP;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::WHEEL;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mousePressedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mousePressedEvent);	 
	  vtkInteractorStyleTrackballCamera::OnMouseWheelForward();
	}

	virtual void OnMouseWheelBackward ()
	{
	  int* mouse_coords = GetMousePosition();
	  if (debug) cout << "MouseWheel down" << endl;
	  int win_height = voreensource->getCanvas()->getHeight();
	  int win_width = voreensource->getCanvas()->getWidth();
	  tgt::MouseEvent::MouseButtons pressedButton = tgt::MouseEvent::MOUSE_WHEEL_DOWN;
	  tgt::MouseEvent::MouseAction action = tgt::MouseEvent::WHEEL;
	  tgt::Event::Modifier tgtModifier = static_cast<tgt::Event::Modifier>(getModifier());
	  tgt::MouseEvent* mousePressedEvent = new tgt::MouseEvent(mouse_coords[0],mouse_coords[1],action,tgtModifier,pressedButton,tgt::ivec2(win_width, win_height));
	  voreensource->getCanvas()->getEventHandler()->broadcast(mousePressedEvent);	
	  vtkInteractorStyleTrackballCamera::OnMouseWheelBackward();
	  
	}

    virtual void OnKeyRelease()
    {
        vtkRenderWindowInteractor *rwi = this->Interactor;		
		char ch = rwi->GetKeyCode();  //KeyCodes for normal keys
		std::string special = rwi->GetKeySym(); //KeySym for special keys like F1-F2 etc
		int tgtModifier = getModifier();
		tgt::KeyEvent::KeyCode tgt_key = getKeyCode(ch);
		if (tgt_key == tgt::KeyEvent::K_UNKNOWN){
			tgt_key = getSpecialKeyCode(special);
		}
		if (debug) cout << "Pressed " << tgt_key << endl;
		tgt::KeyEvent* ke_release = new tgt::KeyEvent(tgt_key, tgtModifier, false);
		voreensource->getCanvas()->getEventHandler()->broadcast(ke_release);
        
        // Define custom key events for released keys here

        vtkInteractorStyleTrackballCamera::OnKeyRelease();

    }

	virtual void OnKeyPress()
	{
		vtkRenderWindowInteractor *rwi = this->Interactor;		
		
		char ch = rwi->GetKeyCode(); 
		std::string special = rwi->GetKeySym(); 
		int tgtModifier = getModifier();
		tgt::KeyEvent::KeyCode tgt_key = getKeyCode(ch);
		if (tgt_key == tgt::KeyEvent::K_UNKNOWN){
			tgt_key = getSpecialKeyCode(special);
		}
		if (debug) cout << "Pressed " << tgt_key << endl;
		tgt::KeyEvent* ke_press = new tgt::KeyEvent(tgt_key, tgtModifier, true);
		voreensource->getCanvas()->getEventHandler()->broadcast(ke_press);
		

		// Define custom key events for pressed keys here
        std::stringstream help;

		switch (ch) {
        case 'h':
            help << "Commands:\n";
            help << " c     : Toggle Clipping Plane\n";
            help << " + / -         : Move Clipping Plane\n";
            help << " Shift+/Shift- : Change Clipping Stepping\n";
            help << " k     : Select transfer function key\n";
            help << " arrows: Move selected tf-key\n";
            help << " t     : Reset Transfer Function\n";
            help << " z     : Clear Transfer Function\n";
            help << " r     : Toggle Voreen-Renderers\n";
            help << " 1     : Toggle Compositing Mode\n";
            help << " 2     : Toggle Shading\n";
            
            help << " 6     : Load TF\n";
            help << " PGUP/DOWN : Increase/Decrease Sampling Rate\n";
            help << " 7-0   : Add Keys to Transfer Function (7+8 grey|9+0 color)\n";
            help << " p     : List Processors\n";
            help << " b     : Create a transfer function (example)\n";
            help << "F1-F12 : List Processor Properties";
            voreensource->PrintMessage(help.str());
            break;
		case '+':
			 if (rwi->GetShiftKey())
                voreensource->setClippingStepping(true);
            else
                voreensource->moveClippingPlane(true);
			break;
		case '-':
            if (rwi->GetShiftKey())
                voreensource->setClippingStepping(false);
            else
                voreensource->moveClippingPlane(false);
			break;
		case 'p':           
			voreensource->listProcessors();
			break;
        case 'a':
            voreensource->loadTransferFunction("all", "../../data/transferfuncs/walnut.tfi");
            break;
        case 'b':
            //TF for bonsai.raw-Example
            voreensource->clearTransferFunction("all");
            voreensource->addTransferFunctionKey("all", 0,0,0,0,0);

            voreensource->addTransferFunctionKey("all", 0.35,0.5,0.25,0.1,0);
            voreensource->addTransferFunctionKey("all", 0.36,0.5,0.25,0.1,1);
            voreensource->addTransferFunctionKey("all", 0.58,0.5,0.25,0.1,1);
            voreensource->addTransferFunctionKey("all", 0.59,0.5,0.25,0.1,0);

            voreensource->addTransferFunctionKey("all", 0.75,0,0,1,0);
            voreensource->addTransferFunctionKey("all", 0.76,0,0,1,1);
            voreensource->addTransferFunctionKey("all", 1,0,0,1,1);
            break;
        case '5': 
            voreensource->enableBoundingBox();
        case 'k':
            voreensource->selectTransferfunctionKey(voreensource->getSelectedKey() + 1);
             break;
        case 's':
            voreensource->statusToggle();
            break;
		case 't':
			voreensource->resetTransferFunction("all");
			break;
        case 'z':
			voreensource->clearTransferFunction("SingleVolumeRaycaster");
			break;
        case '1':
            if (voreensource->getCompositeStyle() == "dvr")
                voreensource->setCompositingMIP();
            else
            if (voreensource->getCompositeStyle() == "mip")
                voreensource->setCompositingISO();
            else
            if (voreensource->getCompositeStyle() == "iso")
                voreensource->setCompositingFHP();
            else
            if (voreensource->getCompositeStyle() == "fhp")
                voreensource->setCompositingFHN();            
            else
            if (voreensource->getCompositeStyle() == "fhn")
                voreensource->setCompositingDVR();
            break;
        case '2':
            voreensource->setShading((voreensource->getShading()+1)%7);
            break;
        case 'c':
            if (voreensource->getSelectedClippingPlane() == "left")
                voreensource->selectClippingPlane("right");
            else
                if (voreensource->getSelectedClippingPlane() == "right")
                voreensource->selectClippingPlane("top");
            else
            if (voreensource->getSelectedClippingPlane() == "top")
                voreensource->selectClippingPlane("bottom");
            else
                if (voreensource->getSelectedClippingPlane() == "bottom")
                voreensource->selectClippingPlane("front");
            else
            if (voreensource->getSelectedClippingPlane() == "front")
                voreensource->selectClippingPlane("back");
            else
            if (voreensource->getSelectedClippingPlane() == "back")
                voreensource->selectClippingPlane("left");
            else
                voreensource->selectClippingPlane("left");
            break;
        case '6': 
                voreensource->loadTransferFunction("all", "../../data/transferfuncs/walnut.tfi");
            break;
        case '7':
             voreensource->addTransferFunctionKey(voreensource->getRaycasterName(),0,0,0,0,0);
            break;
        case '8':
             voreensource->addTransferFunctionKey(voreensource->getRaycasterName(),0.05,1,1,1,1);
            break;
        case '9':
             voreensource->addTransferFunctionKey(voreensource->getRaycasterName(),0.05,1,0,0,0.2);
            break;
        case '0':
             voreensource->addTransferFunctionKey(voreensource->getRaycasterName(),0.1,0.5,1,0,0.4);
            break;
		case 'r':
            if (voreensource->isSingleVolumeRaycaster())
                voreensource->setRendererToSimpleRayCaster();
            else
                voreensource->setRendererToSingleVolumeRayCaster();
			break;
        }

        if (special.compare("Left") == 0){
            voreensource->moveTransferFunctionKey(0);
        }

        if (special.compare("Up") == 0){
            voreensource->moveTransferFunctionKey(1);
        }

        if (special.compare("Right") == 0){
            voreensource->moveTransferFunctionKey(2);   
        }

        if (special.compare("Down") == 0){
            voreensource->moveTransferFunctionKey(3);
        }

        if (special.compare("Next") == 0){voreensource->decreaseSamplingRate();}
        if (special.compare("Prior") == 0){voreensource->increaseSamplingRate();}

        if (special.compare("F1") == 0) voreensource->listProperties(0);
		if (special.compare("F2") == 0) voreensource->listProperties(1);
		if (special.compare("F3") == 0) voreensource->listProperties(2);
		if (special.compare("F4") == 0) voreensource->listProperties(3);
		if (special.compare("F5") == 0) voreensource->listProperties(4);
		if (special.compare("F6") == 0) voreensource->listProperties(5);
        if (special.compare("F7") == 0) voreensource->listProperties(6);
        if (special.compare("F8") == 0) voreensource->listProperties(7);
        if (special.compare("F9") == 0) voreensource->listProperties(8);
        if (special.compare("F10") == 0) voreensource->listProperties(9);
        if (special.compare("F11") == 0) voreensource->listProperties(10);
        if (special.compare("F12") == 0) voreensource->listProperties(11);
		
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}
 
};

vtkCxxRevisionMacro(vtkVoreenInteractor, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkVoreenInteractor);
