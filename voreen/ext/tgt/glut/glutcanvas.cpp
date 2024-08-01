/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/glut/glutcanvas.h"

#include "tgt/glcontextmanager.h"

namespace tgt {

/// used to hold the canvases registered in our app.
/// must be global for the glut wrapper.
/// see GLUTApplication for more details
enum { MAX_NUMBER_OF_WINDOWS = 256 };

GLUTCanvas * GLUTCanvas::Canvases_[MAX_NUMBER_OF_WINDOWS];

GLUTCanvas::GLUTCanvas(const std::string& title,
                       const ivec2& size,
                       const GLCanvas::Buffers buffers)
    : GLCanvas(title, size, buffers)
{
    fullscreen_ = false;
    glutMouse_ = 0;
    holdButton_ = MouseEvent::MOUSE_BUTTON_NONE;
}

GLUTCanvas::~GLUTCanvas() {
    Canvases_[windowID_] = 0;
    glutDestroyWindow(windowID_);
}

void GLUTCanvas::init() {
    glutInitDisplayMode(getDisplayMode());

    glutInitWindowSize(size_.x, size_.y);
    glutCreateWindow(title_.c_str());

    windowID_ = glutGetWindow();
    Canvases_[windowID_] = this;
    registerCallbacks();

    rgbaSize_ = ivec4(glutGet(GLUT_WINDOW_RED_SIZE),
                      glutGet(GLUT_WINDOW_GREEN_SIZE),
                      glutGet(GLUT_WINDOW_BLUE_SIZE),
                      glutGet(GLUT_WINDOW_ALPHA_SIZE)) ;
    stencilSize_ = glutGet(GLUT_WINDOW_STENCIL_SIZE);
    depthSize_ = glutGet(GLUT_WINDOW_DEPTH_SIZE);
    doubleBuffered_ = (0 != glutGet(GLUT_WINDOW_DOUBLEBUFFER));
    quadBuffered_ =  (0 != glutGet(GLUT_WINDOW_STEREO));

    GLCanvas::init();
}

void GLUTCanvas::swap() {
    glutSwapBuffers();
}

void GLUTCanvas::repaint() {
    GLContextStateGuard guard(this);
    display();
}

void GLUTCanvas::update() {
    GLContextStateGuard guard(this);
    glutPostRedisplay();
}

void GLUTCanvas::toggleFullScreen() {
    if (!fullscreen_) {
        noFullScreenSize_ = size_;
        windowPosition_ = ivec2(glutGet(GLUT_WINDOW_X), glutGet(GLUT_WINDOW_Y));
        glutFullScreen();
    } else {
        // reshaping calls sizeChanged by reshape
        // this updates size_
        glutReshapeWindow(noFullScreenSize_.x, noFullScreenSize_.y);
        glutPositionWindow(windowPosition_.x, windowPosition_.y);
    }
    fullscreen_ = !fullscreen_;
}

bool GLUTCanvas::isFullScreen() {
    return fullscreen_;
}

GLUTMouse* GLUTCanvas::getMouse() {
    return glutMouse_;
}

void GLUTCanvas::setMouse(GLUTMouse* glutMouse) {
    glutMouse_ = glutMouse;
}

int GLUTCanvas::getWindowID() {
    return windowID_;
}

void GLUTCanvas::activate() {
    glutSetWindow(windowID_);
}

bool GLUTCanvas::isActive() {
    return glutGetWindow() == windowID_;
}

//------------------------------------------------------------------------------
// Callback functions for GLUT.

void GLUTCanvas::display() {
    paint();
}

void GLUTCanvas::mouseMotion(const int& x, const int& y) {
    // maybe there is a glutMouse that needs update
    keepMouseUpdated(x, y);

    // Create tgt and broadcast it
    MouseEvent* moveEvent = new MouseEvent(x, y, MouseEvent::MOTION, Event::MODIFIER_NONE, holdButton_, getSize());
    getEventHandler()->broadcast(moveEvent);
    paint();
}

void GLUTCanvas::passiveMouseMotion(const int& x, const int& y) {
    // maybe there is a glutMouse that needs update
    keepMouseUpdated(x, y);

    // Create tgt and broadcast it
    MouseEvent* moveEvent = new MouseEvent(x, y, MouseEvent::MOTION, Event::MODIFIER_NONE, MouseEvent::MOUSE_BUTTON_NONE, getSize());
    getEventHandler()->broadcast(moveEvent);
    paint();
}

void GLUTCanvas::mousePressed(const int& button, const int& state, const int& x, const int& y, const int& modifier) {
    // maybe there is a glutMouse that needs update
    keepMouseUpdated(x, y);

    // map glut event to tgt event
    MouseEvent::MouseButtons pressedButton;
    if (button ==  GLUT_LEFT_BUTTON)
        pressedButton = MouseEvent::MOUSE_BUTTON_LEFT;
    else if (button == GLUT_RIGHT_BUTTON)
        pressedButton = MouseEvent::MOUSE_BUTTON_RIGHT;
    else
        pressedButton = MouseEvent::MOUSE_BUTTON_MIDDLE;

    // Create tgt action and update holdButton.
    MouseEvent::MouseAction action;
    if (state == GLUT_DOWN) {
        action = MouseEvent::PRESSED;
        holdButton_ = pressedButton;
    }
    else{
        action = MouseEvent::RELEASED;
        holdButton_ = MouseEvent::MOUSE_BUTTON_NONE;
    }

    // calculate right TGT modifier
    Event::Modifier tgtModifier = static_cast<Event::Modifier>(getModifier(modifier));

    // Create and broadcast event
    MouseEvent* mousePressedEvent = new MouseEvent(x, y, action, tgtModifier, pressedButton, getSize());
    getEventHandler()->broadcast(mousePressedEvent);
    paint();
}

void GLUTCanvas::keyboard(const unsigned char& key, const int& x, const int& y, const int& modifier)
{
    // glutKeyrepeat
    // glutKeyUpFunc()

    // maybe there is a GLUTMouse that needs update
    keepMouseUpdated(x, y);

    // get TGT keyCode and modifier
    KeyEvent::KeyCode tgtKey = getKeyCode(key);
    int tgtModifier = getModifier(modifier);

    // TODO: encapsulate this by shift press and release events?
    tgt::KeyEvent* ke_press = new tgt::KeyEvent(tgtKey, tgtModifier, false, true);
    tgt::KeyEvent* ke_release = new tgt::KeyEvent(tgtKey, tgtModifier, false, false);
    eventHandler_->broadcast(ke_press);
    eventHandler_->broadcast(ke_release);
    paint();
}

void GLUTCanvas::keyboardSpecial(const int& key, const int& x, const int& y, const int& modifier) {
    // maybe there is a glutMouse that needs update
    keepMouseUpdated(x, y);

    // get TGT keyCode and modifier
    KeyEvent::KeyCode tgtKey = getSpecialKeyCode(key);
    int tgtModifier = getModifier(modifier);

    // TODO: encapsulate this by shift press and release events?
    tgt::KeyEvent* ke_press = new tgt::KeyEvent(tgtKey, tgtModifier, false, true);
    tgt::KeyEvent* ke_release = new tgt::KeyEvent(tgtKey, tgtModifier, false, false);
    eventHandler_->broadcast(ke_press);
    eventHandler_->broadcast(ke_release);
    paint();
}

void GLUTCanvas::reshape(int width, int height) {
    sizeChanged(ivec2(width, height));
    paint();
}

void GLUTCanvas::visibility(const int& /*visible*/) {
}

//------------------------------------------------------------------------------
// Helpers


/**
 * Helper to keep track of the mouse. Mainly used to have an optional GLUTMouse updated.
 */
void GLUTCanvas::keepMouseUpdated(const int& x, const int& y) {
    if (glutMouse_)
        glutMouse_->setPosition(ivec2(x, y));
}

/**
 * Helper to map from GLUT modifier bitmap to TGT modifier bitmap.
 */
int GLUTCanvas::getModifier(const int& glutModifier) {
    int result = 0;

    if (glutModifier & GLUT_ACTIVE_SHIFT)
        result |= tgt::Event::SHIFT;

    if (glutModifier & GLUT_ACTIVE_ALT)
        result |= tgt::Event::ALT;

    if (glutModifier & GLUT_ACTIVE_CTRL)
        result |= tgt::Event::CTRL;

    return result;
}



/**
 * Helper to convert plain ascii (from GLUT) to TGT keys.
 */
KeyEvent::KeyCode GLUTCanvas::getKeyCode(const int& key) {
    // Attention!!! This resides on tgt internal (ascii at the moment) settings which could change.
    // Normaly it should be handled in a switch statement.
    // this is faster :)

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



/**
 * Helper to convert "special" Keys from GLUT to TGT keys.
 */
KeyEvent::KeyCode GLUTCanvas::getSpecialKeyCode(const int& key) {
    // Note: According to its API, this is everything GLUT has to offer.
    switch (key) {
    case GLUT_KEY_F1:           return tgt::KeyEvent::K_F1;
    case GLUT_KEY_F2:           return tgt::KeyEvent::K_F2;
    case GLUT_KEY_F3:           return tgt::KeyEvent::K_F3;
    case GLUT_KEY_F4:           return tgt::KeyEvent::K_F4;
    case GLUT_KEY_F5:           return tgt::KeyEvent::K_F5;
    case GLUT_KEY_F6:           return tgt::KeyEvent::K_F6;
    case GLUT_KEY_F7:           return tgt::KeyEvent::K_F7;
    case GLUT_KEY_F8:           return tgt::KeyEvent::K_F8;
    case GLUT_KEY_F9:           return tgt::KeyEvent::K_F9;
    case GLUT_KEY_F10:          return tgt::KeyEvent::K_F10;
    case GLUT_KEY_F11:          return tgt::KeyEvent::K_F11;
    case GLUT_KEY_F12:          return tgt::KeyEvent::K_F12;

    case GLUT_KEY_LEFT:         return tgt::KeyEvent::K_LEFT;
    case GLUT_KEY_UP:           return tgt::KeyEvent::K_UP;
    case GLUT_KEY_RIGHT:        return tgt::KeyEvent::K_RIGHT;
    case GLUT_KEY_DOWN:         return tgt::KeyEvent::K_DOWN;
    case GLUT_KEY_PAGE_UP:      return tgt::KeyEvent::K_PAGEUP;
    case GLUT_KEY_PAGE_DOWN:    return tgt::KeyEvent::K_PAGEDOWN;
    case GLUT_KEY_HOME:         return tgt::KeyEvent::K_HOME;
    case GLUT_KEY_END:          return tgt::KeyEvent::K_END;
    case GLUT_KEY_INSERT:       return tgt::KeyEvent::K_INSERT;

    default:                    return tgt::KeyEvent::K_UNKNOWN;
    }
}


/**
 * Helper to convert buffer specitication from tgt buffer enum to
 * glut display mode.
 */
unsigned int GLUTCanvas::getDisplayMode() {
  return ((RGB_BUFFER     & buffers_ ? GLUT_RGB         : 0) |
          (RGBA_BUFFER    & buffers_ ? GLUT_RGBA        : 0) |
          (ALPHA_BUFFER   & buffers_ ? GLUT_ALPHA       : 0) |
          (DEPTH_BUFFER   & buffers_ ? GLUT_DEPTH       : 0) |
          (DOUBLE_BUFFER  & buffers_ ? GLUT_DOUBLE      : 0) |
          (STENCIL_BUFFER & buffers_ ? GLUT_STENCIL     : 0) |
          (ACCUM_BUFFER   & buffers_ ? GLUT_ACCUM       : 0) |
          (QUAD_BUFFER & buffers_ ? GLUT_STEREO      : 0) |
          (MULTISAMPLING  & buffers_ ? GLUT_MULTISAMPLE : 0));
}

//------------------------------------------------------------------------------
// callback wrapper functions
//
// remember: these are wrapping the member functions of GLUTCanvas.

void GLUTCanvas::registerCallbacks() {
    glutDisplayFunc(displayFunc);
    glutReshapeFunc(reshapeFunc);
    glutKeyboardFunc(keyboardFunc);
    glutSpecialFunc(keyboardSpecialFunc);
    glutMouseFunc(mousePressedFunc);
    glutMotionFunc(mouseMotionFunc);
    glutPassiveMotionFunc(mouseMotionFunc);
    glutVisibilityFunc(visibilityFunc);
}

void GLUTCanvas::displayFunc(void) {
    Canvases_[glutGetWindow()]->display();
}

void GLUTCanvas::keyboardFunc(unsigned char key, int x, int y) {
    Canvases_[glutGetWindow()]->keyboard(key, x, y, glutGetModifiers());
}

void GLUTCanvas::mouseMotionFunc(int x, int y) {
    Canvases_[glutGetWindow()]->mouseMotion(x, y);
}

void GLUTCanvas::mousePressedFunc(int button, int state, int x, int y) {
    Canvases_[glutGetWindow()]->mousePressed(button, state, x, y, glutGetModifiers());
}

void GLUTCanvas::passiveMouseMotionFunc(int x, int y) {
    Canvases_[glutGetWindow()]->passiveMouseMotion(x, y);
}

void GLUTCanvas::reshapeFunc(int width, int height) {
    Canvases_[glutGetWindow()]->reshape(width, height);
}

void GLUTCanvas::keyboardSpecialFunc(int key, int x, int y) {
    Canvases_[glutGetWindow()]->keyboardSpecial(key, x, y, glutGetModifiers());
}

void GLUTCanvas::visibilityFunc(int visible) {
    Canvases_[glutGetWindow()]->visibility(visible);
}

} //namespace tgt
