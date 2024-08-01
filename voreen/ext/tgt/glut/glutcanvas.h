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

#ifndef TGT_GLUTCANVAS_H
#define TGT_GLUTCANVAS_H

#include "tgt/glcanvas.h"
#include "tgt/event/keyevent.h"
#include "tgt/glut/glutmouse.h"

#include <GL/glut.h>

#include <stack>

namespace tgt {

class GLUTCanvas : public GLCanvas {
public:
    /// Used to hold the canvases for glut callback functions.
    /// Has to be static due to glut restrictions
    static GLUTCanvas* Canvases_[];


    /// Constructor
    GLUTCanvas(const std::string& title = "",
               const ivec2& size = ivec2(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT),
               const GLCanvas::Buffers buffers = RGBADD );

    ~GLUTCanvas();

public:
    /// create a glut window and stuff for this canvas
    void init();

    /// swap buffers
    void swap();

    /// Toggles between fullscreen an no fullsreen.
    /// Saves window-states (size, position) before switching to fullscreen to be able
    /// to switch back to them later
    void toggleFullScreen();

    /// if canvas is displayed in fullscreen mode
    bool isFullScreen();

    /// to be called by application to cause immidiate (re)painting on the canvas
    void repaint();

    /// to be called by application to cause (re)painting on the canvas when idle next time
    void update();

protected:
    /** @override GLContextBase */
    virtual void activate();
    virtual bool isActive();

public:
    /// getter
    GLUTMouse* getMouse();
    /// setter
    void setMouse(GLUTMouse* glutMouse);

    /// returns the number glut assigned to this canvas' window
    inline int getWindowID();


    /**
     * These are the functions used to translate GLUT events into TGT events.
     * They are by far not complete.
     *
     * See the GLUT API for detailed explanation:
     * http://www.opengl.org/resources/libraries/glut/spec3/spec3.html
     */

    /**
     * Tells the painter to repaint.
     *
     * A general note on callback functions:
     *
     * They should be called from matching GLUT callback functions in your program.
     *
     * Planned is to register some static (they have to be for GLUT) dummy functions
     * as Callback functions with GLUT which then call the following "translation"
     * functions.
     *
     * These "translation" functions handle the actual translations from GLUT events
     * into TGT events. Some of them need an additional GLUT modifier which has to
     * be fetched from GLUT in the according dummy functions.
     *
     */
    void display();

    /**
     * Mouse moved with button pressed.
     * @see display()
     */
    void mouseMotion(const int& x, const int& y);

    /**
     * Mouse moved without button pressed.
     * @see display()
     */
    void passiveMouseMotion(const int& x, const int& y);

    /**
     * Handles button presses and releases.
     * @see display()
     * @param modifier Is the GLUT Modifier Bitmap returned by glutGetModifiers()
     */
    void mousePressed(const int& button, const int& state, const int& x, const int& y, const int& modifier );

    /**
     * Handles "normal" keys that are passed by GLUT in plain Ascii.
     * Note: It is very likely that your OS/Window Manager will not pass some special keys and/or modifiers.
     *
     * @see display()
     * @param modifier Is the GLUT Modifier Bitmap returned by glutGetModifiers()
     */
    void keyboard(const unsigned char& key, const int& x, const int& y, const int& modifier);

    /**
     * Handles Special Keys UP, DOWN, LEFT, RIGHT, INSERT, HOME, END, PGDOWN, PGUP, F1, ..., F12.
     * It is very likely that your OS/Window Manager will not pass every key and modifier.
     * According to GLUT API backspace, delete and escape characters are ASCII.
     *
     * Note: Actually generates two (maybe four) KeyEvents for pressing and releasing
     * (maybe shifting). This is an inconvinience of GLUT. Which only announces presses.
     *
     * @see display()
     * @param modifier Is the GLUT Modifier Bitmap returned by glutGetModifiers()
     */
    void keyboardSpecial( const int& key, const int& x, const int& y, const int& modifier );

    /**
     * Tells TGT that window size has changed.
     * @see display()
     */
    void reshape( int width, int height );

    /**
     *
     * @see display()
     */
    void visibility( const int& visible );

    //    void idle(); deprecated


//------------------------------------------------------------------------------
// callback wrapper functions
//
// remember: these are wrapping the member functions of GLUTCanvas.


    /*
     * These are the functions used to translate GLUT events into TGT events.
     *
     * See the GLUT API for detailed explanation:
     * http://www.opengl.org/resources/libraries/glut/spec3/spec3.html
     */

    static void displayFunc();
    static void mouseMotionFunc( int x, int y );
    static void passiveMouseMotionFunc( int x, int y);
    static void mousePressedFunc( int button, int state, int x, int y );
    static void keyboardFunc( unsigned char key, int x,  int y );
    static void keyboardSpecialFunc( int key, int x,  int y );
    static void reshapeFunc( int width,  int height );
    static void idleFunc();
    static void visibilityFunc(int visible);

private:
    
    // glut only gives position by callback functions
    // so we have to set them from here if wanted so
    GLUTMouse* glutMouse_;

    // GLUT state
    int windowID_;     ///< store the id of our OpenGL window
    ivec2 windowPosition_;
    ivec2 noFullScreenSize_;
    bool fullscreen_;
    int millisec_; /// timer intervall
    bool timer_; /// timer in use
    MouseEvent::MouseButtons holdButton_;   ///< needed to bypass glut restriction

    /// kepp optional GLUTMouse up to date
    void keepMouseUpdated( const int& x, const int& y);


/// helper functions ///////////////////////////////////////////////////////////////////

    /**
     * Helper to map from GLUT modifier bitmap to TGT modifier bitmap.
     */
    KeyEvent::KeyCode getKeyCode(const int& key);

    /**
     * Helper to convert plain ascii (from GLUT) to TGT keys.
     */
    KeyEvent::KeyCode getSpecialKeyCode(const int& key);

    /**
     * Helper to convert "special" Keys from GLUT to TGT keys.
     */
    int getModifier(const int& glutModifier);

    /// make a glut display mode bitmask out of the buffer specification in buffers_
    inline unsigned int getDisplayMode();

    void registerCallbacks();
};

}
#endif // TGT_GLUTCANVAS_H
