/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/qt/widgets/keydetectorwidget.h"

namespace voreen {

ModifierDetectorWidget::ModifierDetectorWidget(QWidget* parent)
    : QLineEdit(parent)
{
    setReadOnly(true);
    setAlignment(Qt::AlignHCenter);
}

void ModifierDetectorWidget::keyPressEvent(QKeyEvent* event) {
    event->accept();

    setText(getStringForModifier(event->modifiers()));

    emit modifierChanged(event->modifiers());
}

QString ModifierDetectorWidget::getStringForModifier(Qt::KeyboardModifiers modifier) {
    QString modifierString = "";

    if (modifier & Qt::ShiftModifier)
        modifierString += "Shift";

    if (modifier & Qt::ControlModifier)
        modifierString += ((modifierString == "") ? "Ctrl" : " + Ctrl");

    if (modifier & Qt::AltModifier)
        modifierString += ((modifierString == "") ? "Alt" : " + Alt");

    if (modifier & Qt::MetaModifier)
        modifierString += ((modifierString == "") ? "Meta" : " + Meta");

    return modifierString;
}

tgt::Event::Modifier ModifierDetectorWidget::getTGTModifierFromQt(Qt::KeyboardModifiers qtKey) {
    int result = 0;
    if (qtKey & Qt::ShiftModifier)
        result |= tgt::Event::SHIFT;
    if (qtKey & Qt::ControlModifier)
        result |= tgt::Event::CTRL;
    if (qtKey & Qt::MetaModifier)
        result |= tgt::Event::META;
    if (qtKey & Qt::AltModifier)
        result |= tgt::Event::ALT;

    return tgt::Event::Modifier(result);
}

Qt::KeyboardModifiers ModifierDetectorWidget::getQtModifierFromTGT(tgt::Event::Modifier tgtKey) {
    int result = 0;

    if (tgtKey & tgt::Event::SHIFT)
        result |=  Qt::ShiftModifier;
    if (tgtKey & tgt::Event::CTRL)
        result |=  Qt::ControlModifier;
    if (tgtKey & tgt::Event::ALT)
        result |=  Qt::AltModifier;
    if (tgtKey & tgt::Event::META)
        result |=  Qt::MetaModifier;

    return Qt::KeyboardModifier(result);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------

KeyDetectorWidget::KeyDetectorWidget(bool ignoreModifierKeys, QWidget* parent)
    : QLineEdit(parent)
    , ignoreModifierKeys_(ignoreModifierKeys)
{
    setReadOnly(true);
    setAlignment(Qt::AlignHCenter);
}

void KeyDetectorWidget::keyPressEvent(QKeyEvent* event) {
    event->accept();

    setText(getStringForKey(event->key(), ignoreModifierKeys_));

    emit keyChanged(event->key());
}

QString KeyDetectorWidget::getStringForKey(int key, bool ignoreModifierKeys) {
    switch (key) {
    case Qt::Key_Escape: return "Esc";
    case Qt::Key_Tab: return "Tab";
    case Qt::Key_Backspace: return "Backspace";
    case Qt::Key_Return: return "Return";
    case Qt::Key_Enter: return "Enter";
    case Qt::Key_Insert: return "Insert";
    case Qt::Key_Delete: return "Delete";
    case Qt::Key_Pause: return "Pause";
    case Qt::Key_Print: return "Print";
    case Qt::Key_Clear: return "Clear";
    case Qt::Key_Home: return "Home";
    case Qt::Key_End: return "End";
    case Qt::Key_Left: return "Left Arrow";
    case Qt::Key_Up: return "Up Arrow";
    case Qt::Key_Right: return "Right Arrow";
    case Qt::Key_Down: return "Down Arrow";
    case Qt::Key_PageUp: return "Page Up";
    case Qt::Key_PageDown: return "Page Down";
    case Qt::Key_ScrollLock: return "ScrollLock";
    case Qt::Key_F1: return "F1";
    case Qt::Key_F2: return "F2";
    case Qt::Key_F3: return "F3";
    case Qt::Key_F4: return "F4";
    case Qt::Key_F5: return "F5";
    case Qt::Key_F6: return "F6";
    case Qt::Key_F7: return "F7";
    case Qt::Key_F8: return "F8";
    case Qt::Key_F9: return "F9";
    case Qt::Key_F10: return "F10";
    case Qt::Key_F11: return "F11";
    case Qt::Key_F12: return "F12";
    case Qt::Key_Space: return "Space";
    case Qt::Key_Exclam: return "!";
    case Qt::Key_Dollar: return "$";
    case Qt::Key_Percent: return "%";
    case Qt::Key_Ampersand: return "&";
    case Qt::Key_Apostrophe: return "Key_Apostrophe";
    case Qt::Key_ParenLeft: return "(";
    case Qt::Key_ParenRight: return ")";
    case Qt::Key_Asterisk: return "*";
    case Qt::Key_Plus: return "+";
    case Qt::Key_Comma: return ",";
    case Qt::Key_Minus: return "-";
    case Qt::Key_Period: return ".";
    case Qt::Key_Slash: return "/";
    case Qt::Key_0: return "0";
    case Qt::Key_1: return "1";
    case Qt::Key_2: return "2";
    case Qt::Key_3: return "3";
    case Qt::Key_4: return "4";
    case Qt::Key_5: return "5";
    case Qt::Key_6: return "6";
    case Qt::Key_7: return "7";
    case Qt::Key_8: return "8";
    case Qt::Key_9: return "9";
    case Qt::Key_Colon: return ":";
    case Qt::Key_Semicolon: return ";";
    case Qt::Key_Less: return "<";
    case Qt::Key_Equal: return "=";
    case Qt::Key_Greater: return ">";
    case Qt::Key_Question: return "?";
    case Qt::Key_At: return "@";
    case Qt::Key_A: return "A";
    case Qt::Key_B: return "B";
    case Qt::Key_C: return "C";
    case Qt::Key_D: return "D";
    case Qt::Key_E: return "E";
    case Qt::Key_F: return "F";
    case Qt::Key_G: return "G";
    case Qt::Key_H: return "H";
    case Qt::Key_I: return "I";
    case Qt::Key_J: return "J";
    case Qt::Key_K: return "K";
    case Qt::Key_L: return "L";
    case Qt::Key_M: return "M";
    case Qt::Key_N: return "N";
    case Qt::Key_O: return "O";
    case Qt::Key_P: return "P";
    case Qt::Key_Q: return "Q";
    case Qt::Key_R: return "R";
    case Qt::Key_S: return "S";
    case Qt::Key_T: return "T";
    case Qt::Key_U: return "U";
    case Qt::Key_V: return "V";
    case Qt::Key_W: return "W";
    case Qt::Key_X: return "X";
    case Qt::Key_Y: return "Y";
    case Qt::Key_Z: return "Z";
    case Qt::Key_BracketLeft: return "Key_BracketLeft";
    case Qt::Key_Backslash: return "\\";
    case Qt::Key_BracketRight: return "Key_BracketRight";
    case Qt::Key_Underscore: return "Key_Underscore";
    case Qt::Key_Shift:
        if (!ignoreModifierKeys)
            return "Shift";
        else
            return "";
        break;
    case Qt::Key_Control:
        if (!ignoreModifierKeys)
            return "Ctrl";
        else
            return "";
        break;
    case Qt::Key_Meta:
        if (!ignoreModifierKeys)
            return "Meta";
        else
            return "";
        break;
    case Qt::Key_Alt:
        if (!ignoreModifierKeys)
            return "Alt";
        else
            return "";
        break;
    case Qt::Key_AltGr:
        if (!ignoreModifierKeys)
            return "AltGr";
        else
            return "";
        break;
    default:
        return "unknown key";
}
}

tgt::KeyEvent::KeyCode KeyDetectorWidget::getTGTKeyFromQt(int qtKey) {
    switch (qtKey) {
    case Qt::Key_Escape: return tgt::KeyEvent::K_ESCAPE;
    case Qt::Key_Tab: return tgt::KeyEvent::K_TAB;
    case Qt::Key_Backspace: return tgt::KeyEvent::K_BACKSPACE;
    case Qt::Key_Return: return tgt::KeyEvent::K_RETURN;
    case Qt::Key_Enter: return tgt::KeyEvent::K_KP_ENTER;
    case Qt::Key_Insert: return tgt::KeyEvent::K_INSERT;
    case Qt::Key_Delete: return tgt::KeyEvent::K_DELETE;
    case Qt::Key_Pause: return tgt::KeyEvent::K_PAUSE;
    case Qt::Key_Print: return tgt::KeyEvent::K_PRINT;
    case Qt::Key_Clear: return tgt::KeyEvent::K_CLEAR;
    case Qt::Key_Home: return tgt::KeyEvent::K_HOME;
    case Qt::Key_End: return tgt::KeyEvent::K_END;
    case Qt::Key_Left: return tgt::KeyEvent::K_LEFT;
    case Qt::Key_Up: return tgt::KeyEvent::K_UP;
    case Qt::Key_Right: return tgt::KeyEvent::K_RIGHT;
    case Qt::Key_Down: return tgt::KeyEvent::K_DOWN;
    case Qt::Key_PageUp: return tgt::KeyEvent::K_PAGEUP;
    case Qt::Key_PageDown: return tgt::KeyEvent::K_PAGEDOWN;
    case Qt::Key_Shift: return tgt::KeyEvent::K_LSHIFT;
    case Qt::Key_Control: return tgt::KeyEvent::K_LCTRL;
    case Qt::Key_Meta: return tgt::KeyEvent::K_RMETA;
    case Qt::Key_Alt: return tgt::KeyEvent::K_LALT;
    case Qt::Key_AltGr: return tgt::KeyEvent::K_MODE;
    case Qt::Key_CapsLock: return tgt::KeyEvent::K_CAPSLOCK;
    case Qt::Key_NumLock: return tgt::KeyEvent::K_NUMLOCK;
    case Qt::Key_ScrollLock: return tgt::KeyEvent::K_SCROLLOCK;
    case Qt::Key_F1: return tgt::KeyEvent::K_F1;
    case Qt::Key_F2: return tgt::KeyEvent::K_F2;
    case Qt::Key_F3: return tgt::KeyEvent::K_F3;
    case Qt::Key_F4: return tgt::KeyEvent::K_F4;
    case Qt::Key_F5: return tgt::KeyEvent::K_F5;
    case Qt::Key_F6: return tgt::KeyEvent::K_F6;
    case Qt::Key_F7: return tgt::KeyEvent::K_F7;
    case Qt::Key_F8: return tgt::KeyEvent::K_F8;
    case Qt::Key_F9: return tgt::KeyEvent::K_F9;
    case Qt::Key_F10: return tgt::KeyEvent::K_F10;
    case Qt::Key_F11: return tgt::KeyEvent::K_F11;
    case Qt::Key_F12: return tgt::KeyEvent::K_F12;
    case Qt::Key_Space: return tgt::KeyEvent::K_SPACE;
    case Qt::Key_Exclam: return tgt::KeyEvent::K_EXCLAIM;
    case Qt::Key_QuoteDbl: return tgt::KeyEvent::K_QUOTEDBL;
    case Qt::Key_NumberSign: return tgt::KeyEvent::K_HASH;
    case Qt::Key_Dollar: return tgt::KeyEvent::K_DOLLAR;
    case Qt::Key_Ampersand: return tgt::KeyEvent::K_AMPERSAND;
    case Qt::Key_ParenLeft: return tgt::KeyEvent::K_LEFTPAREN;
    case Qt::Key_ParenRight: return tgt::KeyEvent::K_RIGHTPAREN;
    case Qt::Key_Asterisk: return tgt::KeyEvent::K_ASTERISK;
    case Qt::Key_Plus: return tgt::KeyEvent::K_PLUS;
    case Qt::Key_Comma: return tgt::KeyEvent::K_COMMA;
    case Qt::Key_Minus: return tgt::KeyEvent::K_MINUS;
    case Qt::Key_Period: return tgt::KeyEvent::K_PERIOD;
    case Qt::Key_Slash: return tgt::KeyEvent::K_SLASH;
    case Qt::Key_0: return tgt::KeyEvent::K_0;
    case Qt::Key_1: return tgt::KeyEvent::K_1;
    case Qt::Key_2: return tgt::KeyEvent::K_2;
    case Qt::Key_3: return tgt::KeyEvent::K_3;
    case Qt::Key_4: return tgt::KeyEvent::K_4;
    case Qt::Key_5: return tgt::KeyEvent::K_5;
    case Qt::Key_6: return tgt::KeyEvent::K_6;
    case Qt::Key_7: return tgt::KeyEvent::K_7;
    case Qt::Key_8: return tgt::KeyEvent::K_8;
    case Qt::Key_9: return tgt::KeyEvent::K_9;
    case Qt::Key_Colon: return tgt::KeyEvent::K_COLON;
    case Qt::Key_Semicolon: return tgt::KeyEvent::K_SEMICOLON;
    case Qt::Key_Less: return tgt::KeyEvent::K_LESS;
    case Qt::Key_Equal: return tgt::KeyEvent::K_EQUALS;
    case Qt::Key_Greater: return tgt::KeyEvent::K_GREATER;
    case Qt::Key_Question: return tgt::KeyEvent::K_QUESTION;
    case Qt::Key_At: return tgt::KeyEvent::K_AT;
    case Qt::Key_A: return tgt::KeyEvent::K_A;
    case Qt::Key_B: return tgt::KeyEvent::K_B;
    case Qt::Key_C: return tgt::KeyEvent::K_C;
    case Qt::Key_D: return tgt::KeyEvent::K_D;
    case Qt::Key_E: return tgt::KeyEvent::K_E;
    case Qt::Key_F: return tgt::KeyEvent::K_F;
    case Qt::Key_G: return tgt::KeyEvent::K_G;
    case Qt::Key_H: return tgt::KeyEvent::K_H;
    case Qt::Key_I: return tgt::KeyEvent::K_I;
    case Qt::Key_J: return tgt::KeyEvent::K_J;
    case Qt::Key_K: return tgt::KeyEvent::K_K;
    case Qt::Key_L: return tgt::KeyEvent::K_L;
    case Qt::Key_M: return tgt::KeyEvent::K_M;
    case Qt::Key_N: return tgt::KeyEvent::K_N;
    case Qt::Key_O: return tgt::KeyEvent::K_O;
    case Qt::Key_P: return tgt::KeyEvent::K_P;
    case Qt::Key_Q: return tgt::KeyEvent::K_Q;
    case Qt::Key_R: return tgt::KeyEvent::K_R;
    case Qt::Key_S: return tgt::KeyEvent::K_S;
    case Qt::Key_T: return tgt::KeyEvent::K_T;
    case Qt::Key_U: return tgt::KeyEvent::K_U;
    case Qt::Key_V: return tgt::KeyEvent::K_V;
    case Qt::Key_W: return tgt::KeyEvent::K_W;
    case Qt::Key_X: return tgt::KeyEvent::K_X;
    case Qt::Key_Y: return tgt::KeyEvent::K_Y;
    case Qt::Key_Z: return tgt::KeyEvent::K_Z;
    case Qt::Key_BracketLeft: return tgt::KeyEvent::K_LEFTBRACKET;
    case Qt::Key_Backslash: return tgt::KeyEvent::K_BACKSLASH;
    case Qt::Key_BracketRight: return tgt::KeyEvent::K_RIGHTBRACKET;
    case Qt::Key_Underscore: return tgt::KeyEvent::K_UNDERSCORE;
    default: return tgt::KeyEvent::K_UNKNOWN;
    }
}

int KeyDetectorWidget::getQtKeyFromTGT(tgt::KeyEvent::KeyCode tgtKey) {
    switch (tgtKey) {
    case tgt::KeyEvent::K_ESCAPE: return Qt::Key_Escape;
    case tgt::KeyEvent::K_TAB: return Qt::Key_Tab;
    case tgt::KeyEvent::K_BACKSPACE: return Qt::Key_Backspace;
    case tgt::KeyEvent::K_RETURN: return Qt::Key_Return;
    case tgt::KeyEvent::K_KP_ENTER: return Qt::Key_Enter;
    case tgt::KeyEvent::K_INSERT: return Qt::Key_Insert;
    case tgt::KeyEvent::K_DELETE: return Qt::Key_Delete;
    case tgt::KeyEvent::K_PAUSE: return Qt::Key_Pause;
    case tgt::KeyEvent::K_PRINT: return Qt::Key_Print;
    case tgt::KeyEvent::K_CLEAR: return Qt::Key_Clear;
    case tgt::KeyEvent::K_HOME: return Qt::Key_Home;
    case tgt::KeyEvent::K_END: return Qt::Key_End;
    case tgt::KeyEvent::K_LEFT: return Qt::Key_Left;
    case tgt::KeyEvent::K_UP: return Qt::Key_Up;
    case tgt::KeyEvent::K_RIGHT: return Qt::Key_Right;
    case tgt::KeyEvent::K_DOWN: return Qt::Key_Down;
    case tgt::KeyEvent::K_PAGEUP: return Qt::Key_PageUp;
    case tgt::KeyEvent::K_PAGEDOWN: return Qt::Key_PageDown;
    case tgt::KeyEvent::K_LSHIFT: return Qt::Key_Shift;
    case tgt::KeyEvent::K_RSHIFT: return Qt::Key_Shift;
    case tgt::KeyEvent::K_LCTRL: return Qt::Key_Control;
    case tgt::KeyEvent::K_RCTRL: return Qt::Key_Control;
    case tgt::KeyEvent::K_LMETA: return Qt::Key_Meta;
    case tgt::KeyEvent::K_RMETA: return Qt::Key_Meta;
    case tgt::KeyEvent::K_LALT: return Qt::Key_Alt;
    case tgt::KeyEvent::K_RALT: return Qt::Key_Alt;
    case tgt::KeyEvent::K_MODE: return Qt::Key_AltGr;
    case tgt::KeyEvent::K_CAPSLOCK: return Qt::Key_CapsLock;
    case tgt::KeyEvent::K_NUMLOCK: return Qt::Key_NumLock;
    case tgt::KeyEvent::K_SCROLLOCK: return Qt::Key_ScrollLock;
    case tgt::KeyEvent::K_F1: return Qt::Key_F1;
    case tgt::KeyEvent::K_F2: return Qt::Key_F2;
    case tgt::KeyEvent::K_F3: return Qt::Key_F3;
    case tgt::KeyEvent::K_F4: return Qt::Key_F4;
    case tgt::KeyEvent::K_F5: return Qt::Key_F5;
    case tgt::KeyEvent::K_F6: return Qt::Key_F6;
    case tgt::KeyEvent::K_F7: return Qt::Key_F7;
    case tgt::KeyEvent::K_F8: return Qt::Key_F8;
    case tgt::KeyEvent::K_F9: return Qt::Key_F9;
    case tgt::KeyEvent::K_F10: return Qt::Key_F10;
    case tgt::KeyEvent::K_F11: return Qt::Key_F11;
    case tgt::KeyEvent::K_F12: return Qt::Key_F12;
    case tgt::KeyEvent::K_SPACE: return Qt::Key_Space;
    case tgt::KeyEvent::K_EXCLAIM: return Qt::Key_Exclam;
    case tgt::KeyEvent::K_QUOTEDBL: return Qt::Key_QuoteDbl;
    case tgt::KeyEvent::K_HASH: return Qt::Key_NumberSign;
    case tgt::KeyEvent::K_DOLLAR: return Qt::Key_Dollar;
    case tgt::KeyEvent::K_AMPERSAND: return Qt::Key_Ampersand;
    case tgt::KeyEvent::K_LEFTPAREN: return Qt::Key_ParenLeft;
    case tgt::KeyEvent::K_RIGHTPAREN: return Qt::Key_ParenRight;
    case tgt::KeyEvent::K_ASTERISK: return Qt::Key_Asterisk;
    case tgt::KeyEvent::K_PLUS: return Qt::Key_Plus;
    case tgt::KeyEvent::K_COMMA: return Qt::Key_Comma;
    case tgt::KeyEvent::K_MINUS: return Qt::Key_Minus;
    case tgt::KeyEvent::K_PERIOD: return Qt::Key_Period;
    case tgt::KeyEvent::K_SLASH: return Qt::Key_Slash;
    case tgt::KeyEvent::K_0: return Qt::Key_0;
    case tgt::KeyEvent::K_1: return Qt::Key_1;
    case tgt::KeyEvent::K_2: return Qt::Key_2;
    case tgt::KeyEvent::K_3: return Qt::Key_3;
    case tgt::KeyEvent::K_4: return Qt::Key_4;
    case tgt::KeyEvent::K_5: return Qt::Key_5;
    case tgt::KeyEvent::K_6: return Qt::Key_6;
    case tgt::KeyEvent::K_7: return Qt::Key_7;
    case tgt::KeyEvent::K_8: return Qt::Key_8;
    case tgt::KeyEvent::K_9: return Qt::Key_9;
    case tgt::KeyEvent::K_COLON: return Qt::Key_Colon;
    case tgt::KeyEvent::K_SEMICOLON: return Qt::Key_Semicolon;
    case tgt::KeyEvent::K_LESS: return Qt::Key_Less;
    case tgt::KeyEvent::K_EQUALS: return Qt::Key_Equal;
    case tgt::KeyEvent::K_GREATER: return Qt::Key_Greater;
    case tgt::KeyEvent::K_QUESTION: return Qt::Key_Question;
    case tgt::KeyEvent::K_AT: return Qt::Key_At;
    case tgt::KeyEvent::K_A: return Qt::Key_A;
    case tgt::KeyEvent::K_B: return Qt::Key_B;
    case tgt::KeyEvent::K_C: return Qt::Key_C;
    case tgt::KeyEvent::K_D: return Qt::Key_D;
    case tgt::KeyEvent::K_E: return Qt::Key_E;
    case tgt::KeyEvent::K_F: return Qt::Key_F;
    case tgt::KeyEvent::K_G: return Qt::Key_G;
    case tgt::KeyEvent::K_H: return Qt::Key_H;
    case tgt::KeyEvent::K_I: return Qt::Key_I;
    case tgt::KeyEvent::K_J: return Qt::Key_J;
    case tgt::KeyEvent::K_K: return Qt::Key_K;
    case tgt::KeyEvent::K_L: return Qt::Key_L;
    case tgt::KeyEvent::K_M: return Qt::Key_M;
    case tgt::KeyEvent::K_N: return Qt::Key_N;
    case tgt::KeyEvent::K_O: return Qt::Key_O;
    case tgt::KeyEvent::K_P: return Qt::Key_P;
    case tgt::KeyEvent::K_Q: return Qt::Key_Q;
    case tgt::KeyEvent::K_R: return Qt::Key_R;
    case tgt::KeyEvent::K_S: return Qt::Key_S;
    case tgt::KeyEvent::K_T: return Qt::Key_T;
    case tgt::KeyEvent::K_U: return Qt::Key_U;
    case tgt::KeyEvent::K_V: return Qt::Key_V;
    case tgt::KeyEvent::K_W: return Qt::Key_W;
    case tgt::KeyEvent::K_X: return Qt::Key_X;
    case tgt::KeyEvent::K_Y: return Qt::Key_Y;
    case tgt::KeyEvent::K_Z: return Qt::Key_Z;
    case tgt::KeyEvent::K_LEFTBRACKET: return Qt::Key_BracketLeft;
    case tgt::KeyEvent::K_BACKSLASH: return Qt::Key_Backslash;
    case tgt::KeyEvent::K_RIGHTBRACKET: return Qt::Key_BracketRight;
    case tgt::KeyEvent::K_UNDERSCORE: return Qt::Key_Underscore;
    default: return Qt::Key_Cancel;
    }
}


} // namespace
