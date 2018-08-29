/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreenveapplication.h"

#include "tgt/filesystem.h"
#include "tgt/logmanager.h"
#include "voreen/core/utils/commandlineparser.h"
#include "voreen/core/version.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/qt/voreenapplicationqt.h"
#include <QMessageBox>

#if !defined(WIN32) && !defined(APPLE) && defined(USE_XINPUT2)

#include <QX11Info>
#include <QTouchEvent>

#include <X11/Xlib.h>
#include <X11/extensions/XInput.h>
#include <X11/extensions/XInput2.h>
#include <X11/Xutil.h>

#include <linux/input.h>
#include <xorg/xserver-properties.h>

#endif

namespace voreen {

VoreenVEApplication* VoreenVEApplication::veApp_ = 0;
const std::string VoreenVEApplication::loggerCat_("voreenve.VoreenVEApplication");

VoreenVEApplication::VoreenVEApplication(int& argc, char** argv)
    : QApplication(argc, argv)
    , VoreenApplicationQt("voreenve", "VoreenVE", "Voreen Visualization Environment", argc, argv, VoreenApplication::APP_ALL)
#if !defined(WIN32) && !defined(APPLE) && defined(USE_XINPUT2)
    , display_(0)
    , xi_opcode_(-1)
    , minMaxRetrieved_(false)
#endif
{
    veApp_ = this;
}

VoreenVEApplication::~VoreenVEApplication() {
}

VoreenVEApplication* VoreenVEApplication::veApp() {
    return veApp_;
}

bool VoreenVEApplication::notify(QObject* receiver, QEvent* event) {
    bool result = false;
    try {
        result = QApplication::notify(receiver, event);
    }
    catch (const VoreenException& e) {
        LERRORC("voreenve.main", "Caught unhandled VoreenException: " << e.what());
#ifndef TGT_NON_INTERACTIVE_ASSERT
        int choice = QMessageBox::critical(0, tr("VoreenVE"), tr("Caught unhandled VoreenException:\n\"")
                                           + e.what() + +"\"\n" + tr("Continue?"),
                                           QMessageBox::Ok | QMessageBox::Cancel);
        if (choice == QMessageBox::Cancel) {
  #ifdef VRN_DEBUG
            TGT_THROW_BREAKPOINT;
  #else
            exit(1);
  #endif
        }
#else
        exit(1);
#endif // TGT_NON_INTERACTIVE_ASSERT
    }
    catch (const std::exception& e) {
        LERRORC("voreenve.main", "Caught unhandled std::exception: " << e.what());
#ifndef TGT_NON_INTERACTIVE_ASSERT
        int choice = QMessageBox::critical(0, tr("VoreenVE"), tr("Caught unhandled std::exception:\n\"")
                                           + e.what() + "\"\n" + tr("Continue?"),
                                           QMessageBox::Ok | QMessageBox::Cancel);
        if (choice == QMessageBox::Cancel) {
#ifdef VRN_DEBUG
            TGT_THROW_BREAKPOINT;
#else
            exit(1);
#endif
        }
#else
            exit(1);
#endif // TGT_NON_INTERACTIVE_ASSERT
    }
    catch (...) {
        LERRORC("voreenve.main", "Caught unhandled unknown exception!");
#ifndef TGT_NON_INTERACTIVE_ASSERT
        int choice = QMessageBox::critical(0, tr("VoreenVE"), tr("Caught unhandled unknown exception!\nContinue?"),
                                           QMessageBox::Ok | QMessageBox::Cancel);
        if (choice == QMessageBox::Cancel) {
#ifdef VRN_DEBUG
            TGT_THROW_BREAKPOINT;
#else
            exit(1);
#endif
        }
#else
        exit(1);
#endif // TGT_NON_INTERACTIVE_ASSERT
        throw;
    }
    return result;
}

std::string VoreenVEApplication::getApplicationResourcePath(const std::string& filename) const {
    return tgt::FileSystem::cleanupPath(getBasePath() + "/resource/voreenve" + (filename.empty() ? "" : "/" + filename));
}

#if !defined(WIN32) && !defined(APPLE) && defined(USE_XINPUT2)

// These functions add experimental Linux Multitouch support.  The code was largely adapted from
// http://gabrbedd.wordpress.com/2012/07/10/getting-multitouch-qt-on-ubuntu-12-04/.
void VoreenVEApplication::initXinput() {
    if (display_)
        return;

    int event, error;
    display_ = QX11Info::display();

    if (!XQueryExtension(display_, "XInputExtension", &xi_opcode_, &event, &error))
        throw std::runtime_error("Could not find XInputExtension");
}

void VoreenVEApplication::sendTouchEventsTo(QWidget *w) {
    XIEventMask mask;
    Window win;
    Display *dpy;

    dpy = QX11Info::display();
    win = w->winId();

    memset(&mask, 0, sizeof(XIEventMask));
    mask.deviceid = XIAllMasterDevices;
    mask.mask_len = XIMaskLen(XI_LASTEVENT);
    mask.mask = (unsigned char*) calloc(mask.mask_len, sizeof(char));

    XISetMask(mask.mask, XI_TouchBegin);
    XISetMask(mask.mask, XI_TouchUpdate);
    XISetMask(mask.mask, XI_TouchEnd);

    XISelectEvents(dpy, win, &mask, 1);

    free(mask.mask);
}

bool VoreenVEApplication::handleX11TouchEvent(XEvent *ev) {
    XGenericEventCookie *cookie = &ev->xcookie;

    QEvent::Type eventType;
    switch (cookie->evtype) {
        case XI_TouchBegin: eventType = QEvent::TouchBegin; break;
        case XI_TouchUpdate: eventType = QEvent::TouchUpdate; break;
        case XI_TouchEnd: eventType = QEvent::TouchEnd; break;
        default: return false;
    }

    XIDeviceEvent *de = reinterpret_cast<XIDeviceEvent*>(cookie->data);
    const char *evname = 0;
    QWidget *target;

    const double *value;
    const unsigned char *mask;
    int touchID, nbyte, nbit;

    value = de->valuators.values;
    mask = de->valuators.mask;
    touchID = de->detail;
    nbyte = de->valuators.mask_len;

    double posX = 0.0;
    double posY = 0.0;
    bool foundCoords = false;

    // TODO This assumes that the coords are always at position 0 and 1, which seems to be the case for the magic trackpad,
    // but not necessarily for every device...
    if(nbyte > 0 && XIMaskIsSet(mask, 0) && XIMaskIsSet(mask, 1)) {
        posX = value[0];
        posY = value[1];
        foundCoords = true;
    }

    QTouchEvent::TouchPoint touchPoint;
    if(!minMaxRetrieved_) {
        int ndevices;
        int id = de->deviceid;
        XIDeviceInfo* info = XIQueryDevice(display_, XIAllDevices, &ndevices);
        for (int i = 0; i < ndevices; i++) {
            if (info[i].deviceid == id) {
                XIDeviceInfo curInf = info[i];
                bool foundXAxis = false;
                bool foundYAxis = false;
                for(int j = 0; j < curInf.num_classes; j++) {
                    if(curInf.classes[j]->type == XIValuatorClass) {
                        XIValuatorClassInfo* val = reinterpret_cast<XIValuatorClassInfo*>(curInf.classes[j]);
                        char* name = XGetAtomName(display_, val->label);
                        if(!strcmp(name, AXIS_LABEL_PROP_ABS_MT_POSITION_X)) {
                            foundXAxis = true;
                            minEvX_ = val->min;
                            maxEvX_ = val->max;
                        }
                        if(!strcmp(name, AXIS_LABEL_PROP_ABS_MT_POSITION_Y)) {
                            foundYAxis = true;
                            minEvY_ = val->min;
                            maxEvY_ = val->max;
                        }
                    }
                }
                if(foundXAxis && foundYAxis)
                    minMaxRetrieved_ = true;
                break;
            }
        }

        XIFreeDeviceInfo(info);
        if(!minMaxRetrieved_)
            LWARNING("Could not initialize touch event device coordinate mapping");
    }

    target = QWidget::find(de->event);

    // TODO this re-maps the event coords to the target widget, which will not always be necessary
    if(foundCoords && minMaxRetrieved_) {
        QSize targetSize = target->size();
        posX = ((double)targetSize.width())  * ((posX - minEvX_) / (maxEvX_ - minEvX_));
        posY = ((double)targetSize.height()) * ((posY - minEvY_) / (maxEvY_ - minEvY_));
        touchPoint.setPos(QPoint((int)posX, (int)posY));
    }

    // TODO
    touchPoint.setPressure(1.0);

    // TODO
    //touchPoint.setRect(target->rect());
    //touchPoint.setScreenRect(screenRect);
    //touchPoint.setScreenPos(screenPos);

    Qt::TouchPointStates touchPointStates = 0;
    bool pointRemoved = false;

    switch (eventType) {
        case QEvent::TouchBegin: {
            touchPointStates = Qt::TouchPointPressed;

            // TODO
            //touchPoint.setStartNormalizedPos(normPos);
            //touchPoint.setStartPos(touchPoint.pos());
            //touchPoint.setStartScreenPos(screenPos);
            //touchPoint.setStartScenePos(touchPoint.scenePos());

            //touchPoint.setLastNormalizedPos(normPos);
            //touchPoint.setLastPos(touchPoint.pos());
            //touchPoint.setLastScreenPos(screenPos);
            //touchPoint.setLastScenePos(touchPoint.scenePos());

            touchPoint.setState(Qt::TouchPointPressed);
            touchPointMap_.insert(touchID, touchPoint);
            break;
        }
        case QEvent::TouchUpdate: {
            // TODO
            //if (tcur->getMotionSpeed() > 0)
            //else
                //touchPointStates = Qt::TouchPointStationary;
                //touchPoint.setState(Qt::TouchPointStationary);

            // TODO
            //touchPoint.setStartNormalizedPos(touchPointMap_->value(tcur->getSessionID()).startNormalizedPos());
            //touchPoint.setStartPos(touchPointMap_->value(tcur->getSessionID()).startPos());
            //touchPoint.setStartScreenPos(touchPointMap_->value(tcur->getSessionID()).startScreenPos());
            //touchPoint.setStartScenePos(touchPointMap_->value(tcur->getSessionID()).startScenePos());

            //touchPoint.setLastNormalizedPos(touchPointMap_->value(tcur->getSessionID()).normalizedPos());
            //touchPoint.setLastPos(touchPointMap_->value(tcur->getSessionID()).pos());
            //touchPoint.setLastScreenPos(touchPointMap_->value(tcur->getSessionID()).screenPos());
            //touchPoint.setLastScenePos(touchPointMap_->value(tcur->getSessionID()).scenePos());

            if(touchPointMap_.contains(touchID) && foundCoords) {
                touchPointStates = Qt::TouchPointMoved;
                touchPoint.setState(Qt::TouchPointMoved);
                touchPointMap_[touchID] = touchPoint;
            }
            break;
        }
        case QEvent::TouchEnd: {
            touchPointStates = Qt::TouchPointReleased;

            // TODO
            //touchPoint.setStartNormalizedPos(touchPointMap_->value(tcur->getSessionID()).startNormalizedPos());
            //touchPoint.setStartPos(touchPointMap_->value(tcur->getSessionID()).startPos());
            //touchPoint.setStartScreenPos(touchPointMap_->value(tcur->getSessionID()).startScreenPos());
            //touchPoint.setStartScenePos(touchPointMap_->value(tcur->getSessionID()).startScenePos());

            //touchPoint.setLastNormalizedPos(touchPointMap_->value(tcur->getSessionID()).normalizedPos());
            //touchPoint.setLastPos(touchPointMap_->value(tcur->getSessionID()).pos());
            //touchPoint.setLastScreenPos(touchPointMap_->value(tcur->getSessionID()).screenPos());
            //touchPoint.setLastScenePos(touchPointMap_->value(tcur->getSessionID()).scenePos());

            if(touchPointMap_.contains(touchID)) {
                touchPointMap_[touchID].setState(Qt::TouchPointReleased);
                pointRemoved = true;
            }
            break;
        }
        default: {}
    }

    QList<QTouchEvent::TouchPoint> values = touchPointMap_.values();
    if(pointRemoved)
        touchPointMap_.remove(touchID);
    foreach(QTouchEvent::TouchPoint t, values) {
        touchPointStates |= t.state();
    }
    QEvent* touchEvent = new QTouchEvent(eventType, QTouchEvent::TouchScreen, Qt::NoModifier, touchPointStates, values);
    for(QMap<int, QTouchEvent::TouchPoint>::iterator it =  touchPointMap_.begin(); it != touchPointMap_.end(); it++) {
        if(it->state() == Qt::TouchPointMoved)
            it->setState(Qt::TouchPointStationary);
    }

    postEvent(target, touchEvent);
    return true;
}

bool VoreenVEApplication::x11EventFilter(XEvent *ev) {
    bool handled = false;
    XGenericEventCookie *cookie;

    if (!display_)
        initXinput();

    // All XI2 events are GenericEvents
    if (ev->type != GenericEvent)
        return false;

    cookie = &ev->xcookie;

    if (cookie->extension != xi_opcode_)
        return false;

    if (!cookie->data)
        XGetEventData(cookie->display, cookie);

    switch (cookie->evtype) {
        case XI_TouchBegin:
        case XI_TouchUpdate:
        case XI_TouchEnd:
            handled = handleX11TouchEvent(ev);
            break;
        default:
            break;
    }

    if (handled)
        XFreeEventData(cookie->display, cookie);

    return handled;
}

#endif // NOT WIN32 && NOT APPLE && USE_XINPUT

} // namespace
