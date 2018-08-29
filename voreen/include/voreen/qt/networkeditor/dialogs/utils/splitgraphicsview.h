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

#ifndef VRN_SPLITGRAPHICSVIEW_H
#define VRN_SPLITGRAPHICSVIEW_H

#include <QGraphicsView>

class QWheelEvent;
class QMouseEvent;
class QKeyEvent;

namespace voreen {

    class PropertyGraphicsItem;
    class PropertyListGraphicsItem;
    class PropertyLinkArrowGraphicsItem;
    class PropertyOwnerGraphicsItem;

/**
 * This class is a QGraphicsView to visualize processor properties in the propertylinkdialog.
 * Its main feaure is the connection to another SplitGraphicsView which allows to scroll propertylists of two
 * processors independent.
 */
class SplitGraphicsView : public QGraphicsView {
Q_OBJECT;
public:
    /**
     * Enum which determines the position of the splitgraphicsview.
     */
    enum GlobalViewAlignment{
        GVA_UNKNOWN,    ///< default alignment. Must be changed before use to work properly
        GVA_CENTER,     ///< This implicates no buddy view. Only one propertylist is in the center of the screen.
        GVA_LEFT,       ///< This implicates the need of an buddy view with GVN_RIGHT. This view is on the left side of the screen.
        GVA_RIGHT       ///< This implicates the need of an buddy view with GVN_LEFT. This view is on the right side of the screen.
    };

public:
    //constructor and destructor
    SplitGraphicsView(QWidget* parent = 0);
    ~SplitGraphicsView();

//---------------------------------------------------------------------------------------------------------------
//                  splitgraphicsview functions
//---------------------------------------------------------------------------------------------------------------
    //-------------------------
    //      getter/setter
    //-------------------------
        //buddy view
    SplitGraphicsView* getBuddyView() const;
    void setBuddyView(SplitGraphicsView* view);
        //global view alignment
    GlobalViewAlignment getGlobalViewAlignment() const;
    void setGlobalViewAlignment(GlobalViewAlignment gva);
        //property list
    PropertyListGraphicsItem* getPropertyList() const;
    void setPropertyList(PropertyListGraphicsItem* list, PropertyOwnerGraphicsItem* destOwner);
protected:
    void updatePropertyListPosition(bool resetVerticalPos);

    //-------------------------
    //      buddylist
    //-------------------------
public:
    /**
     * Updates the viewport and the viewport of the buddylist.
     * Must be called since auto viewport updates are disabled.
     */
    void updateViewport();
public slots:
    /**
     * Updates the selected arrow of the buddylist.
     * Is triggered, if an arrow is selected and the counterpart in the
     * buddylist has to be selected too.
     */
    void updateBuddySelection();
    /**
     * @note is done by the PropertyLinkArrowItems themselves.
     */
    //void setBuddyHoverArrow(PropertyLinkArrowGraphicsItem* pItem);
    /**
     * @note is done by the PropertyLinkArrowItems themselves.
     */
    //void unsetBuddyHoverArrow(PropertyLinkArrowGraphicsItem* pItem);

    //-------------------------
    //  current draged arrow
    //-------------------------
public:
    PropertyLinkArrowGraphicsItem* getCurrentDragedArrow();
    void initCurrentDragedArrow(PropertyGraphicsItem* item, QPointF mousePos);
    void deleteCurrentDragedArrow();
    signals:
    void linkProperties(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);

    //-------------------------
    //      member
    //-------------------------
private:
    SplitGraphicsView* buddyView_;              ///< pointre to the buddy view. Can be null
    GlobalViewAlignment globalViewAlignment_;   ///< actual globalviewalignment
    PropertyListGraphicsItem* propertyList_;    ///< pointer to the propertylist which should be visualized
    PropertyLinkArrowGraphicsItem* currentDragedArrow_; ///< pointer to the new arrow been draged

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
protected:
    //wheel
    void wheelEvent(QWheelEvent* event);
    //mouse
    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    //key
    void keyPressEvent(QKeyEvent* event);
    //resize
    void resizeEvent(QResizeEvent* event);
};

} //namespace voreen

#endif // VRN_SPLITGRAPHICSVIEW_H

