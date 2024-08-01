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

#ifndef VRN_TOUCHTABLEMENU_H
#define VRN_TOUCHTABLEMENU_H

#include "tgt/bounds.h"

namespace voreen {
    // widget menu
class  TouchTableMenuFrame {

public:
     TouchTableMenuFrame(tgt::ivec2 anchor, tgt::vec4 color);
    virtual tgt::ivec2 getAnchor() const;
    virtual void setAnchor(tgt::ivec2 anchor);
    virtual tgt::ivec2 getUR() const;
    virtual void setUR(tgt::ivec2 ur);
    virtual tgt::ivec2 getLL() const;
    virtual void setLL(tgt::ivec2 ll);
    virtual tgt::vec4 getColor() const;
    virtual void setColor(tgt::vec4 color);

    /**
    * check if bounds contain point
    *
    * @param tp touchpoint to check
    */
    virtual bool contains(tgt::ivec2 tp);

    /**
    * check if bounds intersect with rectangle
    *
    * @param ll lower left corner of rectangle
    * @param ur upper right corner of rectangle
    */
    virtual bool intersects(tgt::ivec2 ll, tgt::ivec2 ur) const;


protected:
    //position
    tgt::Bounds bounds_;

    //position of the button that triggered the menu
    tgt::ivec2 anchor_;

    //menu color
    tgt::vec4 color_;

};

class TouchTableOverlayMenu : public  TouchTableMenuFrame {
public:
    TouchTableOverlayMenu(tgt::ivec2 anchor, tgt::vec4 color);
    int getID();
    void setID(int id);
    tgt::ivec2 getButton();
    void setButton(tgt::ivec2 button);
    std::vector<tgt::ivec2>* getButtons();
    void addButton(tgt::ivec2 button);
    void clearButtons();

protected:
    //id of touch point currently opening menu
    int id_;

    //position of the button that triggered the menu
    tgt::ivec2 button_;

    //positions of all buttons that activate the menu
    std::vector<tgt::ivec2> buttons_;
};

class TouchTableClippingMenu: public  TouchTableMenuFrame {
public:
    TouchTableClippingMenu(tgt::ivec2 anchor, tgt::vec4 color);
protected:

};

}



#endif

