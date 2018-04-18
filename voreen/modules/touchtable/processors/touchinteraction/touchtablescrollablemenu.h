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

#ifndef VRN_TOUCHTABLESCROLLABLEMENU_H
#define VRN_TOUCHTABLESCROLLABLEMENU_H

#include "touchtablemenuframe.h"
#include "touchtableslider.h"
#include "voreen/core/properties/vectorproperty.h"
#include "tgt/font.h"
#include "tgt/texture.h"

namespace voreen{

class TouchTableScrollableMenu : public  TouchTableMenuFrame {

public:
    TouchTableScrollableMenu(tgt::ivec2 anchor, tgt::vec4 color, tgt::Font* font = 0);

    virtual void setUR(tgt::ivec2 ur);
    virtual void setLL(tgt::ivec2 ll);
    virtual void setContent(std::vector<std::pair<std::string, tgt::Texture*> > content);
    virtual std::vector<std::pair<std::string, tgt::Texture*> > getContent() const;
    virtual void setSelectionHandler(boost::function<void (std::string)>);    ///< handle gets called if content menu got hit with corresponding string as argument
    virtual  TouchTableMenuFrame& getContentMenu();
    virtual TouchTableSlider& getSlider();
    virtual tgt::Font* getFont() const;
    virtual void setFont(tgt::Font* font);
    virtual int getContentOffset() const;
    virtual void setContentOffset(int offset);
    virtual int getMaxElementsInMenu() const;
    virtual int getTextColumnWidth() const;
    virtual void setTextColumnWidth(int width);
    virtual int getRowHeight() const;
    virtual void setRowHeight(int height);
    virtual void setVisualizingSelection(bool selection);
    virtual bool isVisualizingSelection();
    virtual int getSelectedRow();

    /**
    * scrollablemenu handles itself and checks if menu or slider was hit
    *
    * @param tp touchpoints associated with scrollablemenu
    */
    virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

    //places slider depending on margins
    virtual void placeSlider();

    //places content menu depending on margins
    virtual void placeContentMenu();

    ///sets selected row to -1 (ie. no selected row is displayed)
    virtual void clearSelection();

protected:

    virtual void setSelectedRow(int row);

    void handleContentMenu(tgt::vec2 pos);
    void computeContentMenuProperties();

    int maxElementsInMenu_;                    ///< number of elements that can be displayed in menu
    int contentOffset_;                        ///< number of elements above those shown in content menu
    TouchTableSlider menuSlider_;            ///< slider to scroll through menu
    tgt::ivec2 sliderMargin_;                ///< distance of slider to border of  main menu
    tgt::ivec2 menuMargin_;                    ///< distance of menu with content to border of  main menu
     TouchTableMenuFrame contentMenu_;        ///< menu with content to scroll through
    int sliderID_;                            ///< if slider got hit by Touchpoint, this is the TPs ID. Otherwise -1
    std::vector<std::pair<std::string, tgt::Texture*> > content_;        ///< content of menu
    boost::function<void (std::string)> handle_; ///< method that gets called if content in menu got hit by TP
    tgt::Font* font_;
    int textColumnWidth_;                    ///< width of text columns in content menu
    int rowHeight_;                            ///< height of rows in content menu
    bool isVisualizingSelection_;            ///< if true selected row is highlighted
    int selectedRow_;                        ///< selected row that might get highlighted
};

}

#endif
