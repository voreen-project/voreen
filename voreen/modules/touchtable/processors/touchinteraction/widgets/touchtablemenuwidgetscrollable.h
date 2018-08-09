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

#ifndef VRN_TOUCHTABLEMENUWIDGETSCROLLABLE_H
#define VRN_TOUCHTABLEMENUWIDGETSCROLLABLE_H

#include "touchtablemenuwidget.h"
#include "voreen/core/properties/fontproperty.h"

namespace voreen {

    class TouchTableMenuWidgetScrollable : public TouchTableMenuWidget {

    public:

        /**
         * Constructor
         *
         * @see TouchTableMenuWidget
         */
        TouchTableMenuWidgetScrollable(std::string textureFile = "nosymbol.png");

        virtual std::string getClassName() const {
            return "TouchTableMenuWidgetScrollable";
        }

        virtual void initialize();

        virtual void deinitialize();

        /**
         * Do not overwrite this method! Calls handleTouchPoints(...), which has to be implemented in subclasses.
         * This automatically handles touchpoints of the scrollabe menu.
         * @see TouchTableMenuWidget
         */
        virtual void handleWidgetTouchPoints(const std::deque<tgt::TouchPoint>&);

        /**
         * Renders the menu and then calls renderComponents().
         */
        virtual void render();

        /**
         * Has to be overwritten in subclasses to handle touch points for the widget.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>&) = 0;

        /**
         * Has to be overwritten in subclasses. Is called from ScrollableMenu if an entry has been selected.
         */
        virtual void handleScrollableMenuSelection(std::string file) = 0;

    protected:

        /**
         * @see TouchTableMenuWidget
         */
        virtual void updateComponentAttributes();

        /**
         * Is called if properties change to update the attributes of the menu's components (e.g. control elements, submenus, etc.).
         * Has to be implemented in subclasses to place the control elements inside the menu, update their status etc.
         * Properties of subclasses should call this method by using the properties' "onChange" method.
         *
         * Note that component coordinates should be set in menu coordinate space (ie. image space with menu LL being (0,0)).
         * The scrollable menu is updated separately (method updateScrollableMenuPosition).
         */
        virtual void updateComponents() = 0;

        /**
         * Has to be overwritten in subclasses, should contain rendering of control elements, submenus, sliders etc.
         * Should NOT contain redering of scrollable menu as this is handled automatically.
         */
        virtual void renderComponents() = 0;

        /**
         * Updates the position of the scrollable menu and attributes like lineWidth etc.
         * Has to be implemented in subclasses.
         */
        virtual void updateScrollableMenuPosition() = 0;

        /// true by default, may be used to close the scrollable menu if needed in subclasses
        BoolProperty scrollableMenuIsOpen_;

        FontProperty fontProp_;                             ///< property for font rendering of scrollable menu

        TouchTableScrollableMenu scrollableMenu_;

        std::vector<int> scrollableMenuIDs_;                ///< remember touch point ids that are associated with the scrollable menu

    };

} // namespace

#endif // VRN_TOUCHTABLEMENUWIDGETSCROLLABLE_H
