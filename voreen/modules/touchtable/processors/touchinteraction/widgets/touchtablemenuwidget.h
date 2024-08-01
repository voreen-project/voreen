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

#ifndef VRN_TOUCHTABLEMENUWIDGET_H
#define VRN_TOUCHTABLEMENUWIDGET_H

#include "../touchtableoverlay.h"

namespace voreen {

    class TouchTableMenuWidget : public TouchTableWidget {

        /**
         * struct for keeping track of the menu placement of the widget, which changes if parts of the menu would be placed outside the screen.
         */
        struct MenuPlacement {

            enum MenuPlacementHorizontal {
                LEFT, RIGHT
            };

            enum MenuPlacementVertical {
                UP, DOWN
            };

            MenuPlacement() {
                horizontal_ = RIGHT;
                vertical_ = UP;
            }

            MenuPlacementHorizontal horizontal_;
            MenuPlacementVertical vertical_;
        };

    public:

        /**
         * Constructor
         *
         * @see TouchTableWidget
         * @param textureFile the file (has to be located in voreen/modules/touchtable/textures) that should be used as the widget's symbol.
         */
        TouchTableMenuWidget(std::string textureFile = "nosymbol.png");

        virtual std::string getClassName() const {
            return "TouchTableMenuWidget";
        }

        virtual void initialize();

        virtual void deinitialize();

        /**
         * Sets the current mode widget in the associated overlay processor and calls parent method.
         * @see TouchTableWidget
         * @see TouchTableOverlay
         */
        virtual void pressed();

        /**
         * Checks if the outport size of the overlay has changed and updates menu coordinates if necessary.
         * Then calls parent method.
         * When overwriting this method in subclasses, be sure call the parent's process()-method.
         *
         * @see TouchTableWidget
         */
        virtual void process();


        /**
         * Calls handleTouchPoints(...), which has to be implemented in subclasses.
         * This should not be overwritten!
         */
        virtual void handleWidgetTouchPoints(const std::deque<tgt::TouchPoint>&);


        /**
         * Has to be overwritten in subclasses to handle touch points for the widget.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>&) = 0;

        /**
         * Renders the menu and then calls renderComponents()
         */
        virtual void render();

        /**
         * Called by TouchTableOverlay to check if a touch point (given by its screen coordinates) is located in the menu of the widget.
         * Returns false if the menu is currently not open or the touch point is outside the menu area, true if it is inside.
         *
         * @param tp screen coordinates of the touch point
         */
        virtual bool isInMenu(tgt::ivec2 tp);

        /**
         * Additionally to setting the overlay, this sets itself as the current mode widget of the overlay if the widget is active.
         * If the overlay is 0, the "active" flag is deactivated.
         *
         * @see TouchTableWidget
         */
        virtual void setOverlay(TouchTableOverlay* overlay);

        /// True if the menu orientation is inverted (upside-down)
        virtual bool isMenuInverted() const;

        /// returns the size of the menu
        virtual tgt::ivec2 getMenuDimensions() const;

        /// returns the lower left corner of the menu in canvas coordinates
        virtual tgt::ivec2 getMenuLL() const;

    protected:


        /**
         * Check if a touch point (ie. its coordinates) hits a control element.
         *
         * @param tp the coordinates of the touch point (in menu coordinates)
         * @param element the control element to be checked
         */
        bool hitsControlElement(tgt::ivec2 tp, const TouchTableControlElement& element) const;

        /**
         * This converts coordinates of a touch point from canvas coordinates to menu coordinate space.
         * The method also takes into account if the menu orientation is inverted so that touch points can always be handled in the same coordinate space the components are set.
         */
        virtual tgt::ivec2 convertToMenuCoordinates(tgt::ivec2 p) const;

        /**
         * This takes a list of touch points in screen coordinates and converts them into menu coordinates, also taking into account the menu rotation.
         */
        std::deque<tgt::TouchPoint> convertTouchPointsToMenuCoordinates(const std::deque<tgt::TouchPoint>& points) const;

        /**
         * Is called if properties change to update the attributes of the menu's components (e.g. control elements, submenus, etc.).
         * Has to be implemented in subclasses to place the control elements inside the menu, update their status etc.
         * Properties of subclasses should call this method by using the properties' "onChange" method.
         *
         * Note that component coordinates should be set in menu coordinate space (ie. image space with menu LL being (0,0)).
         */
        virtual void updateComponents() = 0;

        /**
         * Has to be overwritten in subclasses, should contain rendering of control elements, submenus, sliders etc.
         */
        virtual void renderComponents() = 0;

        /**
         * Is called if the widget position or menu dimensions change or the widget is placed inside the overlay menu.
         * Also updates the coordinates of the control elements within the menu.
         * Sets the menu coordinates to the middle of the screen if the widget has been placed in the menu.
         */
        virtual void updateMenuCoordinates();

         TouchTableMenuFrame menu_;                               ///< menu of the mode widget

        IntVec2Property menuDimensions_;                    ///< size of the whole menu of the widget
        IntProperty controlElementRadius_;                  ///< radius of the touch control elements
        BoolProperty invertMenuOrientation_;                ///< if enabled, the menu and its components are rendered inverted regarding their orientation

        tgt::ivec2 lastOutportSize_;                        ///< store the last outport size of the TouchTableOverlay to check if it changed

        std::string symbolFileString_;                      ///< set by constructor, determines the file used as symbol texture

        virtual void updateComponentAttributes();           ///< updates the rotate button, then calls updateComponents()

    private:

        TouchTableMenuFrame rotationFrame_;                 ///< a menu that contains one control element: the rotation button
        TouchTableControlElement rotate_;                   ///< rotates (inverts) the menu orientation
        TouchTableControlElement close_;                    ///< closes the menu widget (alternative to pressing the widget itself)

        tgt::Texture* rotateSymbol_;                        ///< symbol texture for rotate button
        tgt::Texture* closeSymbol_;                         ///< symbol texture for closing the menu

        MenuPlacement menuPlacement_;                       ///< placement of the mode widget menu
        int rotateButtonID_;                                ///< the id of a touch point associated with the rotate button
        int closeButtonID_;                                 ///< the id of a touch point associated with the close button
    };

} // namespace

#endif // VRN_TOUCHTABLEMENUWIDGET_H
