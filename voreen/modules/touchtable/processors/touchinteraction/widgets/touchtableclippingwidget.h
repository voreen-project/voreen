/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_TOUCHTABLECLIPPINGWIDGET_H
#define VRN_TOUCHTABLECLIPPINGWIDGET_H

#include "touchtablemenuwidget.h"
#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "../../../touchtablemodule.h"

namespace voreen {


    /**
    * A struct that contains all properties according to a ClippingSet
    *
    * @memb color_: color property of the whole ClippingSet
    * @memb renderBoolProperty_: bool property for whether or not the clipping plane should be rendered
    * @memb toggleBoolProperty_: bool property for whether or not the clipping is active
    * @memb planeNormalProperty_: FloatVec3 property for the normal of the clipping plane
    * @memb planePositionProperty_: Float property for the distance of the clipping plane from the center
    */
    struct ClippingPropertySet {

        ClippingPropertySet(int i)
           : color_(std::string("clipplane" + itos(i) + "color"), std::string("Clip Plane " + itos(i) + " Color"), tgt::vec4(1.f, 0.f, 0.f, 1.f))
           , renderBoolProperty_(std::string("clipplane" + itos(i) + "render"), std::string("Clip Plane " + itos(i) + " Render Geometry"), false)
           , renderManipulation_(std::string("clipplane" + itos(i) + "manipulation"), std::string("Clip Plane " + itos(i) + " render manipulation"), false)
           , toggleBoolProperty_(std::string("clipplane" + itos(i) + "toggle"), std::string("Clip Plane " + itos(i) + " Toogle Plane"), false)
           , planeNormalProperty_(std::string("clipplane" + itos(i) + "planeNormal"), std::string("Clip Plane " + itos(i) + "Plane Normal"), tgt::vec3(0.f,1.f,0.f), tgt::vec3(-1.f), tgt::vec3(1.f), false)
           , planePositionProperty_(std::string("clipplane" + itos(i) + "planeposition"), std::string("Clip Plane " + itos(i) + "Plane Position"), 0.f,-FLT_MAX , FLT_MAX, false)
           , invertProperty_(std::string("clipplane" + itos(i) + "invert"), std::string("Clip Plane " + itos(i) + "Invert"))
        {
            planePositionProperty_.setVisibleFlag(false);
            planeNormalProperty_.setVisibleFlag(false);
        }

        FloatVec4Property color_;
        BoolProperty renderBoolProperty_;
        BoolProperty renderManipulation_;
        BoolProperty toggleBoolProperty_;
        FloatVec3Property planeNormalProperty_;
        BoolProperty invertProperty_;
        FloatProperty planePositionProperty_;
    };

    /**
    * A struct that contains all control elements for one clipping plane
    *
    * @memb toggleClippingPlane_: control element to (de)activate the clipping
    * @memb renderClippingPlaneGeometry_: control element to (de)activate the rendering of the clippnig plane
    * @memb invertClippingPlane_: control element to invert the clipped geometry on the clipping plane
    * @memb freeClipping_: control element to activate the drawing of a free clipping plane
    * @memb deletePlane_: control element to delet the clipping set
    * @memb slider_: slider to change the plane distance by sliding the indicator
    * @memb menu_: menu that contains the clipping set
    * @memb propertySet_: PropertySet containing associated clipping properties
    */

    struct ClippingSet{
        /**
        * ClippingSet containig all control elements associated with one clipping plane
        *
        * @param anchor: lower left position of the menu the clipping set is located in
        * @param posLL: starting point at which the control elements are being rendered
        * @param color: initial color of the control elements
        * @param radius: initial raduis of the control elements
        * @param propertySet: PropertySet with the linked clipping properties
        * @param posInMenu: initial position in the clipping menu containing all smaller clippingset menus
        * @param sliderLength: initial slider length
        */
        ClippingSet(tgt::ivec2 anchor, tgt::ivec2 posLL, tgt::vec4 color, int radius, ClippingPropertySet* propertySet, int posInMenu, int sliderLength, int id)
            : menu_(anchor, color + tgt::vec4(0.f, 0.f, 0.f, -0.5f))
            , toggleClippingPlane_(radius, posLL+tgt::ivec2(5+radius, 5+radius), 0, 1.f, color, color, id, true)
            , renderClippingPlaneGeometry_(radius, posLL+tgt::ivec2(10+3*radius, 5+radius), 0, 1.f, color, color,id, true)
            , renderManipulation_(radius, posLL + tgt::ivec2(15+5 * radius, 5+radius), 0 , 1.f, color, color, id, true)
            , invertClippingPlane_(radius, posLL+tgt::ivec2(20+7*radius, 5+radius), 0, 1.f, color, color, id, true)
            , freeClipping_(radius, posLL+tgt::ivec2(25+9*radius, 5+radius),0, 1.f, color, color, id, true, true)
            , deletePlane_(radius, posLL+tgt::ivec2(30+11*radius, 5+radius),0, 1.f, color, color,id , true, false)
            , slider_(posLL+tgt::ivec2(35+12*radius, 5+radius), sliderLength, 13, HORIZONTAL, 30, 10)
            , propertySet_(propertySet)
        {
            menu_.setLL(posLL + tgt::ivec2(0,40));
            menu_.setUR(menu_.getLL()+tgt::ivec2(40+6*radius,10+2*radius));
        }

        TouchTableControlElement toggleClippingPlane_;
        TouchTableControlElement renderClippingPlaneGeometry_;
        TouchTableControlElement renderManipulation_;
        TouchTableControlElement invertClippingPlane_;
        TouchTableControlElement freeClipping_;
        TouchTableControlElement deletePlane_;
        TouchTableSlider slider_;

         TouchTableMenuFrame menu_;
        ClippingPropertySet* propertySet_;

    };

    /*
    * A struct that contains all parameters for a new clipping plane drawn by free clipping
    *
    * @memb: id_: -1 if clipping plane has not been drawn otherwise equals the value of the associated touchpoint
    * @memb: startPos_: starting position of the drawn line
    * @memb: currenPos_: current position of the touchpoint drawing the line
    * @memb: activeClippingSet_: clipping set calling the free clipping
    */
    struct freeClippingParameters{

        freeClippingParameters(int id, tgt::ivec2 startPos, tgt::ivec2 currentPos)
            :id_(id)
            , startPos_(startPos)
            , currentPos_(currentPos)
            , activeClippingSet_(0)
        {

        }

        int id_;
        tgt::ivec2 startPos_;
        tgt::ivec2 currentPos_;
        ClippingSet* activeClippingSet_;
    };


    /**
    * TouchTableClippingWidget inherits from TouchTableMenuWidget. When active through this widget
    * settings of a linked clipping processor can be configured
    */
    class TouchTableClippingWidget : public TouchTableMenuWidget {

    public:

        TouchTableClippingWidget();

        ~TouchTableClippingWidget();

        Processor* create() const {
            return new TouchTableClippingWidget();
        }

        std::string getClassName() const {
            return "TouchTableClippingWidget";
        }

        virtual void initialize();
        virtual void deinitialize();

        /**
         * Handles touch points within the clipping overlay menu
         *
         * @param tp: a std::deque containing all the touch points that should be handled by the widget.
         */
        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>&);

    protected:

        void setDescriptions(){
            setDescription("When active through this widget settings of a linked clipping processor can be configured, such as drawing and configuring clipping planes and toggling the bounding box");
        }

         /**
         * Renders the clipping sets, boundingbox control element, etc.
         */
        virtual void renderComponents();

        /**
        * sets the textures for the control elements in a clipping set
        *@param texClipping texture pointer of the clipping control element
        *@param texRenderGeometry texture pointer of the render geometry control element
        *@param texInvert texture pointer of the insert control element
        *@param freeClipping texture pointer of the free clipping control element
        *@param texDelete texture pointer of the delete control element
        *@param texRenderManipulation pointer of the render manipulation control element
        */
        virtual void setSymbolTexClippingSet(tgt::Texture* texClipping, tgt::Texture* texRenderGeometry, tgt::Texture* texInvert, tgt::Texture* freeClipping, tgt::Texture* texDelete, tgt::Texture* texRenderManipulation);

        virtual void updateComponents();

        /**
         * overwritten from TouchTableMenuWidget
         */
        virtual void updateMenuCoordinates();

        /**
        *sets all properties and control elements on false, that are not supposed to be shown when widget inactive
        */
        virtual void deactivateWidget();

        /*****Functions called when properties change******/
        /**
        * updates the clipping set color accoring to the color property
        */
        void updateClippingSetColor();
        /**
        * updates the render control elements according to the render geometry properties from each clipping set
        */
        void updateRenderControlElements();
        /**
        * updates the toggle control elements accoring to the toggle properties from eacht clipping set
        */
        void updateToggleControlElements();
        /**
        * updates the slider position of every slider in each clipping set
        */
        void updateSliderPosition();
        /**************************************************/

        /*****Functions called when touchpoints occur in the clipping menu of the widget*****/
        /*
        * handles released touch points in normal mode(no free clipping active)
        *@param tpPos position of current touch point converted to clipping menu coordinates
        *@param tpId id of the current touch point
        */
        void handleTouchPointReleased(tgt::ivec2 tpPos, int tpId);

        /**
        * handles pressed touch points in normal mode(no free clipping active)
        *@param tp current touch point
        */
        void handleTouchPointPressed(const tgt::TouchPoint tp);

        /**
        * handles moved touch points in normal mode(no free clipping active)
        *@param current touch point
        */
        void handleTouchPointMoved(const tgt::TouchPoint tp);

        /**
        * handles pressed touch points in exculsive mode, when free clipping line is to be drawn
        *@param tp curren touch point
        */
        void handleExclusiveTouchPointPressed(const tgt::TouchPoint& tp);

        /**
        * handles moved touch points in exculsive mode, when free clipping line is to be drawn
        *@param tp curren touch point
        */
        void handleExclusiveTouchPointMoved(const tgt::TouchPoint& tp);

        /**
        * handles released touch points in exculsive mode, when free clipping line is to be drawn
        *@param tp curren touch point
        */
        void handleExclusiveTouchPointReleased(const tgt::TouchPoint& tp);
        /****************************************************************************************/

        /**
        * updates the position of the clipping set within the clipping menu, when a set is added/deleted
        *
        * @param clippingSet: clippingSet whose position is to be updated
        * @param posInMenu: new position for the clippningSet in the clipping menu
        */
        void updateClippingSetPosition(ClippingSet& clippingSet, int posInMenu);

        /**
        * updates the properties according to the free clipping line that has been drawn
        *
        * @param startPos: starting position of the free clipping line
        * @param endPos: end position of the free clipping line
        * @param planePosition: plane distance of the clipping plane that needs to be updated
        * @param planeNormal: plane normal of the clipping plane that needs to be updated
        */
        void setClippingPropsFreeClipping(tgt::ivec2 startPos, tgt::ivec2 endPos,FloatProperty& planePosition,  FloatVec3Property& planeNormal);


        /**
        * renders the clipping line with its given properties
        *@param currentPos position of the touch point drawing the clipping line
        *@param startingPos position of the starting touch point drawing the clipping line
        *@param lineColor color of the clipping line
        */
        void renderClippingLine(tgt::ivec2 currentPos, tgt::ivec2 startingPos, tgt::vec4 lineColor);

        /**
         * Checks the scene bounds of the camera and computes the minimal and maximal distance from origin for clipping planes.
         */
        void updateSliderRange();

        std::map<int, ClippingSet*> touchPointSliderMapping_;            ///< map binding clipping set and the touchpoints associated with the slider in the clipping set

        BoolProperty boundingBoxOn_;                                    ///< bool: true if boundingbox is shown, false otherwise
        CameraProperty cameraProp_;                                        ///< camera: managing position etc. of the camera
        std::vector<ClippingPropertySet*> propertySets_;                ///< vector of all properties according to one clipping set

        std::vector<ClippingSet> clippingSets_;                            ///< list of all clipping sets with the control elements, shown for the user
        std::vector<ClippingSet> inactiveClippingSets_;                    ///< list of all clipping sets with the control elements, hidden for the user
        TouchTableControlElement contrelemBoundingBox_;                    ///< control element to (de)activate the bounding box
        TouchTableControlElement contrelemAdd_;                            ///< control element to add a new clipping plane by free clipping
        TouchTableControlElement contrelemAddXY_;                        ///< control element to add a new clipping plane in center of the geometry with Z-Achsis as planeNormal
        TouchTableControlElement contrelemAddYZ_;                        ///< control element to add a new clipping plane in center of the geometry with X-Achsis as planeNormal
        TouchTableControlElement contrelemAddXZ_;                        ///< control element to add a new clipping plane in center of the geometry with Y-Achsis as planeNormal

        freeClippingParameters freeClippingLineParams_;                    ///< saves touchpoint id, current position of tp and starting point where the touchpoint was pressed

        bool exclusiveModeSet_;                                            ///< bool true if exclusive mode is set, false otherwi

        tgt::Texture* symbolTexClipping_;                                ///< texture pointer of clipping control element
        tgt::Texture* symbolTexBounding_;                                ///< texture pointer of bounding box control element
        tgt::Texture* symbolTexRenderGeometry_;                            ///< texture pointer of render geometry control element
        tgt::Texture* symbolTexRenderManipulation_;                        ///< texture pointer of render manipulation control element
        tgt::Texture* symbolTexInvert_;                                    ///< texture pointer of invert control element
        tgt::Texture* symbolTexDelete_;                                    ///< texture pointer of delete control element
        tgt::Texture* symbolTexFreeClipping_;                            ///< texture pointer of free clipping control element
        tgt::Texture* symbolTexAdd_;                                    ///< texture pointer of add control element
        tgt::Texture* symbolTexAddXY_;                                    ///< texture pointer of addXY control element
        tgt::Texture* symbolTexAddXZ_;                                    ///< texture pointer of addXZ control element
        tgt::Texture* symbolTexAddYZ_;                                    ///< texture pointer of addYZ control element

        int sliderLength_;                                                ///< length of slider to modify plane positions

        FloatVec2Property sliderRange_;                                    ///< the min and max values of the slider for clip plane positions

        IntProperty maxPlanes_;                                         ///< the maximum number of planes that are available in the GUI

        private:

        static const std::string loggerCat_;
    };

} // namespace

#endif // VRN_TOUCHTABLECLIPPINGWIDGET_H
