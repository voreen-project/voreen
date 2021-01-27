/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_TOUCHTABLEOVERLAY_H
#define VRN_TOUCHTABLEOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "voreen/core/ports/genericcoprocessorport.h"
#include "voreen/core/voreenapplication.h"

#include "modules/touchtable/datastructures/geometry/trianglemeshgeometrysinglecolor.h"

#include "tgt/texturemanager.h"

#include "widgets/touchtablewidget.h"
#include "touchtablemenuframe.h"
#include "touchtablecontrolelement.h"
#include "touchtableslider.h"
#include "touchtablescrollablemenu.h"
#include "touchtablescrollabletimelinemenu.h"

namespace voreen {

class TouchTableMenuWidget;

/**
 * Handle touch events and render touch widget overlays ("pucks").
 */
class TouchTableOverlay : public ImageProcessor {
public:
    TouchTableOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "TouchTableOverlay";     }
    virtual std::string getCategory() const   { return "Touch Table"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

    bool isReady() const;

    void setMenuWidget(TouchTableMenuWidget* widget);

    virtual void onEvent(tgt::Event* e);

    void renderMenuFrame(const  TouchTableMenuFrame& menu, bool useMenuWidgetOrientationInvert = true, tgt::ivec2 offset = tgt::ivec2(0,0));

    /**
     * Renders a control element
     *
     * @param controlelement the control element to be rendered, texture has to be set to avoid crashes
     * @param menuLL offset (lower left corner of the menu the control element is placed in)
     * @param greyOutFactor factor that determines if (and how much) the control element is greyed out (1: completely grey, 0: not greyed out, between is linearly interpolated)
     */
    void renderControlElement(TouchTableControlElement* controlelement, tgt::ivec2 menuLL, float greyOutFactor = 0.f, bool useMenuWidgetOrientationInvert = true);

    /**
     * Renders a slider with two indicators
     *
     * @param slider slider to be rendered
     * @param menuLL offset (lower left corner of the menu the control element is placed in)
     * @param greyOutFactor factor that determines if (and how much) the control element is greyed out (1: completely grey, 0: not greyed out, between is linearly interpolated)
     */
    void renderSlider(TouchTableTwoIndicatorSlider* slider, tgt::ivec2 menuLL, float greyOutFactor = 0.f);

    /**
     * Renders a slider with one indicators
     *
     * @param slider slider to be rendered
     * @param menuLL offset (lower left corner of the menu the control element is placed in)
     * @param greyOutFactor factor that determines if (and how much) the control element is greyed out (1: completely grey, 0: not greyed out, between is linearly interpolated)
     */
    void renderSlider(TouchTableSlider* slider, tgt::ivec2 menuLL, float greyOutFactor = 0.f);

    /**
     * Renders a TouchTableScrollabeMenu
     *
     * @param menu the menu to be rendered
     * @param offset the offset in screen coordinates (lower left corner of the widget menu where the scrollable menu is placed in
     */
    void renderScrollableMenu(TouchTableScrollableMenu* menu, tgt::ivec2 offset = tgt::ivec2(0,0));


    void setExclusiveMode(bool mode=true);

    tgt::ivec2 getOutportSize() const;

    int getWidgetRadius() const;

    static const tgt::vec4 getGlobalMenuColor();

    bool getCameraMovement() const;

    void setCameraMovement(bool move = true);


    /**
     * Renders a quad with the input texture
     *
     *
     * @param ll the lower left of the quad (in screen coordinates)
     * @param ur the upper right (in screen coordinates)
     * @param tex the texture
     * @param colorMod a color modificator that is simply added to the texture color
     * @param opacity an opcity factor that is multiplied with the opacity of the color (after adding the colorMod)
     * @param greyOutFactor sets a float (clamped to a value between 0 and 1) for greying out the quad (1: completely greyed out, 0: not greyed out, between is linearly interpolated)
     * @param useMenuWidgetOrientationInvert if set to true menu and is rendered upside-down
     */
    void renderTexturedQuad(tgt::ivec2 ll, tgt::ivec2 ur, const tgt::Texture* tex, tgt::vec4 colorMod = tgt::vec4(0.f), float opacity = 1.f, const tgt::mat4& textureMatrix = tgt::mat4::identity, float greyOutFactor = 0.f, bool useMenuWidgetOrientationInvert = true);

    /**
     * Renders a one dimensional quad with the input texture
     *
     *
     * @param ll the lower left of the quad (in screen coordinates)
     * @param ur the upper right (in screen coordinates)
     * @param tex the texture
     * @param colorMod a color modificator that is simply added to the texture color
     * @param opacity an opcity factor that is multiplied with the opacity of the color (after adding the colorMod)
     * @param greyOutFactor sets a float (clamped to a value between 0 and 1) for greying out the quad (1: completely greyed out, 0: not greyed out, between is linearly interpolated)
     * @param useMenuWidgetOrientationInvert if set to true menu and is rendered upside-down
     */
    void renderTexturedQuad1DTexture(tgt::ivec2 ll, tgt::ivec2 ur, const tgt::Texture* tex, tgt::vec4 colorMod = tgt::vec4(0.f), float opacity = 1.f, float greyOutFactor = 0.f, bool useMenuWidgetOrientationInvert = true);

    /**
     * Renders a quad in screen coordinates with the input color
     */
    void renderColoredQuad(tgt::ivec2 ll, tgt::ivec2 ur, tgt::vec4 color, bool useMenuWidgetOrientationInvert = true) const;

    /**
     * Renders the message using the font object of the mode widget
     */
    void renderMenuWidgetMessage(tgt::ivec2 position, tgt::ivec2 menuLL, std::string message, tgt::vec4 color, tgt::Font* font);

    /**
     *Renders the geometry and returns the texture
     */
    tgt::Texture* renderGeometryTexture(TriangleMeshGeometrySingleColor* geometry, tgt::Shader* shader, tgt::Camera* camera); //bp3d

    /**
     *
     */
    void renderSphere(float radius, float shininess, const tgt::vec4& lightPosition, const tgt::vec4& diffuseLight, const tgt::vec4& specularLight, const tgt::vec4& ambientLight, const tgt::ivec2& position);
protected:
    virtual void setDescriptions() {
        setDescription("Doing cool stuff. All the time.");
    }

    virtual void beforeProcess();

    virtual void process();

    virtual bool isWidgetInMenu(tgt::vec2 pos, float radius) const;

    void renderInputImage();

    //renders menu buttons in lower left and upper right corner
    void renderMenuButtons();

    //renders widgets not in menu
    void renderWidgetsOnTable();

    /**
     * renders widget
     *
     * @param widget widget to be rendered
     */
    void renderWidget(TouchTableWidget* widget);

    //renderes widgets menu if menu is open
    void renderWidgetsInMenu();

    //renders menu
    void renderMenu();


    /**
     * Update the vertices of the menu for the current list of widgets that are included in the menu.
     *
     * @param buttonPos center position of the button that triggered the menu
     */
    void updateMenuCoordinates(tgt::vec2 buttonPos);

    virtual void initialize();
    virtual void deinitialize();

    /**
     * handles touchevent (e.g. got menu or widget hit?, if touchevent is associated with widget pass it on,...)
     *
     * @param e touchevent to handle
     */
    virtual void handleTouchEvent(tgt::TouchEvent* e);

    /**
     * if widget collides with other widgets finds position to place widget avoiding collision
     *
     * @param destPos position to which widget is supposed to move to
     * @param widget widget that is moved
     * @param collisionWidgetVector widgets to check collision with
     */
    virtual void handleWidgetCollision(tgt::vec2 destPos, TouchTableWidget* widget, std::vector<TouchTableWidget*> collisionWidgetVector);


private:

    //When in ExclusiveMode all TouchPoints are handled by currentRenderer
    BoolProperty isInExclusiveMode_;

    RenderPort privateRenderPort_; //bp3d

    RenderPort inport_;
    RenderPort outport_;

    IntProperty widgetRadius_;          ///< property for radius of the touch widgets
    IntProperty menuRadius_;            ///< property for radius of the menu buttons
    IntProperty movementThreshold_;     ///< property for determining the threshold to determine if the motion of a touch point is interpreted as motion

    BoolProperty menuOpen_;             ///< property that is set to true while the menu is open
    BoolProperty moveView_;                ///< property that is set to true while view is moved

    FloatProperty shiftDiameter_;       ///< when using 3-finger shift this controls the maximum diameter of the bounding box of the 3 touch points for using them as a single camera shifting event

    GenericCoProcessorPort<TouchTableWidget> widgetPort_;

    std::vector<TouchTableWidget*> widgetsOnTable_; ///< list of widgets that are currently on the table
    std::vector<TouchTableWidget*> widgetsInMenu_; ///< list of widgets that are currently placed in menu

    /*
     * shader and textures for rendering widgets, control elements, sliders
     */
    tgt::Shader* controlElementShader_;
    tgt::Shader* quadShader_;
    tgt::Shader* quadShader1D_;

    tgt::Texture* puckOnTex_;
    tgt::Texture* puckOffTex_;
    tgt::Texture* innerTex_;
    tgt::Texture* alternateControlTex_;
    tgt::Texture* menuButtonTex_;

    tgt::Texture* sliderMidTex_;
    tgt::Texture* sliderIndicatorTex_;
    tgt::Texture* sliderLeftEndTex_;
    tgt::Texture* sliderRightEndTex_;

    tgt::Texture* checkerBoardTex_;

    std::map<int, TouchTableWidget*> touchPointMapping_;    ///< touch point IDs that have been associated with a connected widget (e.g. for moving)
    std::vector<int> menuWidgetTouchPoints_;                ///< touch point IDs that have been associated with the current mode widget menu
    std::vector<int> shiftTouchPoints_;                     ///< touch point IDs that have been associated with the 3-finger shift gesture

    static const std::string loggerCat_;

    TouchTableOverlayMenu overlayMenu_;                     ///< menu the widgets can be placed in

    TouchTableMenuWidget* currentMenuWidget_;               ///< 0 by default, pointer to the mode widget that is currently active (eg. clipping, tf)
};

} // namespace

#endif
