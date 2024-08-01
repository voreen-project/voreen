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

#ifndef VRN_TABBEDVIEW_H
#define VRN_TABBEDVIEW_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

/**
 * Allows you to switch between multiple views using a tab-bar.
 */
class TabbedView : public RenderProcessor {
public:
    TabbedView();
    ~TabbedView();

    virtual bool isReady() const;
    virtual void process();

    virtual std::string getCategory() const { return "View"; }
    virtual std::string getClassName() const { return "TabbedView"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }
    virtual Processor* create() const;

    virtual void initialize();

    virtual void invalidate(int inv = INVALID_RESULT);

    virtual void onEvent(tgt::Event* e);
protected:
    virtual void setDescriptions() {
        setDescription("Allows you to switch between multiple views using a tab-bar.");
    }

    int getBorderWidth();
    tgt::ivec2 getInternalSize();

    void toggleMaximization(tgt::MouseEvent* me);

    void handleMouseEvent(tgt::MouseEvent* e);

    //--------------
    //  Callbacks
    //--------------
    /** Updates the inport render sizes. */
    void updateSizes();
    /** Updates the options of the currentView_ property. */
    void updateNumOptions();

    //--------------
    //  Members
    //--------------
        //ports
    RenderPort inport1_;
    RenderPort inport2_;
    RenderPort inport3_;
    RenderPort inport4_;
    RenderPort outport_;
        //properties
    IntOptionProperty currentView_;
            //tab-bar settings
    BoolProperty hideTabbar_;
    BoolProperty renderAtBottom_;
    ColorProperty borderColor_;
    ColorProperty buttonColor_;
    ColorProperty buttonHoverColor_;
    ColorProperty buttonActiveColor_;
    ColorProperty textColor_;
    ColorProperty textHoverColor_;
    ColorProperty textActiveColor_;
    FontProperty fontProp_;
            //tab configuration
    StringProperty tabText1_;
    StringProperty tabText2_;
    BoolProperty enableTab3_;
    StringProperty tabText3_;
    BoolProperty enableTab4_;
    StringProperty tabText4_;
        // event helper
    bool insideViewPort_;
    int mouseOverButton_;  // counts the tabs from left to right NOT the ports
    bool isDragging_;
};

} // namespace voreen

#endif
