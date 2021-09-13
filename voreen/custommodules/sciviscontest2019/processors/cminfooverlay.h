/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#ifndef VRN_CMINFOOVERLAY_H
#define VRN_CMINFOOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include "tgt/texture.h"


/**
 * An overlay processor that displays a provided (informative) image when the user
 * clicks the otherwise overlayed info icon (i).
 */
namespace voreen {

    class TriangleMeshGeometryBase;
    class TriangleMeshGeometryUInt16IndexedColorNormal;
class VRN_CORE_API CMInfoOverlay : public ImageProcessor {
public:
    CMInfoOverlay();
    ~CMInfoOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "CMInfoOverlay"; }
    virtual std::string getCategory() const     { return "Image Processing"; }
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Adds a simple icon which renders an info image when clicked.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void onEvent(tgt::Event *ev) override;

private:
    /**
     * Tries to load the info image and outputs a warning if the image cannot be loaded
     */
    void tryLoadInfoImageTex();

    /**
     * Tries to load the icon and outputs a warning if the icon cannot be loaded
     */
    void tryLoadIconTex();

    /**
     * Gets the bounds of the icon in opengl image coordinates [-1,1]^2
     * @return 2d-bounds in vec4: [ll.x, ll.y, ur.x, ur.y]
     */
    tgt::vec4 getIconBounds();

    /**
     * Possible positions of the icon
     */
    enum POSITION {
        TOP_LEFT,
        TOP_RIGHT,
        BOTTOM_LEFT,
        BOTTOM_RIGHT
    };
    //ports
    /// Image the icon or the image will be drawn on
    RenderPort inport_;
    /// Output for input image overlayed with icon or info image.
    RenderPort outport_;

    //basic
    /// Used to disable this overlay
    BoolProperty enableProp_;
    //image
    /// Dialog to choose an info image for this overlay
    FileDialogProperty infoImageFileProp_;
    //position
    /// Property to choose the position of the icon
    OptionProperty<POSITION> positionProp_;

    /// Size of the icon in opengl window coordinates
    FloatProperty sizeProp_;
    /// Alpha of the icon
    FloatProperty alphaProp_;

    /// Whether or not the info image is shown (true: info image, false: only icon)
    BoolProperty activatedProp_;

    /// Texture of the icon
    tgt::Texture* iconTexture_;
    /// Texture of the info image
    tgt::Texture* infoTexture_;

    /// Used in logging
    static const std::string loggerCat_;
};

} // namespace

#endif
