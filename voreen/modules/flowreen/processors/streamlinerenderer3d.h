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

#ifndef VRN_STREAMLINERENDERER3D_H
#define VRN_STREAMLINERENDERER3D_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "modules/flowreen/ports/streamlinelistport.h"
#include "voreen/core/ports/renderport.h"

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "modules/flowreen/utils/colorcodingability.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/matrixproperty.h"

namespace voreen {

    class VolumeBase;

/**
 * Performs rendering of streamlines from a stationary input flow volume by
 * using geometric primitives like lines, tubes or arrows.
 */
class StreamlineRenderer3D : public RenderProcessor {

public:
    StreamlineRenderer3D();
    virtual ~StreamlineRenderer3D();

    virtual Processor* create() const { return new StreamlineRenderer3D(); }

    virtual std::string getCategory() const { return "Flow Visualization"; }
    virtual std::string getClassName() const { return "StreamlineRenderer3D"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }
protected:
    virtual void setDescriptions() {
        setDescription("Renders streamline objects from 3D vector fields for visualizing flow. " \
                       "As input a StreamlineList is used. See <b>StreamlineCreator</b>.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    //------------------
    // Enums
    //------------------
    enum RenderingOption {
        OPTION_STREAMLINES,
        OPTION_STREAMLINEBUNDLES_NOISE_OFF,
        OPTION_STREAMLINEBUNDLES_NOISE_ON,
    };
    friend class OptionProperty<RenderingOption>;

    enum StreamlineBundleStyle {
        BUNDLE_STYLE_LINES,
        BUNDLE_STYLE_TUBES,
        BUNDLE_STYLE_ARROWS,
    };
    friend class OptionProperty<StreamlineBundleStyle>;

    enum StreamlineStyle {
        STYLE_LINES,
        STYLE_ARROWS,
    };
    friend class OptionProperty<StreamlineStyle>;

    enum StreamlineColorCoding {
        COLOR_VELOCITY,
        COLOR_DIRECTION,
    };
    friend class OptionProperty<StreamlineColorCoding>;

    friend class PBLinkControl;

    //------------------
    // Callbacks
    //------------------
    void onStreamlineDataChange();   ///< called from inport_
    void onStyleChange();            ///< called from streamlineStyleProp_
    void onColorChange();            ///< called from colorProp_
    void computeDirectionColorRotationMatrix(); ///< called from rotateAroundXYZ_

    //------------------
    // Helper
    //------------------
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    void compile();
    void rebuild();
    void buildStreamlineData(const std::vector<Streamline>& streamlines = std::vector<Streamline>());
    void buildStreamlineBundleData(const std::vector<StreamlineBundle>& bundles = std::vector<StreamlineBundle>());
    GlMeshGeometryUInt32Color* createStreamTube(const StreamlineBundle& bundle) const;
    GlMeshGeometryUInt32Color* createArrowPath(const StreamlineBundle& bundle) const;
    tgt::mat4 createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity) const;

    //------------------
    //  Members
    //------------------
        //ports
    StreamlineListPort streamlineInport_;
    RenderPort imgOutport_;
        //style
    OptionProperty<StreamlineStyle> streamlineStyleProp_;                ///< used to change streamline representations
    OptionProperty<StreamlineBundleStyle> streamlineBundleStyleProp_;    ///< used to change streamline bundle representations
        //options
    OptionProperty<RenderingOption> renderingOptionProp_;                ///< used to determine rendering output
        //color
    OptionProperty<StreamlineColorCoding> colorProp_;   ///< color encoding
    TransFunc1DKeysProperty tfProp_;                    ///< tf for velocity color coding
        VolumeBase* tfVolume_;                          ///< volume for tf fit to data
    FloatMat4Property colorRotationMatrix_;             ///< rotation matrix for the direction encoding color
    FloatOptionProperty rotateAroundX_;                 ///< rotate around x to get nice color coding
    FloatOptionProperty rotateAroundY_;                 ///< rotate around y to get nice color coding
    FloatOptionProperty rotateAroundZ_;                 ///< rotate around z to get nice color coding

        //must haves
    ShaderProperty streamlineShaderProp_;               ///< used for rendering
        bool requiresRecompileShader_;
    CameraProperty cameraProp_;                         ///< the camera
        CameraInteractionHandler* cameraHandler_;       ///< handler to allow navigation

        //rendering
    std::vector<GlMeshGeometryUInt32Color* > meshes_;   ///< all other meshes to be rendered
    size_t bundleMeshStartIndex_;                       ///< determines the position of the first bundle mesh
    bool requiresRebuild_;
};

}   // namespace

#endif  // VRN_STREAMLINERENDERER3D_H
