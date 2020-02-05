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

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

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

    enum StreamlineStyle {
        STYLE_LINES,
        STYLE_TUBES,
        STYLE_ARROWS,
    };

    enum StreamlineColorCoding {
        COLOR_VELOCITY,
        COLOR_DIRECTION,
    };

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
    GlMeshGeometryUInt32Color* createLineGeometry(const std::vector<Streamline>& streamlines);
    GlMeshGeometryUInt32Color* createTubeGeometry(const Streamline& streamline) const;
    GlMeshGeometryUInt32Color* createArrowGeometry(const Streamline& streamline) const;

    //------------------
    //  Members
    //------------------
        //ports
    StreamlineListPort streamlineInport_;
    RenderPort imgOutport_;
        //style
    OptionProperty<StreamlineStyle> streamlineStyle_;   ///< used to change streamline representations
        //color
    OptionProperty<StreamlineColorCoding> color_;       ///< color encoding
    TransFunc1DKeysProperty transferFunction_;          ///< tf for velocity color coding
    std::unique_ptr<VolumeBase> tfVolume_;              ///< volume for tf fit to data
    FloatMat4Property colorRotationMatrix_;             ///< rotation matrix for the direction encoding color
    FloatOptionProperty rotateAroundX_;                 ///< rotate around x to get nice color coding
    FloatOptionProperty rotateAroundY_;                 ///< rotate around y to get nice color coding
    FloatOptionProperty rotateAroundZ_;                 ///< rotate around z to get nice color coding
    BoolProperty enableShading_;                        ///< enables phong shading
    BoolProperty enableMaximumIntensityProjection_;     ///< enables maximum intensity projection (MIP)

        //must haves
    ShaderProperty streamlineShader_;                   ///< used for rendering
    bool requiresRecompileShader_;
    CameraProperty camera_;                             ///< the camera
    CameraInteractionHandler* cameraHandler_;           ///< handler to allow navigation

        //rendering
    std::vector<std::unique_ptr<GlMeshGeometryUInt32Color>> meshes_;   ///< all other meshes to be rendered
    bool requiresRebuild_;

    static std::string loggerCat_;
};

}   // namespace

#endif  // VRN_STREAMLINERENDERER3D_H
