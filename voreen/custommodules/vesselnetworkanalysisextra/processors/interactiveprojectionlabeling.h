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

#ifndef VRN_INTERACTIVEPROJECTIONLABELING_H
#define VRN_INTERACTIVEPROJECTIONLABELING_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include "tgt/vector.h"

namespace voreen {

struct LabelProjection;
struct LabelGuard {
public:
    LabelGuard(LabelProjection& labelProjection);
    ~LabelGuard();
    tgt::vec2& at(tgt::svec2);
private:
    LabelProjection& labelProjection_;
};

struct LabelProjection {
    LabelProjection(tgt::svec2 dimensions);

    const VolumeAtomic<tgt::vec2>& projection() const {
        return projection_;
    }
    LabelGuard projection_mut() {
        return LabelGuard { *this };
    }
    void bindTexture();
private:
    friend struct LabelGuard;
    void ensureTexturesPresent();

    VolumeAtomic<tgt::vec2> projection_;
    boost::optional<tgt::Texture> projectionTexture_;
};

struct ProjectionLabels {
    std::vector<std::deque<tgt::vec2>> foreground_;
    std::vector<std::deque<tgt::vec2>> background_;
    void clear();
};

struct LabelUnit : public Serializable {
    LabelUnit();
    LabelUnit(const LabelUnit&) = default;
    LabelUnit& operator=(const LabelUnit&) = default;
    LabelUnit(LabelUnit&&) = default;
    LabelUnit& operator=(LabelUnit&&) = default;

    // Configuration
    tgt::Camera camera_;
    tgt::IntBounds clippingRegion_;

    std::deque<tgt::vec2> displayLine_;
    ProjectionLabels projectionLabels_;


    // result:
    std::vector<std::vector<tgt::vec3>> backgroundLabels_;
    std::vector<std::vector<tgt::vec3>> foregroundLabels_;

    void setZoomRegion(tgt::vec2 newRegion);
    tgt::vec2 getZoomRegion() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    tgt::vec2 zoomRegion_; // subset range von [0,1]
};

class InteractiveProjectionLabeling : public RenderProcessor {
public:
    InteractiveProjectionLabeling();
    virtual ~InteractiveProjectionLabeling();

    virtual std::string getClassName() const         { return "InteractiveProjectionLabeling"; }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription( "TODO"); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();
    virtual void adjustPropertiesToInput();
    virtual void onPortEvent(tgt::Event* e, Port* port);
    virtual void initialize();
    virtual void deinitialize();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    LabelUnit& currentUnit();
    void resetCurrentUnit();
    boost::optional<size_t> currentUnitIndex();
    void startNewUnit();

    enum State {
        LABELING,
        FREE,
    };

    enum InitializationMode {
        NONE,
        BRIGHT_LUMEN,
        BRIGHT_WALL,
    };

    void updateSizes();
    void renderOverlay();
    void renderProjection();
    void updateProjection();
    void finishProjection();
    void initializeProjectionLabels();
    void synchronizeUnitIndex();

    void projectionEvent(tgt::MouseEvent* e);
    void overlayEvent(tgt::MouseEvent* e);

    boost::optional<VolumeAtomic<tgt::vec4>> getFhp();
    boost::optional<VolumeAtomic<tgt::vec4>> getLhp();

    VolumePort inport_;
    bool seedsChanged_;
    GeometryPort foregroundLabelGeometry_;
    GeometryPort backgroundLabelGeometry_;
    RenderPort overlayOutput_;
    RenderPort projectionOutput_;

    RenderPort fhp_;
    RenderPort lhp_;

    CameraProperty camera_;
    OptionProperty<InitializationMode> initializationMode_;
    TransFunc1DKeysProperty projectionTransfunc_;
    FloatProperty maxLineSimplificationDistance_;
    FloatProperty backgroundLineDistanceMultiplier_;

    static const std::string loggerCat_;

    ShaderProperty projectionShader_;

    boost::optional<LabelProjection> projection_;
    bool projectionLabelsModified_;
    bool projectionRequiresUpdate_;
    LabelUnit currentUnit_;
    boost::optional<tgt::vec2> prevProjectionMousePos_;

    State state_;

    std::vector<LabelUnit> labelUnits_;
    IntProperty currentUnitIndex_;
    IntBoundingBoxProperty clippingRegion_;
};

} // namespace voreen

#endif // VRN_INTERACTIVEPROJECTIONLABELING_H
