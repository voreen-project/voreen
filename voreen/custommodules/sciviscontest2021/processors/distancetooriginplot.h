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

#ifndef VRN_DISTANCETOORIGINPLOT_H
#define VRN_DISTANCETOORIGINPLOT_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"
#include <modules/ensembleanalysis/ports/ensembledatasetport.h>
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "modules/plotting/ports/plotport.h"

namespace voreen {

class CameraInteractionHandler;
class PlotLibrary;

class VRN_CORE_API DistanceToOriginPlot : public RenderProcessor {
public:
    DistanceToOriginPlot();
    virtual ~DistanceToOriginPlot();

    virtual Processor* create() const;
    virtual std::string getClassName() const        { return "DistanceToOriginPlot";   }
    virtual std::string getCategory() const         { return "Plotting";         }
    virtual CodeState getCodeState() const          { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const   { return true;               }

protected:
    virtual void initialize();
    virtual void deinitialize();
    virtual void process();
    virtual bool isReady() const;
    virtual void onEvent(tgt::Event* e);

private:

    enum ColorCoding {
        COLOR_MEMBER,
        COLOR_TIMESTEP,
        COLOR_MEMBER_AND_TIMESTEP,
        COLOR_DURATION,
    };

    void adjustToEnsemble();
    void mouseEvent(tgt::MouseEvent* e);

    void renderingPass(bool picking);
    void render(bool picking);
    void renderAxes();
    void renderTooltip() const;
    void renderSphere(const tgt::vec3& position, const tgt::vec3& color, bool drawBorder) const;
    void renderCircle(tgt::vec2 position, float size, tgt::vec3 color, tgt::vec2 canvasSize);
    tgt::vec3 getColor(size_t memberIdx, size_t timeStepIdx, bool picking) const;
    void calculateDistances();
    void save();
    void load();

    void renderedMembersChanged();
    void saveResult();
    void loadResult();

    tgt::ivec2 getMargins() const;

    ButtonProperty calculateButton_;
    ProgressProperty progressBar_;

    FloatProperty sphereRadius_;
    IntProperty fontSize_;
    BoolProperty showTooltip_;
    BoolProperty renderTimeSelection_;
    OptionProperty<ColorCoding> colorCoding_;

    StringListProperty firstSelectedMember_;
    FloatIntervalProperty firstSelectedTimeInterval_;
    FloatIntervalProperty firstSelectedRadiusInterval_;
    StringListProperty secondSelectedMember_;
    FloatIntervalProperty secondSelectedTimeInterval_;

    StringListProperty selectedFeatureIdProperty_;
    OptionProperty<std::string> distanceMethod_;

    FileDialogProperty saveFileDialog_;
    ButtonProperty saveButton_;
    FileDialogProperty loadFileDialog_;
    ButtonProperty loadButton_;

    CameraProperty camera_;
    CameraInteractionHandler* cameraHandler_;
    std::unique_ptr<PlotLibrary> plotLib_;

    // Color coding
    TransFunc1DKeysProperty transferFunc_;
    OptionProperty<std::string> selectField_;


    /// Inport for the ensemble data structure.
    VolumeListPort volumeListPort_;

    /// Inport for Ensmevle Data.
    EnsembleDatasetPort originalInport_;

    /// The whole output.
    RenderPort outport_;

    /// Actual plot.
    RenderPort privatePort_;

    /// Used for picking in the plot.
    RenderPort pickingBuffer_;

    /// Plotport used to output eigenvalues.
    PlotPort eigenValueOutport_;

    /// Sphere geometry for timestep selection.
    GlMeshGeometryUInt16Simple sphere_;

    /// Rendering order.
    std::deque<int> renderingOrder_;

    /// Selected members (sorted for faster access).
    std::set<int> subSelection_;

    bool isSelecting_;
    tgt::vec2 selectionStart_;
    tgt::vec2 selectionEnd_;

    /// Last picked member and time step (also when hovering).
    struct Hit {
        int x, y;
        int memberIdx;
        int timeStepIdx;
    };
    boost::optional<Hit> lastHit_;

    /// The font being used for rendering the tooltip.
    static const std::string fontName_;

    static const std::string loggerCat_;

    
    std::map<uint32_t, std::vector<tgt::ivec3>> distances_;

    // Additional properties:
    std::map<uint32_t, std::vector<tgt::vec2>> temperature_;
    std::map<uint32_t, std::vector<tgt::vec2>> vectormagnitude_;
    std::map<uint32_t, std::vector<tgt::vec2>> condactivityanamoly_;
    std::map<uint32_t, std::vector<tgt::vec2>> expansivityanamoly_;
    std::map<uint32_t, std::vector<tgt::vec2>> temperatureanomaly_;
    std::map<uint32_t, std::vector<tgt::vec2>> densityanamoly_;

    tgt::vec2 temperatureRange_;
    tgt::vec2 vectormagnitudeRange_;
    tgt::vec2 condactivityanamolyRange_;
    tgt::vec2 expansivityanamolyRange_;
    tgt::vec2 temperatureanamolyRange_;
    tgt::vec2 densityanamolyRange_;

    std::map<uint32_t, tgt::vec3> colors_;
    int largestFeatureVoxelCount_;

    tgt::vec3 getColorForValue(float val) const;
};

} // namespace

#endif
