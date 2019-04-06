/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_SIMILARITYPLOT_H
#define VRN_SIMILARITYPLOT_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "modules/plotting/ports/plotport.h"

#include "../ports/ensembledatasetport.h"
#include "../ports/similaritymatrixport.h"
#include "../properties/stringlistproperty.h"

namespace voreen {

class CameraInteractionHandler;
class PlotLibrary;

class VRN_CORE_API SimilarityPlot : public RenderProcessor {

    class MDSData : public Serializable {
    public:
        /// Actual principal components, to be drawn.
        std::vector<std::vector<float>> nVectors_;

        /// Corresponding eigen values.
        std::vector<float> eigenvalues_;

        /// Serialization
        virtual void serialize(Serializer& s) const;
        /// Deserialization
        virtual void deserialize(Deserializer& s);
    };

public:
    SimilarityPlot();
    virtual ~SimilarityPlot();

    virtual Processor* create() const;
    virtual std::string getClassName() const        { return "SimilarityPlot"; }
    virtual std::string getCategory() const         { return "Plotting";              }
    virtual CodeState getCodeState() const          { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const   { return true;                    }

protected:
    virtual void initialize();
    virtual void deinitialize();
    virtual void process();
    virtual bool isReady() const;
    virtual void onEvent(tgt::Event* e);
    
    void calculate();
    void adjustPropertiesToInput();

    enum ColorCoding {
        COLOR_RUN,
        COLOR_TIMESTEP,
        COLOR_RUN_AND_TIMESTEP,
        COLOR_DURATION,
    };

    ButtonProperty calculateButton_;
    ProgressProperty progressBar_;
    IntProperty numIterations_;
    IntProperty numEigenvalues_;

    IntProperty numDimensions_;
    IntProperty principalComponent_;
    BoolProperty scaleToMagnitude_;
    FloatProperty sphereRadius_;
    IntProperty fontSize_;
    BoolProperty toggleAxes_; //< used for merging plots
    OptionProperty<ColorCoding> colorCoding_;
    OptionProperty<std::string> renderedChannel_;
    StringListProperty renderedRuns_;

    FloatIntervalProperty selectedTimeSteps_;
    StringListProperty selectedRuns_;

    FileDialogProperty saveFileDialog_;
    FileDialogProperty loadFileDialog_;

    CameraProperty camera_;
    CameraInteractionHandler* cameraHandler_;
    std::unique_ptr<PlotLibrary> plotLib_;
    /*
    tgt::ivec2 grabAnchorPosition_;
    bool interaction_;
    bool clicked_;
    */

    /// Inport for the ensemble data structure.
    EnsembleDatasetPort ensembleInport_;

    /// Inport for similarity matrices.
    SimilarityMatrixPort similarityMatrixInport_;

    /// The whole output.
    RenderPort outport_;

    /// Actual plot.
    RenderPort privatePort_;

    /// Used for picking in the plot.
    RenderPort pickingBuffer_;

    /// Plotport used to output eigenvalues.
    PlotPort eigenValueOutport_;

    /// The hash of the ensemble the plot was generated for.
    std::string ensembleHash_;

    /// Actual MDS data for each channel.
    std::vector<MDSData> mdsData_;

    /// Sphere geometry for timestep selection.
    GlMeshGeometryUInt16Simple sphere_;

    /// Rendering order.
    std::deque<int> renderingOrder_;

private:

    void mouseClickEvent(tgt::MouseEvent* e);
    void renderingPass(bool picking);
    void renderAxes();
    void drawTimeStepSelection(size_t runIdx, size_t timeStepIdx, const tgt::vec3& position) const;
    tgt::vec3 getColor(size_t runIdx, size_t timeStepIdx, bool picking) const;
    MDSData computeFromDM(const SimilarityMatrix& matrix, float epsilon = -1.0f);

    void outputEigenValues();
    void renderedChannelsChanged();
    void save();
    void load();

    static const std::string loggerCat_;
};

} // namespace

#endif
