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

#ifndef VRN_TIMESERIESPLOT_H
#define VRN_TIMESERIESPLOT_H

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
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "modules/plotting/ports/plotport.h"

#include "../ports/timeserieslistport.h"
#include "modules/ensembleanalysis/ports/similaritymatrixport.h"

namespace voreen {

class CameraInteractionHandler;
class PlotLibrary;

class VRN_CORE_API TimeseriesPlot : public RenderProcessor {

    class Embedding : public Serializable {
    public:
        /// Actual principal components, to be drawn.
        // * First layer encodes member.
        // * Second layer encodes time step.
        // * Third layer encodes principle component.
        std::map<int, std::vector<std::vector<float>>> nVectors_;

        /// Corresponding eigen values.
        std::vector<float> eigenvalues_;

        /// (Optional) member names.
        std::vector<std::string> names_;

        /// Serialization
        virtual void serialize(Serializer& s) const;
        /// Deserialization
        virtual void deserialize(Deserializer& s);
    };

public:
    TimeseriesPlot();
    virtual ~TimeseriesPlot();

    virtual Processor* create() const;
    virtual std::string getClassName() const        { return "TimeseriesPlot";   }
    virtual std::string getCategory() const         { return "Plotting";         }
    virtual CodeState getCodeState() const          { return CODE_STATE_TESTING; }
    virtual bool usesExpensiveComputation() const   { return true;               }

protected:
    virtual void initialize();
    virtual void deinitialize();
    virtual void process();
    virtual bool isReady() const;
    virtual void onEvent(tgt::Event* e);

    virtual void setDescriptions() {
        setDescription("This processor uses a classical Multi Dimensional Scaling (MDS) appraoch in order to "
                       "create a distance-presevering, low-dimensional embedding for the given time series. "
                       "Distances are stored in the input <br>SimilarityMatrix</br>, e.g. provided by "
                       "<br>SimilarityMatrixCreator</br>. Members with multiple runs are represented by "
                       "curves, members with a single time step by points (in a 2D and 3D embedding) or dotted "
                       "lines (when using a 1D embedding). Use CTRL to add hovered members to the subselection and "
                       "use ALT to remove them. Click the middle mouse button to reset the subselection to all "
                       "currently rendered members. Pressing shift and the middle mouse button will change the "
                       "rendering order of hovered members.");

        timeseriesInport_.setDescription("The input time series list.");
        //similarityMatrixInport_.setDescription("Input similarity matrix. It must have been created for the "
        //                                       "currently connected time series");
        eigenValueOutport_.setDescription("Outputs the eigenvalues of the current embedding. Connect a <br>BarPlot</br> "
                                          "Processor to display them.");

        numIterations_.setDescription("Number of iterations that define accuracy of the embedding, so less steps "
                                      " will calculate faster. The algorithm will, however, terminate early "
                                      "if the result converged");
        numEigenvalues_.setDescription("Number of eigenvalues to be calculated. Since only the largest three "
                                       "principal directions (to which the eigenvalues correspond to) can be displayed, "
                                       "it might be sufficient to only calculate three. However, the significant "
                                       "intrinsic dimensionality might be larger than three!");
        numDimensions_.setDescription("Number of dimensions to be used to visualize the embedding");
        principleComponent_.setDescription("If the 1D Embedding is used, any principle component can be plotted over time");
//        scaleToMagnitude_.setDescription("If the 3D Embeding is used, all axis can be scaled such that they represent "
//                                         "the variance of the data represented by the principle component with respect "
//                                         "to the first and largest principle component");
        sphereRadius_.setDescription("Currently selected time steps and members with only a single time steps both "
                                     "are represented by spheres in the 2D and 3D embedding. This sets their radius.");
        fontSize_.setDescription("Sets the font size");
        showTooltip_.setDescription("Enables/Disables tool tips when hovering over members");
        renderTimeSelection_.setDescription("Enables/Disables indication of the currently selected time step/range");
        colorCoding_.setDescription("Defines the color coding for the time series to be used");
        renderedMembers_.setDescription("Can be used to show all members at the same time or a subset thereof");

        firstSelectedMember_.setDescription("Use left-click to perform selection of a member.");
        firstSelectedTimeInterval_.setDescription("Use left-click to perform selection of a time step.");

        secondSelectedMember_.setDescription("Use right-click to perform selection of another member. ");
        secondSelectedTimeInterval_.setDescription("Use right-click to perform selection of another time step.");
    }

private:

    enum ColorCoding {
        COLOR_MEMBER,
        COLOR_TIMESTEP,
        COLOR_MEMBER_AND_TIMESTEP,
        COLOR_DURATION,
        COLOR_FIELD,
        COLOR_THRESHOLD,
        COLOR_SLAB,
        COLOR_PLUME,
        COLOR_SLAB_AND_PLUME
    };

    void adjustToTimeserieslist();
    void mouseEvent(tgt::MouseEvent* e);

    void renderingPass(bool picking);
    void renderEmbedding1D(bool picking);
    void renderEmbedding2D(bool picking);
    void renderEmbedding3D(bool picking);
    void renderAxes();
    void renderTooltip() const;
    void renderTimeStepSelection(size_t seriesIdx, size_t timeStepIdx, const tgt::vec3& position, const tgt::vec3& color) const;
    tgt::vec3 getColor(size_t seriesIdx, size_t timeStepIdx, bool picking) const;
    tgt::vec2 getRange(size_t fieldIdx);

    void createEmbeddings();
    Embedding createEmbedding(const TimeSeriesList* timeseries, ProgressReporter& progressReporter, float epsilon = 0.0f) const;

    void outputEigenValues();
    void renderedMembersChanged();
    void saveEmbeddings();
    void loadEmbeddings();
    tgt::vec3 getColorForValue(float val) const;

    ButtonProperty calculateButton_;
    BoolProperty autoCalculate_;
    ProgressProperty progressBar_;
    IntProperty numIterations_;
    IntProperty numEigenvalues_;

    IntProperty numDimensions_;
    IntProperty principleComponent_;
    //BoolProperty scaleToMagnitude_;
    FloatProperty sphereRadius_;
    IntProperty fontSize_;
    BoolProperty showTooltip_;
    BoolProperty renderTimeSelection_;
    OptionProperty<ColorCoding> colorCoding_;
    FloatIntervalProperty angleRange_;
    FloatIntervalProperty tempAnomalyRange_;
    TransFunc1DKeysProperty transferFunc_;  ///< the property that controls the transfer function
    //OptionProperty<std::string> renderedField_;
    StringListProperty renderedMembers_;

    StringListProperty firstSelectedMember_;
    FloatIntervalProperty firstSelectedTimeInterval_;
    StringListProperty secondSelectedMember_;
    FloatIntervalProperty secondSelectedTimeInterval_;

    FileDialogProperty saveFileDialog_;
    ButtonProperty saveButton_;
    FileDialogProperty loadFileDialog_;
    ButtonProperty loadButton_;

    CameraProperty camera_;
    CameraInteractionHandler* cameraHandler_;
    std::unique_ptr<PlotLibrary> plotLib_;

    IntProperty volumeDimension_;
    IntVec3Property position_;


    /// Inport for the time series data structure.
    TimeSeriesListPort timeseriesInport_;

    /// Inport for similarity matrices.
    //SimilarityMatrixPort similarityMatrixInport_;

    /// The whole output.
    RenderPort outport_;

    /// Actual plot.
    RenderPort privatePort_;

    /// Used for picking in the plot.
    RenderPort pickingBuffer_;

    /// Plotport used to output eigenvalues.
    PlotPort eigenValueOutport_;

    //GeometryPort pointListOutport_;

    /// Actual MDS embedding for each field (all members selected).
    Embedding embedding_;

    /// Sphere geometry for timestep selection.
    GlMeshGeometryUInt16Simple sphere_;

    /// Rendering order.
    std::deque<int> renderingOrder_;

    /// Selected members (sorted for faster access).
    std::set<int> subSelection_;

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
};

} // namespace

#endif
