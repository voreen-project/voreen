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

#ifndef VRN_STAR_COORDINATES_H
#define VRN_STAR_COORDINATES_H

#include <Eigen/Eigenvalues>

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

    class VRN_CORE_API StarCoordinates : public RenderProcessor {

        class StarCoordinateAxis : public Serializable {
        public:
            tgt::vec3 color_;
            tgt::vec2 direction_;
            float length_;
            std::string label_;
            size_t id_;
            bool enabled_;

            StarCoordinateAxis();
            StarCoordinateAxis(std::string label, tgt::vec3 color, tgt::vec2 direction, float length, size_t id);

            /// Serialization
            virtual void serialize(Serializer& s) const;
            /// Deserialization
            virtual void deserialize(Deserializer& s);
        };

        class StarCoordinatePoint : public Serializable {
        public:
            size_t seriesID_;
            size_t timeStepID_;

            StarCoordinatePoint();
            StarCoordinatePoint(size_t seriesID, size_t timeStepID);

            /// Serialization
            virtual void serialize(Serializer& s) const;
            /// Deserialization
            virtual void deserialize(Deserializer& s);
        };

    public:
        StarCoordinates();
        virtual ~StarCoordinates();

        virtual Processor* create() const;
        virtual std::string getClassName() const { return "StarCoordinates"; }
        virtual std::string getCategory() const { return "Plotting"; }
        virtual CodeState getCodeState() const { return CODE_STATE_TESTING; }
        virtual bool usesExpensiveComputation() const { return true; }

        /// Serialization
        virtual void serialize(Serializer& s) const;
        /// Deserialization
        virtual void deserialize(Deserializer& s);

    protected:
        virtual void initialize();
        virtual void deinitialize();
        virtual void process();
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

        void mouseEvent(tgt::MouseEvent* e);
        void wheelEvent(tgt::MouseEvent* e);

        void renderingPass(bool picking);
        void setPlotStatus();
        void renderAxes();
        void renderSelection();
        void renderProjection(bool picking);
        void renderTooltip() const;

        // utility
        tgt::vec3 getColor(size_t seriesIdx, size_t timeStepIdx, bool picking) const;
        tgt::vec3 getColorForValue(float val) const;
        tgt::vec2 getRange(size_t fieldIdx);

        // input processing
        void adjustToTimeserieslist();
        void createProjection();
        void extractPoints(const TimeSeriesList* timeseries);
        void projectPoints();
        void resetAxes();
        void resetSelection();

        // interaction
        float getZoomFactor(const float& acuteness, const bool& zoomIn) const;

        // utility functions
        bool almostEqual(float a, float b, float eps) { return abs(a - b) < eps; }

        CameraProperty camera_;
        CameraInteractionHandler* cameraHandler_;
        std::unique_ptr<PlotLibrary> plotLib_;

        TimeSeriesListPort timeseriesInport_;         /// Inport for the time series data structure.

        RenderPort outport_;         /// The whole output.
        RenderPort privatePort_;     /// Actual plot.
        RenderPort pickingBuffer_;   /// Used for picking in the plot.

        // gui
        ButtonProperty calculateButton_;
        StringListProperty renderedMembers_;
        StringListProperty selectedMembers_;
        OptionProperty<ColorCoding> colorCoding_;
        FloatProperty labelOffset_;
        FloatIntervalProperty angleRange_;
        FloatIntervalProperty tempAnomalyRange_;
        TransFunc1DKeysProperty transferFunc_;
        ButtonProperty resetAxesButton_;
        ButtonProperty resetSelectionButton_;
        BoolProperty showTooltip_;

        // Linking
        IntProperty volumeDimension_;
        IntVec3Property position_;

        tgt::Shader* bypassProgram_;

        bool updateAxes;

        /// Selected members (sorted for faster access).
        std::set<int> subSelection_;
        int selectedAxis_;
        std::vector<StarCoordinateAxis> starCoordinateAxes_;
        std::vector<StarCoordinatePoint> starCoordinatePoints_;

        // projection
        Eigen::MatrixXf projection_;
        Eigen::MatrixXf nDpoints_;
        Eigen::MatrixXf projectedPoints_;
        bool updatedProjection;

        // user interaction
        float zoom_;
        float mouseWheelZoomAcuteness_;

        /// sphere used for selection
        FloatProperty sphereRadius_;
        GlMeshGeometryUInt16Simple sphere_;
        FloatProperty fontScaling_;
        FloatProperty pointSize_;

        /// Last picked member and time step (also when hovering).
        struct Hit {
            int x, y;
            int memberIdx;
            int timeStepIdx;
        };
        boost::optional<Hit> lastHit_;

        /// The font being used for rendering the tooltip.
        static const std::string fontName_;
        IntProperty fontSize_;

        static const std::string loggerCat_;
    };

} // namespace

#endif
