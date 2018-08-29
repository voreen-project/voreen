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

#ifndef VRN_VESSELGRAPHCOMPARISON_H
#define VRN_VESSELGRAPHCOMPARISON_H

#include "voreen/core/processors/processor.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/processors/geometryrendererbase.h"

#include "modules/plotting/ports/plotport.h"

namespace voreen {

template<class T>
struct Matching {
    std::vector<std::pair<const T*, const T*>> matches_;
    std::vector<const T*> non_matched1_;
    std::vector<const T*> non_matched2_;

    float matchRatio() const;
    bool hasContent() const;
};

class VesselGraphComparison : public GeometryRendererBase {
public:
    VesselGraphComparison();
    virtual ~VesselGraphComparison();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "VesselGraphComparison"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VesselGraphComparison(); }
    virtual bool isReady() const;

    //Render matches
    virtual void render();

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used to compare two vessel graphs and compute similarity measures.");
    }

    virtual void process();

    void compare(const VesselGraph& g1, const VesselGraph& g2);

    Matching<VesselGraphNode> matchNodesMutualNN(const VesselGraph& g1, const VesselGraph& g2) const;
    Matching<VesselGraphEdge> matchEdgesViaNodes(const VesselGraph& g1, const VesselGraph& g2, const Matching<VesselGraphNode>& node_matching) const;

    Matching<VesselGraphEdge> matchEdgesViaBipartiteGraphMatching(const VesselGraph& g1, const VesselGraph& g2) const;

    template<class D>
    Matching<VesselGraphEdge> matchEdgesViaHungarianAlgorithm(const VesselGraph& g1, const VesselGraph& g2, D distance) const;

    template<class D>
    Matching<VesselGraphEdge> matchEdgesLAP(const VesselGraph& g1, const VesselGraph& g2, D distance) const;

    enum MatchingAlgorithm {
        MUTUAL_NN,
        HUNGARIAN_NODES,
        HUNGARIAN_EDGES,
        HUNGARIAN_EDGES_QUANTIL_THRESHOLD,
        LAP
    };

    enum MatchRenderMode {
        NODES,
        EDGES,
        NONE,
    };

    // ports
    VesselGraphPort inport1_;
    VesselGraphPort inport2_;
    //PlotPort plotOutport_;

    // properties
    BoolProperty enabled_;
    OptionProperty<MatchingAlgorithm> matchingAlgorithm_;
    OptionProperty<MatchRenderMode> renderMode_;
    FileDialogProperty statExportFile_;
    StringProperty datasetIdentifier_;
    FloatProperty nodeMatchRatio_;
    FloatProperty edgeMatchRatio_;
    FloatProperty lengthSimilarity_;
    FloatProperty nodeMatchingCost_;

    // (Non-)Match Rendering
    FloatProperty crossRadius_;
    std::unique_ptr<Matching<VesselGraphEdge>> lastEdgeMatching_;
    std::unique_ptr<Matching<VesselGraphNode>> lastNodeMatching_;

    static const std::string loggerCat_;
};

} // namespace voreen
#endif // VRN_VESSELGRAPHCOMPARISON_H
