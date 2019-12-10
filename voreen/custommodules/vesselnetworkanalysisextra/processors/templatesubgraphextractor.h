/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_TEMPLATE_SUBGRAPH_EXTRACTOR_H
#define VRN_TEMPLATE_SUBGRAPH_EXTRACTOR_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/properties/lightsourceproperty.h"

namespace voreen {

struct TemplateBranch;
struct BranchSearchCone;

class TemplateSubgraphExtractor : public GeometryRendererBase {
public:
    TemplateSubgraphExtractor();
    virtual ~TemplateSubgraphExtractor();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "TemplateSubgraphExtractor"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new TemplateSubgraphExtractor(); }

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used extract a subgraph starting from the node closest to a specified starting point and following a specified graph template.");
    }

    enum ConeRenderMode {
        CONES_NONE,
        CONES_MISSING,
        CONES_ALL
    };

    enum SolvingStrategy {
        STRATEGY_GREEDY,
        STRATEGY_GLOBAL
    };

    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

    //Render missing branches
    virtual void render();
    bool usesTransparency() const {
        return true;
    }
    void renderTransparent();
    void renderCones(bool useTransparency);
    void renderBranchSearchCone(const BranchSearchCone& branch, const tgt::vec4& color, bool useTransparency);
    void adjustColorPropertyVisibility();

    VesselGraphPort inport_;
    GeometryPort startingPoint_;
    VesselGraphPort outport_;

    // properties
    BoolProperty enabled_;
    FileDialogProperty templateFile_;
    ButtonProperty reloadTemplate_;
    BoolProperty keepBounds_;
    OptionProperty<SolvingStrategy> solvingStrategy_;
    OptionProperty<ConeRenderMode> coneRenderMode_;
    ColorProperty missingBranchConeColor_;
    ColorProperty locatedBranchConeColor_;
    FloatProperty coneLength_; ///<in mm

    std::vector<BranchSearchCone> missingBranches_;
    std::vector<BranchSearchCone> locatedBranches_;

    LightSourceProperty lightPosition_; ///< The position of the light source in world coordinates
    ColorProperty lightAmbient_;        ///< The light source's ambient color
    ColorProperty lightDiffuse_;        ///< The light source's diffuse color
    ColorProperty lightSpecular_;       ///< The light source's specular color
    FloatProperty materialShininess_;   ///< The material's specular exponent


    static const std::string loggerCat_;

private:
    std::unique_ptr<VesselGraph> extractSubgraph(const VesselGraph& input, const tgt::vec3& starting_point, const TemplateBranch& tpl, bool keepBounds);
    std::unique_ptr<VesselGraph> extractSubgraphGlobal(const VesselGraph& input, const tgt::vec3& starting_point, const TemplateBranch& tpl, bool keepBounds);

};

} // namespace voreen
#endif // VRN_TEMPLATE_SUBGRAPH_EXTRACTOR_H
