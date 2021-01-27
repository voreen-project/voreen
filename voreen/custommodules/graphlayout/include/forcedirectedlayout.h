/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_FORCEDIRECTEDLAYOUT_H
#define VRN_FORCEDIRECTEDLAYOUT_H

#include "forcedirectednodegraph.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "modules/plotting/ports/plotport.h"

namespace voreen {

class PlotBase;
class PlotData;

class ForceDirectedLayout : public Processor {
public:
    ForceDirectedLayout();

    /// default destructor
    virtual ~ForceDirectedLayout();

    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Plotting"; }
    virtual std::string getClassName() const { return "ForceDirectedLayout"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

    //virtual bool usesExpensiveComputation() const { return true; }
    virtual bool isEndProcessor() const;
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Applies force-directed layout algorithms on graphs.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:
    void readInports();

    void calculateForces();
    void moveNodes();

    void initLayout();
    void applyForces();
    void applyForcesStep();

    void writeOutports();


    PlotPort nodesIn_;
    PlotPort nodesOut_;
    PlotPort edgesIn_;
    PlotPort edgesOut_;

    ForceDirectedNodeGraph nodeGraph_;

    FloatProperty springForce_;
    FloatProperty springLength_;
    FloatProperty charge_;

    FloatProperty horizontalMagneticForce_;
    FloatProperty verticalMagneticForce_;
    FloatProperty magneticAlpha_;
    FloatProperty magneticBeta_;
    FloatProperty edgeAngleEqualization_;

    FloatProperty damping_;
    FloatProperty threshold_;

    ButtonProperty initLayout_;
    ButtonProperty applyForcesStep_;
    ButtonProperty applyForces_;

    static const std::string loggerCat_;

};

} // namespace voreen

#endif // VRN_FORCEDIRECTEDLAYOUT_H
