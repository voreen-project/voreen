/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_CELLMIGRATIONSCORE_H
#define VRN_CELLMIGRATIONSCORE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"

// for quantification center output
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

// plotting stuff
#include "../../../modules/plotting/datastructures/plotdata.h"
#include "../../../modules/plotting/ports/plotport.h"
#include "../../../modules/plotting/datastructures/plotcell.h"

//#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {

/**
 * Set a target point and quantify fluorescence in a 4D data set represented by a volume list.
 */
class VRN_CORE_API CellMigrationScore : public Processor {

public:
    CellMigrationScore();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "CellMigrationScore";     }
    virtual std::string getCategory() const  { return "Quantification";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    // progress is shown during quantification
    virtual bool usesExpensiveComputation() const { return true;    }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Set a target point and a radius of influence and quantify fluorescence within the area, computing a score for each time step starting at a specified time.");
        quantificationCenter_.setDescription("The location of the induced wound (in voxel space).");
        quantificationFrame_.setDescription("The time step interval where the quantification should be performed.");
        radius_.setDescription("The influence radius whre quantification should be performend (in micrometers!).");
        startComputation_.setDescription("Starts the quantification");
    }

    virtual void process();
    virtual void initialize();

    virtual void computeQuantification();

    /// check if the VolumeList represents a time series, i.e., if every volume has the same dimensions, spacing, offset, data type, voxel-to-world matrix, and has a VolumeDisk representation
    bool checkVolumeList(const VolumeList* collection) const;

    /// center of the induced wound
    IntVec3Property quantificationCenter_;

    /// radius of the quantification
    IntProperty radius_;
    
    /// time step interval in which the score computation should be performed
    IntIntervalProperty quantificationFrame_;

    /// starts the quantification for all time steps
    ButtonProperty startComputation_;

    /// progress in application mode
    ProgressProperty progressProperty_;

    /// channel of the quantification (usually not necessary for single-channel volumes)
    IntProperty channel_;

    BoolProperty useClipRegion_;
    IntBoundingBoxProperty clipRegion_;

    /// Number of subdivions for the quantification cylinder
    IntProperty cylinderSubdivisions_;

    /// plotting port for quantification results
    PlotPort quantificationPlot_;

    /// plotting port that only contains the normalized weightes quantification, i.e., the actual score
    PlotPort normalizedQuantificationPlot_;

    /// outputs the quantification center as a pointlist containing a single point with its coordinates in voxel space
    GeometryPort centerOutport_;

    /// geometry for the quantified area
    GeometryPort areaOutport_;
    //GeometryPort lineOutport_;

    /// Inport for the volume list.
    VolumeListPort inport_;

    /// The volume port the selected volume is written to.
    //VolumePort outport_;

    static const std::string loggerCat_;

private:
    void adjustToVolumeList();

    void useClipRegionChanged();

};

} // namespace

#endif
