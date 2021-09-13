/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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

#ifndef VRN_CMHALOIDTIMESTEPLINKER_H
#define VRN_CMHALOIDTIMESTEPLINKER_H


#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "../ports/cmhaloport.h"


//use namespace voreen
namespace voreen {

/**
 * Class to provide linking between time step and halo id properties e.g. for halo renderers.
 * When a new halo is selected, a path between a leaf and the root of the current mergertree
 * is selected which contains so newly selected halo. The time step is set accordingly.
 * When the time step changes, the selected halo changes along that path so that its time
 * step equals the time step property's value.
 */
class VRN_CORE_API CMHaloIDTimeStepLinker : public Processor {

public:
    /**
     * Constructor
     */
    CMHaloIDTimeStepLinker();
    ~CMHaloIDTimeStepLinker();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CMHaloIDTimeStepLinker();     }
    virtual std::string getClassName() const { return "CMHaloIDTimeStepLinker";  }
    virtual std::string getCategory() const  { return "Viscontest2015";        }

protected:

    virtual void setDescriptions() { setDescription("Processor that links time steps and halo IDs."); }
    virtual void process();

    /**
     * Overwrites the base implementation of this function.
     * It is used to load the needed shader.
     * @see Processor
     */
    virtual void initialize();

    /**
     * Overwrites the base implementation of this function.
     * It is used to free the used shader.
     * @see Processor
     */
    virtual void deinitialize();

private:

    /**
     * Used to update interal datastructures if the time step property changes.
     */
    void timeStepChanged();

    /**
     * Used to update interal datastructures if the selected halo changes.
     */
    void selectedHaloChanged();

    /**
     * Get a pointer to the currently selected halo.
     */
    const CMHalo* getSelectedHalo();

    /**
     * Calculate and update the halo path if the time step or halo changes.
     */
    void calculateHaloPath();

    //-------------
    //  members
    //-------------

    /// Inport for halo data.
    CMHaloPort inport_;
    /// ID of the currently selected halo
    IntProperty selectedHaloIDProp_;
    /// ID of a hovered over halo
    IntProperty mouseOverHaloIDProp_;
    /// Current time step
    FloatProperty timeStep_;
    /// Whether or not to update the halo id if the time step changes
    BoolProperty updateHaloIDProp_;

    /// Halos that form the selected path
    std::deque<const CMHalo*> selectedPath_;
    /// Lock to prevent cyclic updating of property
    bool targetSelectionLock_;
};

} // namespace

#endif
