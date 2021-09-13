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
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. * *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef VRN_CMLAZYLINKER_H
#define VRN_CMLAZYLINKER_H


#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/numericproperty.h"
#include "../ports/cmhaloport.h"


//use namespace voreen
namespace voreen {

/**
 * A processor to lazily link input and output properties. Values will be copied from
 * input to output only when the *Event_-properties change.
 */
class VRN_CORE_API CMLazyLinker : public Processor {

public:
    /**
     * Constructor
     */
    CMLazyLinker();
    ~CMLazyLinker();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CMLazyLinker();     }
    virtual std::string getClassName() const { return "CMLazyLinker";  }
    virtual std::string getCategory() const  { return "Viscontest2015";        }

protected:

    virtual void setDescriptions() { setDescription("Processor that only transmits values if events occur"); }
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
     * Will be executed to link in- and out-properties (write value from int to out) when intEvent_ changes.
     */
    void eventTriggered();
    //-------------
    //  members
    //-------------
    /// Input int property
    IntProperty intIn_;
    /// Output int property
    IntProperty intOut_;
    /// Input float property
    FloatProperty floatIn_;
    /// Output float property
    FloatProperty floatOut_;
    /// Input camera property
    CameraProperty cameraIn_;
    /// Output cameraint property
    CameraProperty cameraOut_;
    /// Property that triggers eventTriggered on Change
    IntProperty intEvent_;
};

} // namespace

#endif
