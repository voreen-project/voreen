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

#ifndef VRN_VOLUMELISTMODALITYFILTER
#define VRN_VOLUMELISTMODALITYFILTER

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/datastructures/volume/volumelist.h"

namespace voreen {

/**
 * Permits to select one specific modality among those of all volume
 * within the incoming VolumeList. The filtered VolumeList contains
 * only volumes of the specified modality.
 */
class VRN_CORE_API VolumeListModalityFilter : public Processor {

public:
    VolumeListModalityFilter();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "VolumeListModalityFilter"; }
    virtual std::string getCategory() const  { return "Utility";                  }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;          }
    virtual bool isUtility() const           { return true;                       }

    virtual void invalidate(int inv = INVALID_RESULT);

protected:
    virtual void setDescriptions() {
        setDescription("Permits to select one specific modality among those of all volumes within the incoming volume list. The filtered volume list contains only volumes of the specified modality.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    VolumeListPort inport_;
    VolumeListPort outport_;

    OptionProperty<Modality> modalityProp_;
    Modality currentModality_;
    VolumeList filteredList_;

    static const std::string loggerCat_;

private:
    void adjustFilteredList();

};

}   // namespace voreen

#endif
