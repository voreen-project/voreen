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

#ifndef VRN_VOLUMELISTOFFSET_H
#define VRN_VOLUMELISTOFFSET_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

class VRN_CORE_API VolumeListOffset : public VolumeProcessor, public PortObserver {
public:
    VolumeListOffset();
    virtual ~VolumeListOffset();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "VolumeListOffset"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;  }

protected:
    virtual void setDescriptions() {
        setDescription("Sets the volumes' offset for each volume in the list. The output is a new VolumeList, but each of the volumes still references the volume of the original list, so that memory requirements are not doubled. The volumes' transformation matrix is not changed.");
    }

    virtual void initialize();

    /// deletes the volume list and the decorator volumes
    virtual void deinitialize();

    virtual void process();

    virtual void adjustPropertiesToInput();

    virtual void updateCurrentlySelected();

    virtual void clearOutput();

    virtual void afterConnectionAdded(const Port* source, const Port* connectedPort) {}
    virtual void beforeConnectionRemoved(const Port* source, const Port*) {}
    virtual void dataWillChange(const Port* source);
    //virtual void dataHasChanged(const Port* source) {};

private:
    void resetOffset();

    void adjustPropertyVisibility();

    VolumeListPort inport_;
    VolumeListPort outport_;

    BoolProperty enableProcessing_;
    FloatVec3Property offset_;
    ButtonProperty reset_;

    IntProperty currentlySelected_;
    FloatVec3Property offsetDisplay_;

    std::vector<std::unique_ptr<VolumeBase>> decorators_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_VOLUMELISTOFFSET_H
