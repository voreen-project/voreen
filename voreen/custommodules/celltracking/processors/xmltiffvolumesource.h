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

#ifndef VRN_XMLTIFFVOLUMESOURCE_H
#define VRN_XMLTIFFVOLUMESOURCE_H

#include <string>
#include <vector>
#include <map>

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"


namespace voreen {

class VRN_CORE_API XMLTiffVolumeSource : public Processor {

public:
    XMLTiffVolumeSource();
    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;

protected:
    virtual void setDescriptions();
    virtual void process();

private:
    struct SliceDesc{
        SliceDesc(int z, int channel, int timeStep);
        bool operator<(const SliceDesc & other) const;

        int z_;
        int channel_;
        int timeStep_;

    };

    typedef std::map<XMLTiffVolumeSource::SliceDesc, std::string> DocumentDesc;

    VolumePort outport_;

    FileDialogProperty selectedFile_;
    IntProperty timeStep_;
    IntProperty channel_;
    BoolProperty autoload_;
    ButtonProperty loadButton_;

    bool shouldLoadFileOnProcess_;

    void hasSelectedFile(void);
    void loadButtonPressed(void);
    void loadXMLFile(const std::string path);
    Volume* loadVolume();

    

    DocumentDesc documentDesc_;
    int sizeZ_;
    int sizeX_;
    int sizeY_;
};

} // namespace

#endif
