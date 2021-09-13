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

#ifndef VRN_COSMOLOGYDATASOURCE_H
#define VRN_COSMOLOGYDATASOURCE_H

#include "boost/iostreams/device/mapped_file.hpp"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "../ports/cmparticleport.h"
#include "../datastructures/cmparticledata.h"

namespace voreen{
class CosmologyDataSource : public Processor{
public:
    CosmologyDataSource();
    virtual ~CosmologyDataSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "CosmologyDataSource";            }
    virtual std::string getCategory() const     { return "Viscontest2019";                 }
    virtual void setDescriptions()              { setDescription("Cosmology Data Source"); }
    virtual CodeState   getCodeState() const    { return CODE_STATE_EXPERIMENTAL;          }

  
protected:
    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

private:
    CMParticlePort     outport_;
    FileDialogProperty fileFolder_;
	bool read = false;

    CMParticleData *          readParticleData();
    CMParticleDataTimeSlice * tryReadFile(const std::string& filename);
    
    static const std::string loggerCat_; ///< category used in logging
};
}

#endif
