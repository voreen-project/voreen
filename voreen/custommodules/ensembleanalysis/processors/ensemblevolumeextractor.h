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

#ifndef VRN_ENSEMBLEVOLUMEEXTRACTOR_H
#define VRN_ENSEMBLEVOLUMEEXTRACTOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/genericport.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

class EnsembleDataset;
class Volume;
class VolumeList;

/**
 * Allows access to the volumelist contained in an ensemble dataset.
 */
class VRN_CORE_API EnsembleVolumeExtractor : public Processor {

public:
    EnsembleVolumeExtractor();
    virtual ~EnsembleVolumeExtractor();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "EnsembleVolumeExtractor"; }
    virtual std::string getCategory() const   { return "Input";                   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;   }

protected:
    virtual void setDescriptions() {
        setDescription("Extracts the volume list inside an ensemble dataset.");
    }

    void process();

    /// The ensemble data.
    EnsembleDatasetPort inport_;

    /// The extracted volume data.
    VolumeListPort outport_;

    static const std::string loggerCat_;

private:
    void updateOutput();

};

} // namespace

#endif
