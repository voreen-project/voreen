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

#include <voreen/core/datastructures/volume/volumeatomic.h>
#include <voreen/core/processors/processor.h>
#include <modules/ensembleanalysis/ports/ensembledatasetport.h>
#include "voreen/core/properties/numeric/intervalproperty.h"

#ifndef VRN_FEATUREEXTRACTOR_H
#define VRN_FEATUREEXTRACTOR_H

namespace voreen {

    class Volume;

    class VRN_CORE_API FeatureExtractor : public Processor {
    public:
        FeatureExtractor();
        ~FeatureExtractor();
        virtual Processor* create() const;

        virtual std::string getCategory() const   { return "Volume Processing"; }
        virtual std::string getClassName() const  { return "FeatureExtractor";      }
        virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;   }

    protected:
        virtual void setDescriptions() {
            setDescription("Extracts the Features of the EnsambleDataSet in input and puts out the extracted Features in a VolumeList");
        }

        virtual void process();
        virtual void applyChanges();

    private:
        void initProperties();
        /// Inport for the ensemble data structure.
        EnsembleDatasetPort ensembleInport_;

        /// The extracted volume data.
        VolumeListPort outport_;

        std::map<std::string, FloatIntervalProperty> fieldProperties_;

        ButtonProperty apply_;

        static const std::string loggerCat_; ///< category used in logging
    };

}   //namespace

#endif // VRN_FEATUREEXTRCATOR_H
