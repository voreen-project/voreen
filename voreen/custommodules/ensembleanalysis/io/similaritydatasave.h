/***********************************************************************************
*                                                                                 *
* Voreen - The Volume Rendering Engine                                            *
*                                                                                 *
* Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_SIMILARITYDATASAVE_H
#define VRN_SIMILARITYDATASAVE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/progressproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

    class VRN_CORE_API SimilarityDataSave : public Processor {
    public:
        /** Constructior */
        SimilarityDataSave();
        /** Destructor */
        virtual ~SimilarityDataSave();

        virtual Processor* create() const         { return new SimilarityDataSave(); }
        virtual std::string getClassName() const  { return "SimilarityDataSave"; }
        virtual std::string getCategory() const   { return "Output"; }
        virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }
        virtual bool isEndProcessor() const       { return true; }

    protected:
        virtual void setDescriptions() {
        }

        //--------------------------
        //  Override
        //--------------------------
        /** Saves the volume file. */
        virtual void process() override;
        /** Processor should be ready, even if no volume is attached */
        virtual bool isReady() const override;
        /** Triggers a auto-save if the file property has changed values. */
        virtual void invalidate(int inv = 1) override;

        //--------------
        //  Callbacks
        //--------------
        /** Main function used to save the field plot. Triggered by saveButton_. */
        void saveSimilarityData();
        /** Adjusts propertiy ranges to the incoming ensemble dataset. */
        void adjustToEnsemble();

        //--------------
        //  Member
        //--------------
    private:
        EnsembleDatasetPort ensembleInport_;

        FileDialogProperty filenameProp_;          ///< determines the name of the saved file
        ButtonProperty saveButton_;                ///< triggers a save
        FloatProperty temporalResolution_;         ///< temporal resolution, simulation time will be devided into
        ProgressProperty progress_;                ///< progress of saving files

        bool saveSimilarityData_;                  ///< used to determine, if process should save or not

        static const std::string loggerCat_;
    };

}   //namespace

#endif
