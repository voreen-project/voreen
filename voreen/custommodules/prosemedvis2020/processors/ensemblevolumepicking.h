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

#ifndef VRN_ENSEMBLEVOLUMEPICKING_H
#define VRN_ENSEMBLEVOLUMEPICKING_H

#include "voreen/core/ports/genericcoprocessorport.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "../modules/ensembleanalysis/ports/ensembledatasetport.h"
#include "../modules/poi/processors/poistorage.h"
#include "../modules/poi/datastructures/poilist.h"
#include "../ports/ensemblesamplepointmappingport.h"

namespace voreen {

    /**
     *  Processor, which lets the user pick sample-points for different members of an ensemble.
     */
    class VRN_CORE_API EnsembleVolumePicking : public Processor {
    public:

        EnsembleVolumePicking();
        virtual Processor* create() const;
        virtual std::string getClassName() const;
        virtual std::string getCategory() const;

        virtual void serialize(Serializer& s) const;
        virtual void deserialize(Deserializer& s);

    protected:

        virtual void setDescriptions();
        virtual void process();
        virtual void onEvent(tgt::Event* e);
        virtual bool isReady() const;

    private:

        EnsembleSamplePointMapping mapping_;

        std::string selMember_;

        // Ports
        EnsembleDatasetPort ensembleInport_;  ///< inport used to recieve the ensemble
        GeometryPort samplePointOutport_;     ///< outport containing the current sample-points
        EnsembleSamplePointMappingPort mappingOutport_; // Outport containing the SamplePoint mapping
        POIListPort poiListOutport_;            // Outport containing the POIList of the selected Member
        POIListPort poiListInport_;            // Inport which recieves the new, modified, POIList 

        // Propertys
        StringListProperty memberListProp_;    // List used to select a member of the ensemble
        OptionProperty<std::string> fieldProp_;// Property used to select a field from the ensemble
        StringListProperty spListProp_;    // List used to select a samplepoint
        IntProperty selVolumeProp_;      // The selected volume (timestep)
        FileDialogProperty defaultSpFileProp_; // Property used to select a file from which default-samplepoints are created.
        ButtonProperty createDefaultSpButtonProp_; // Button used to created default-points from the loaded file.
        OptionProperty<std::string> outportUpdatePolicyProp_;// Property used to decide when the ensemblemapping-outport is updated
        ButtonProperty updateEnsembleOutportProp_; // Button used to update the ensemblemapping-outport

        FileDialogProperty samplePointFileProp_;    // Property used to select a file which is used to save/load the samplepoints
        ButtonProperty saveSpButtonProp_;           // Button to save the samplepoints to the file specified by samplePointFileProp_
        ButtonProperty loadSpButtonProp_;           // Button to load the samplepoints from the file specified by samplePointFileProp_

        /** Called when the input of ensembleInport_ changes */
        void ensembleChanged();

        /** Called when the user presses the 'Create from file'-Button */
        void createSpFromFile();

        /** Called when an mouse-event happens */
        void mouseEvent(tgt::MouseEvent* e);

        /**
         *  Returns the currently selected Member. 
         *  (Or nullptr if no member is selected)
         */
        const EnsembleMember* getSelectedMember() const;

        /**
         *  Updates the properties which display the samplepoints for the selected member.
         */
        void updateSpSelection();

        /**
         *  Called when the selection in the memberlist is removed
         */
        void memberSelectionRemoved();

        /**
         *  Called when the selection of outportUpdatePolicyProp_ changes
         */
        void outportUpdatePolicyChanged();

        /**
         * Called when the updateEnsembleOutportProp_-Button is clicked
         */
        void updateEnsembleOutport();

        /**
         *  Called when a new member is selected in the memberlist.
         */
        void memberSelectionChanged(const EnsembleMember* member);

        /**
         * Updates the samplePointOutport.
         */
        void updateSamplePointOutport(const EnsembleMember* member, const bool updateEnsembleOutport);

        /**
         *  Called when the POIListInport changed
         */
        void poiInportChanged();

        /**
         *  Checks if the loaded mapping is compatible to the current ensemble.
         *  Use to decide if the deserialized mapping should be kept or replaced.
         */
        bool ensembleMatchesMapping(const EnsembleSamplePointMapping* mapping);

        /**
         *  Save the samplepoints to a file
         */
        void saveSamplePoints();

        /**
         *  Load samplepoints from a file
         */
        void loadSamplePoints();

    };

} // namespace

#endif // VRN_ENSEMBLEVOLUMEPICKING_H
