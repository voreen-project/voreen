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

#include "ensemblevolumepicking.h"
#include "../utils/samplepointconfigloader.h"

#include <set>
#include <vector>
#include <iostream>
#include <fstream>

namespace voreen {

    EnsembleVolumePicking::EnsembleVolumePicking()
        : Processor()
        //, poiStoragePort_(Port::INPORT, "poistorage_cpPort", "POIStorage Coprocessor", false)
        , ensembleInport_(Port::INPORT, "ensemble_inport", "Ensemble-Data Inport")
        , samplePointOutport_(Port::OUTPORT, "samplepoint_outport", "Samplepoint Outport")
        , mappingOutport_(Port::OUTPORT, "mapping_outport", "Samplepoint-Mapping Outport")
        , poiListOutport_(Port::OUTPORT, "poiList_outport", "POIList Outport")
        , poiListInport_(Port::INPORT, "poiList_inport", "POIList Inport")
        , memberListProp_("member_list_prop", "Selected Member")
        , fieldProp_("field_prop", "Selected Field")
        , spListProp_("sp_list_prop", "Samplepoints")
        , selVolumeProp_("selVolume_prop", "Selected volume", -1, -1, std::numeric_limits<int>::max() - 1)
        , defaultSpFileProp_("defaultSpFile_prop", "Configuration file", "Select file", VoreenApplication::app()->getUserDataPath(),
            "CSV (*.csv);;Text (*.txt)", FileDialogProperty::OPEN_FILE)
        , createDefaultSpButtonProp_("create_default_sp_prop", "Create samplepoints from file")
        , outportUpdatePolicyProp_("outport_update_policy_prop", "Mapping-Outport update policy")
        , updateEnsembleOutportProp_("update_ens_outport_prop", "Update Mapping-Outport")
        , samplePointFileProp_("sample_point_file_prop", "Samplepoint file", "Select file", VoreenApplication::app()->getUserDataPath(),
            "XML (*.xml);;Text (*.txt)", FileDialogProperty::SAVE_FILE)
        , saveSpButtonProp_("save_sp_button_prop", "Save samplepoints")
        , loadSpButtonProp_("load_sp_button_prop", "Load samplepoints")
        , selMember_("")
    {
        // Ports
        //addPort(poiStoragePort_);
        addPort(ensembleInport_);
        ON_CHANGE(ensembleInport_, EnsembleVolumePicking, ensembleChanged);
        addPort(samplePointOutport_);
        addPort(mappingOutport_);
        addPort(poiListOutport_);
        addPort(poiListInport_);

        // Properties
        addProperty(defaultSpFileProp_);
        addProperty(createDefaultSpButtonProp_);
        createDefaultSpButtonProp_.onClick(MemberFunctionCallback<EnsembleVolumePicking>(this, &EnsembleVolumePicking::createSpFromFile));
        addProperty(memberListProp_);
        addProperty(fieldProp_);
        addProperty(selVolumeProp_);
        addProperty(spListProp_);

        addProperty(outportUpdatePolicyProp_);
        outportUpdatePolicyProp_.addOption("Auto", "Auto", "Auto");
        outportUpdatePolicyProp_.addOption("Manual", "Manual", "Manual");
        outportUpdatePolicyProp_.setDefaultValue("Auto");
        ON_CHANGE(outportUpdatePolicyProp_, EnsembleVolumePicking, outportUpdatePolicyChanged);

        addProperty(updateEnsembleOutportProp_);
        updateEnsembleOutportProp_.setVisibleFlag(false);
        ON_CHANGE(updateEnsembleOutportProp_, EnsembleVolumePicking, updateEnsembleOutport);

        addProperty(samplePointFileProp_);

        addProperty(saveSpButtonProp_);
        saveSpButtonProp_.onClick(MemberFunctionCallback<EnsembleVolumePicking>(this, &EnsembleVolumePicking::saveSamplePoints));

        addProperty(loadSpButtonProp_);
        loadSpButtonProp_.onClick(MemberFunctionCallback<EnsembleVolumePicking>(this, &EnsembleVolumePicking::loadSamplePoints));

    }

    Processor* EnsembleVolumePicking::create() const {
        return new EnsembleVolumePicking();
    }

    std::string EnsembleVolumePicking::getClassName() const {
        return "EnsembleVolumePicking";
    }

    std::string EnsembleVolumePicking::getCategory() const {
        return "Ensemble Processing";
    }

    bool EnsembleVolumePicking::isReady() const {
        return ensembleInport_.isReady();
    }

    void EnsembleVolumePicking::setDescriptions() {
        setDescription("Processor which allows the user select an ensemblemember.  "
            " For the selected member a POIList is suppied at the POIList outport. The POI outport should be"
            " connected to the inport of a POIStorage-Processor. The Output of this POIStorage-Processor should be"
            " connected to the POIList inport of the processor. Additionally the Picked points of the selected member are"
            " avaiulable at the geometry-outport. "
            " The output is given in world-coordinates");
    }

    void voreen::EnsembleVolumePicking::serialize(Serializer& s) const {
        Processor::serialize(s);

        // Serialize picked samplepoints:
        s.serialize("mapping", mapping_);
    }

    void voreen::EnsembleVolumePicking::deserialize(Deserializer& s) {
        Processor::deserialize(s);

        // Deserialize picked samplepoints
        s.deserialize("mapping", mapping_);
    }

    void EnsembleVolumePicking::process() {
        tgtAssert(ensembleInport_.hasData(),  "inports not ready");
        bool updatedEnsembleOutport = false;

        const EnsembleMember* member = getSelectedMember();

        if (member == nullptr && selMember_ == "") {
            // Still no member selected
            return;
        }

        if (member == nullptr && selMember_ != "") {
            // Selection has been removed --> Clear member and samplepoint selection
            memberSelectionRemoved();
            return;
        }

        if (member != nullptr && selMember_ != member->getName()) {
            // A new member has been selected 
            selMember_ = member->getName();

            //clear samplepoint selection and reset outports
            memberSelectionChanged(member);
            
        }else {
            // Member selection did not change --> Update the POIList
            poiInportChanged();
            updatedEnsembleOutport = true;
        }

        // Update samplePointOutport_
        updateSamplePointOutport(member, !updatedEnsembleOutport);

    }

    void EnsembleVolumePicking::memberSelectionRemoved() {
        memberListProp_.setSelectedRowIndices(std::vector<int>());
        selMember_ = "";
        updateSpSelection();

        // Adjust volume-selector
        selVolumeProp_.reset();

        // Clear Outports
        samplePointOutport_.clear();
        poiListOutport_.setData(nullptr, true);
    }

    void voreen::EnsembleVolumePicking::memberSelectionChanged(const EnsembleMember* member) {
        updateSpSelection();

        // Adjust volume-selector
        const size_t numVolumes = member->getTimeSteps().size();
        selVolumeProp_.setMinValue(std::min(0, static_cast<int>(numVolumes) - 1));
        selVolumeProp_.setMaxValue(static_cast<int>(numVolumes) - 1);

        // Update POIList-Outport
        const POIList* list = new POIList(*mapping_.getPOIList(selMember_));
        poiListOutport_.setData(list, true);
    }

    void voreen::EnsembleVolumePicking::poiInportChanged() {
        const POIList *inList = poiListInport_.getData();
        if (inList == nullptr) {
            return;
        }
        // Replace POIList of selected member:
        mapping_.setPoiList(selMember_, inList);

        // Update ensemble-outport
        if (outportUpdatePolicyProp_.get() != "Manual") {
            mappingOutport_.setData(&mapping_, false);
        }

        // Update GUI
        updateSpSelection();
    }

    void voreen::EnsembleVolumePicking::updateSamplePointOutport(const EnsembleMember* member, const bool updateEnsembleOutport) {
        const POIList* samplepoints = mapping_.getPOIList(member->getName());
        std::vector<tgt::vec3> positions;
        for (POIPoint sp : samplepoints->getPoints()) {
            tgt::vec3 pos = sp.position_;
            positions.push_back(pos);
        }

        PointListGeometryVec3* result = new PointListGeometryVec3();
        result->setData(positions);
        samplePointOutport_.setData(result, true);

        // Update ensemble-outport
        if (outportUpdatePolicyProp_.get() != "Manual" && updateEnsembleOutport) {
            mappingOutport_.setData(&mapping_, false);
        }
    }

    void EnsembleVolumePicking::createSpFromFile() {
        if (defaultSpFileProp_.get() == "") {
            return;
        }

        SamplePointConfigLoader loader;
        std::vector<SamplePointConfig> config =  loader.loadSamplepointConfigFile(defaultSpFileProp_.get());

        // Create Samplepoints
        mapping_.addDefaultPOIPoint(config);

        // Reset member-selection
        memberListProp_.setSelectedRowIndices(std::vector<int>());

        // Invalidate processor
        invalidate();

    }

    void EnsembleVolumePicking::ensembleChanged() {

        // Clear Propertys
        memberListProp_.reset();
        fieldProp_.reset();
        spListProp_.reset();

        if (!ensembleInport_.isReady())
            return;

        const EnsembleDataset* dataset = ensembleInport_.getData();

        // Populate field-property
        fieldProp_.blockCallbacks(true);
        for (const std::string& fieldName : dataset->getCommonFieldNames()) {
            if (!fieldProp_.hasKey(fieldName))
                fieldProp_.addOption(fieldName, fieldName, fieldName);
        }
        fieldProp_.blockCallbacks(false);

        // Populate member-list
        memberListProp_.blockCallbacks(true);
        for (const EnsembleMember& member : dataset->getMembers()) {
            memberListProp_.addRow(member.getName(), member.getColor());
        }
        memberListProp_.blockCallbacks(false);

        // Check if deserialized mapping is compatible with loaded ensemble
        if (ensembleMatchesMapping(&mapping_)) {
            return;
        }

        // Create new SamplePoint mapping and assign to outport
        mapping_ = EnsembleSamplePointMapping(*dataset);
        mappingOutport_.setData(&mapping_, false);
    }

    const EnsembleMember* EnsembleVolumePicking::getSelectedMember() const {
        if(!ensembleInport_.hasData() || memberListProp_.getSelectedRowIndices().empty()){
            return nullptr;
        }
        const int idx = memberListProp_.getSelectedRowIndices().front();
        const EnsembleDataset* dataset = ensembleInport_.getData();
        return &(dataset->getMembers()[idx]);
    }

    void EnsembleVolumePicking::updateSpSelection() {
        spListProp_.reset();

        if (selMember_ != "") {
            // Fill samplepoint-list
            const POIList* poiList = mapping_.getPOIList(selMember_);
            for (const POIPoint& sp : poiList->getPoints()) {
                // Build label: GROUPNAME::ID
                const std::string  id = std::to_string(sp.id_);
                const std::string group = poiList->getGroupName(sp.group_);
                spListProp_.addRow(group + "::" + id);
            }
        }
    }

    void voreen::EnsembleVolumePicking::outportUpdatePolicyChanged() {
        if (outportUpdatePolicyProp_.get() == "Manual") {
            updateEnsembleOutportProp_.setVisibleFlag(true);
        }else {
            updateEnsembleOutportProp_.setVisibleFlag(false);
            updateEnsembleOutport();
        }
    }

    void voreen::EnsembleVolumePicking::updateEnsembleOutport() {
        mappingOutport_.setData(&mapping_, false);
    }

    bool voreen::EnsembleVolumePicking::ensembleMatchesMapping(const EnsembleSamplePointMapping* mapping) {
        if (!ensembleInport_.isReady() || mapping == nullptr) {
            return false;
        }

        const EnsembleDataset* dataset = ensembleInport_.getData();

        // Create set containing all names of the ensemble
        std::set<std::string> datasetNames;
        for (const EnsembleMember& member : dataset->getMembers()) {
            datasetNames.insert(member.getName());
        }

        // Get list of all names in the mapping
        const std::vector<std::string> mappingNames = mapping->getMappedMemberNames();

        if (datasetNames.size() != mappingNames.size()) {
            return false;
        }

        for (const std::string& mappingName : mappingNames) {
            if (datasetNames.find(mappingName) == datasetNames.end()) {
                // Mapped member not part of ensemble
                return false;
            }
        }

        return true;
    }

    void voreen::EnsembleVolumePicking::saveSamplePoints() {
        if (samplePointFileProp_.get().empty()) {
            return;
        }

        // Serialize mapping
        XmlSerializer serializer(samplePointFileProp_.get());
        serializer.setUseAttributes(true);
        serializer.serialize("mapping", &mapping_);
        if (!serializer.getErrors().empty()) {
            LERROR("Failed to serialize samplepoints!");
            return;
        }

        // Write json to file
        try {
            std::ofstream outStream;
            outStream.open(serializer.getDocumentPath());
            serializer.write(outStream);
            outStream.close();
        } catch (std::ofstream::failure e) {
            LERROR("Failed to save samplepoints!");
            return;
        }

        LINFO("Saved samplepoints to " << serializer.getDocumentPath());
    }

    void voreen::EnsembleVolumePicking::loadSamplePoints() {
        if (samplePointFileProp_.get().empty()) {
            return;
        }

        // Derialize mapping
        EnsembleSamplePointMapping* mapping;
        try {
            std::ifstream inStream(samplePointFileProp_.get());
            XmlDeserializer deserializer(samplePointFileProp_.get());
            deserializer.setUseAttributes(true);
            deserializer.read(inStream);
            deserializer.deserialize("mapping", mapping);
        }
        catch (...) {
            LERROR("Could not load Samplepoints");
            return;
        }

        // Check if mapping is compatible with loaded ensemble
        if (ensembleMatchesMapping(mapping)) {
            mapping_ = *mapping;
            // Remove member-selaction. Otherwise the update from the POI-Inport 
            // will overwirte the deserialized mapping
            memberListProp_.setSelectedRowIndices(std::vector<int>());
            LINFO("Loaded samplepoints from " << samplePointFileProp_.get());
        }else {
            LERROR("The loaded samplepoint-mapping is not compatible with the loaded ensemble.");
        }
    }

    void voreen::EnsembleVolumePicking::onEvent(tgt::Event* e) {
        tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
        if (event) {
            mouseEvent(event);
        }
        Processor::onEvent(e);
    }

    void voreen::EnsembleVolumePicking::mouseEvent(tgt::MouseEvent* e) {
        if (e->button() == tgt::MouseEvent::MOUSE_WHEEL_UP 
            && e->action() == tgt::MouseEvent::WHEEL
            && e->modifiers() == tgt::Event::Modifier::SHIFT){
            // Switch to next timestep
            const int cur_timestep = selVolumeProp_.get();
            if (cur_timestep < selVolumeProp_.getMaxValue()) {
                selVolumeProp_.set(cur_timestep + 1);
            }
        }

        if (e->button() == tgt::MouseEvent::MOUSE_WHEEL_DOWN
            && e->action() == tgt::MouseEvent::WHEEL
            && e->modifiers() == tgt::Event::Modifier::SHIFT) {
            // Switch to previous timestep
            const int cur_timestep = selVolumeProp_.get();
            if (cur_timestep > selVolumeProp_.getMinValue()) {
                selVolumeProp_.set(cur_timestep - 1);
            }
        }

    }
} // namespace
