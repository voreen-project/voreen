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

#include <modules/ensembleanalysis/utils/ensemblehash.h>
#include "featureextractor.h"
#include <include/voreen/core/datastructures/volume/volumeminmax.h>

namespace voreen {

    const std::string FeatureExtractor::loggerCat_("voreen.FeatureExtractor");

    FeatureExtractor::FeatureExtractor()
        : Processor()
        , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
        , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", true)
        , apply_("apply", "Apply changes", Processor::INVALID_RESULT)
    {
        addPort(outport_);
        addPort(ensembleInport_);
        ON_CHANGE(ensembleInport_, FeatureExtractor, initProperties);

        addProperty(apply_);
        ON_CHANGE(apply_, FeatureExtractor, applyChanges);

        std::vector<std::string> fieldNames;
        fieldNames.push_back("angle");
        fieldNames.push_back("magnitude");
        fieldNames.push_back("radius");
        fieldNames.push_back("spin transition-induced density anomaly");
        fieldNames.push_back("temperature");
        fieldNames.push_back("temperature anomaly");
        fieldNames.push_back("thermal conductivity");
        fieldNames.push_back("thermal expansivity");
        fieldNames.push_back("vx");
        fieldNames.push_back("vy");
        fieldNames.push_back("vz");

        //loop durch alle fieldnames
        for (const std::string& fieldName : fieldNames) {
            //add Properties
            fieldProperties_[fieldName] = FloatIntervalProperty(fieldName, fieldName, 0.0, std::numeric_limits<float>::max(), std::numeric_limits<float>::min());
            fieldProperties_[fieldName].setNumDecimals(10);
            fieldProperties_[fieldName].setVisibleFlag(false);

            addProperty(fieldProperties_[fieldName]);
        }
    }

    FeatureExtractor::~FeatureExtractor() {

    }

    Processor* FeatureExtractor::create() const {
        return new FeatureExtractor();
    }

    void FeatureExtractor::process() {
    }

    void FeatureExtractor::initProperties() {
        if (!ensembleInport_.hasData())
            outport_.setData(nullptr);
        else {
            for (auto &pair : fieldProperties_) {
                pair.second.setVisibleFlag(false);
            }

            //propertiy max min values setzten und Properties sichtbar machen
            for(auto& name : ensembleInport_.getData()->getCommonFieldNames()) {
                fieldProperties_[name].setVisibleFlag(true);
                fieldProperties_[name].setMinValue(ensembleInport_.getData()->getValueRange(name).x);
                fieldProperties_[name].setMaxValue(ensembleInport_.getData()->getValueRange(name).y);
            }

            std::cout << "<Feature Extractor> Initialising finished.\n";
        }
    }

    void FeatureExtractor::applyChanges() {
        if (!ensembleInport_.hasData())
            outport_.setData(nullptr);
        else {
            //daten bearbeitbar machen
            std::unique_ptr<EnsembleDataset> ensemble_(new EnsembleDataset(*ensembleInport_.getData()));
            
            const EnsembleDataset& ensemble = *ensemble_;

            //outputlist erstellen
            VolumeList* volumeList = new VolumeList();

            //dimensionen zwischenspeichern
            tgt::svec3 dimensions = ensemble.getVolumes()[0]->getDimensions();
            std::cout << "<Feature Extractor> Dimensions: " << dimensions << "\n";
            std::cout << "<Feature Extractor> Volumesize: " << ensemble.getVolumes().size() << "\n";

            //durch alle member und timesteps loopen
            for (const EnsembleMember& member : ensemble.getMembers()) {
                for (size_t t = 0; t < member.getTimeSteps().size(); t++) {
                    std::cout << "<Feature Extractor> Timestep: " << t << "\n";

                    //vereinfachungsvariablen fuer den zugriff
                    size_t fieldNameSize = member.getTimeSteps()[t].getFieldNames().size();
                    size_t i = t * fieldNameSize;

                    //map die als key ein fieldname nimmt und als value ein Volumen in der RAM_Float representation
                    std::map<std::string, VolumeRAM_Float*> volumeMap;

                    std::cout << "<Feature Extractor> Fieldnamesize: " << fieldNameSize << "\n";
                    //fill map von oben mit den konvertierten Volumen und passenden fieldnames
                    for (int f = 0; f < fieldNameSize; f++) {
                        std::string fieldName = member.getTimeSteps()[t].getFieldNames()[f];
                        VolumeRAMRepresentationLock data(ensemble.getVolumes()[i + f]);
                        volumeMap[fieldName] = dynamic_cast<VolumeRAM_Float*>(data->clone());
                    }

                    std::cout << "<Feature Extractor> VolumeMap Size: " << volumeMap.size() << "\n";

                    //erstell das outputVolume
                    VolumeRAM* outputVolume = 0;
                    VolumeRAM_Float* target = new VolumeRAM_Float(dimensions);
                    outputVolume = target;

                    //hilfsvariablen
                    float voxelValue = 0;
                    float propertyMin = 0;
                    float propertyMax = 0;


                    bool firstRun = true;
                    //die hauptschleife
                    //durch alle fields laufen, indem durch die volumeMap gelooped wird
                    for (auto const& pair : volumeMap) {
                        propertyMin = fieldProperties_[pair.first].getInterval().get().x;
                        propertyMax = fieldProperties_[pair.first].getInterval().get().y;
                        //hier wird durch die dimensionen gelooped
                        //paralellisieren durch #pragma omp parallel for
//#pragma omp parallel for
                        //for (long z = 0; z < static_cast<long>(dimensions.z); z++) {

                        for (size_t z = 0; z < dimensions.z; z++) {
                            for (size_t y = 0; y < dimensions.y; y++) {
                                for (size_t x = 0; x < dimensions.x; x++) {
                                    //setzte Hilfsvariablen
                                    voxelValue = pair.second->voxel(x, y, z);


                                    //guck ob value innerhalb des Intervalls ist
                                    if ((firstRun || target->voxel(x, y, z) != 0) && voxelValue >= propertyMin && voxelValue <= propertyMax) {
                                        target->voxel(x, y, z) = 1;
                                    }
                                    else {
                                        target->voxel(x, y, z) = 0;
                                    }

                                }    
                            }
                        }
                        firstRun = false;

                        std::cout << "<Feature Extractor> Done with Volume: " << pair.first << "\n";
                    }
                    //output Volume zur output Liste hinzufuegen
                    volumeList->add(new Volume(outputVolume, ensemble.getVolumes()[0]->getSpacing(), ensemble.getVolumes()[0]->getOffset()));
                    
                    std::cout << "<Feature Extractor> Done with Timestep: " << t << "\n";
                }
            }
            std::cout << "<Feature Extractor> Done." << "\n";
            //output setzten
            outport_.setData(volumeList, true);
        }
    }
}   // namespace
