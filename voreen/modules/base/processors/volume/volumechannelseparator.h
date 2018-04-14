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

#ifndef VRN_VOLUMECHANNELSEPARATOR_H
#define VRN_VOLUMECHANNELSEPARATOR_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VolumeChannelSeparator : public CachingVolumeProcessor {
public:

    VolumeChannelSeparator();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeChannelSeparator";   }
    virtual std::string getCategory() const   { return "Volume Processing";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Separates a multi-channel volume (up to four channels) into several single-channel volumes.");
    }

    virtual void process();

    /**
     * Has to be overwritten from CachingVolumeProcessor to handle empty outputs before pressing the refresh button.
     */
    virtual void beforeProcess();

    /**
     * Has to be overwritten from CachingVolumeProcessor to handle empty outputs before pressing the refresh button.
     */
    virtual void afterProcess();

    virtual void initialize();
    virtual void deinitialize();

    /**
     * Separates the channels of a multi-channel volume into several single-channel volumes that are newly created.
     */
    std::vector<Volume*> separateChannels(const VolumeBase* volume);

    template<class T>
    std::vector<VolumeRAM*> separateVolumeChannels(const VolumeBase* volume) {

        size_t numChannels = volume->getNumChannels();
        tgt::svec3 dimensions = volume->getDimensions();
        std::vector<VolumeRAM*> separatedData;

        //create new VolumeAtomic representations
        std::vector<VolumeAtomic<T>*> volumeAtomics;
        for (size_t i = 0; i < numChannels; ++i)
            volumeAtomics.push_back(new VolumeAtomic<T>(dimensions));

        //get VolumeRAM representation
        const VolumeRAM* ram = volume->getRepresentation<VolumeRAM>();

        if (!ram) {
            for (size_t i = 0; i < numChannels; ++i) {
                delete volumeAtomics.at(i);
            }
            LERROR("Could not get VolumeRAM representation from input volume");
            return separatedData;
        }

        //try to cast into volume atomic type
        if (numChannels == 2) {
            const VolumeAtomic<tgt::Vector2<T> >* inputAtomic = dynamic_cast<const VolumeAtomic<tgt::Vector2<T> >* >(ram);
            if (!inputAtomic) {
                LERROR("Could not cast VolumeRAM to VolumeAtomic for input volume.");
                for (size_t i = 0; i < numChannels; ++i) {
                    delete volumeAtomics.at(i);
                }
                return separatedData;
            }

            //copy data of the input volume into the separated data representations
            for (tgt::svec3 i = tgt::svec3(0,0,0); i.z < dimensions.z; ++i.z) {
                //set progress after each slice
                setProgress(static_cast<float>(i.z) / static_cast<float>(dimensions.z));

                for (i.y = 0; i.y < dimensions.y; ++i.y) {
                    for (i.x = 0; i.x < dimensions.x; ++i.x) {

                        for (size_t c = 0; c < numChannels; ++c) {
                            volumeAtomics.at(c)->voxel(i) = inputAtomic->voxel(i)[c];
                        }
                    }
                }
            }

            for (size_t c = 0; c < numChannels; ++c) {
                separatedData.push_back(dynamic_cast<VolumeRAM*>(volumeAtomics.at(c)));
            }
        }
        else if (numChannels == 3) {
            const VolumeAtomic<tgt::Vector3<T> >* inputAtomic = dynamic_cast<const VolumeAtomic<tgt::Vector3<T> >* >(ram);
            if (!inputAtomic) {
                LERROR("Could not cast VolumeRAM to VolumeAtomic for input volume.");
                for (size_t i = 0; i < numChannels; ++i) {
                    delete volumeAtomics.at(i);
                }
                return separatedData;
            }

            //copy data of the input volume into the separated data representations
            for (tgt::svec3 i = tgt::svec3(0,0,0); i.z < dimensions.z; ++i.z) {
                //set progress after each slice
                setProgress(static_cast<float>(i.z) / static_cast<float>(dimensions.z));

                for (i.y = 0; i.y < dimensions.y; ++i.y) {
                    for (i.x = 0; i.x < dimensions.x; ++i.x) {

                        for (size_t c = 0; c < numChannels; ++c) {
                            volumeAtomics.at(c)->voxel(i) = inputAtomic->voxel(i)[c];
                        }
                    }
                }
            }

            for (size_t c = 0; c < numChannels; ++c) {
                separatedData.push_back(dynamic_cast<VolumeRAM*>(volumeAtomics.at(c)));
            }
        }
        else if (numChannels == 4) {
            const VolumeAtomic<tgt::Vector4<T> >* inputAtomic = dynamic_cast<const VolumeAtomic<tgt::Vector4<T> >* >(ram);
            if (!inputAtomic) {
                LERROR("Could not cast VolumeRAM to VolumeAtomic for input volume.");
                for (size_t i = 0; i < numChannels; ++i) {
                    delete volumeAtomics.at(i);
                }
                return separatedData;
            }

            //copy data of the input volume into the separated data representations
            for (tgt::svec3 i = tgt::svec3(0,0,0); i.z < dimensions.z; ++i.z) {
                //set progress after each slice
                setProgress(static_cast<float>(i.z) / static_cast<float>(dimensions.z));

                for (i.y = 0; i.y < dimensions.y; ++i.y) {
                    for (i.x = 0; i.x < dimensions.x; ++i.x) {

                        for (size_t c = 0; c < numChannels; ++c) {
                            volumeAtomics.at(c)->voxel(i) = inputAtomic->voxel(i)[c];
                        }
                    }
                }
            }

            for (size_t c = 0; c < numChannels; ++c) {
                separatedData.push_back(dynamic_cast<VolumeRAM*>(volumeAtomics.at(c)));
            }
        }
        else {
            LERROR("invalid number of channels");
            for (size_t i = 0; i < numChannels; ++i) {
                delete volumeAtomics.at(i);
            }
            return separatedData;
        }

        setProgress(1.f);
        return separatedData;
    }

    virtual void refresh();

private:

    VolumePort volumeInport_;

    VolumePort volumeOutport_;
    VolumePort volumeOutport2_;
    VolumePort volumeOutport3_;
    VolumePort volumeOutport4_;

    BoolProperty autoRefresh_;          ///< if enabled, output volumes are re-computed as soon as the input changes. Else the refresh-button has to be used.
    ButtonProperty refresh_;            ///< if pressed this clears the outports and invalidates the processor so that the output volumes are newly computed

    bool refreshNecessary_;             ///< if set to true, process will compute the new output volumes

    static const std::string loggerCat_;
};

} //namespace

#endif
