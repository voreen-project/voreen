/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMECHANNELMERGER_H
#define VRN_VOLUMECHANNELMERGER_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class VolumeChannelMerger : public CachingVolumeProcessor {
public:

    VolumeChannelMerger();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeChannelMerger";   }
    virtual std::string getCategory() const   { return "Volume Processing";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Merges up to four single-channel volumes into one multi-channel volume.");
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
     * Merges several single-channel volumes into one multi-channel volume that is newly created.
     * The input volumes should have the same dimensions, spacing, offset and transformation matrices, since only the meta data of the first volume is used.
     */
    Volume* mergeVolumes(std::vector<const VolumeBase*> volumes);

    /**
     * Tries to cast VolumeBase representations into VolumeAtomic<T> representation and creates a new VolumeAtomic<tgt::Vector"N"<T> > representation where "N" is the number of channels.
     */
    template<class T>
    VolumeRAM* mergeVolumeRepresentations(std::vector<const VolumeBase*> volumes) {

        size_t numVolumes = volumes.size();
        tgt::svec3 dimensions = volumes.at(0)->getDimensions();
        VolumeRAM* mergedData = 0;

        //try to cast into VolumeAtomic representation
        const VolumeAtomic<T>** atomicArray = new const VolumeAtomic<T>*[numVolumes];

        for (size_t i = 0; i < numVolumes; ++i) {

            const VolumeRAM* ram = volumes.at(i)->getRepresentation<VolumeRAM>();

            if (!ram) {
                delete[] atomicArray;
                LERROR("Could not get VolumeRAM representation from input volume " << i);
                return 0;
            }

            atomicArray[i] = dynamic_cast<const VolumeAtomic<T>* >(ram);
            if (!(atomicArray[i])) {
                delete[] atomicArray;
                LERROR("Could not cast VolumeRAM to VolumeAtomic type for input volume " << i);
                return 0;
            }
        }

        // create VolumeRAM for merged data (directly allocates memory)
        if (numVolumes == 2) {
            VolumeAtomic<tgt::Vector2<T> >* vecData = new VolumeAtomic<tgt::Vector2<T> >(dimensions);
            //copy data of the input volumes into the merged data representation
            for (tgt::svec3 i = tgt::svec3(0,0,0); i.z < dimensions.z; ++i.z) {
                //set progress after each slice
                setProgress(static_cast<float>(i.z) / static_cast<float>(dimensions.z));

                for (i.y = 0; i.y < dimensions.y; ++i.y) {
                    for (i.x = 0; i.x < dimensions.x; ++i.x) {
                        tgt::Vector2<T> data = tgt::Vector2<T>(atomicArray[0]->voxel(i), atomicArray[1]->voxel(i));
                        vecData->voxel(i) = data;
                    }
                }
            }
            mergedData = vecData;
        }
        else if (numVolumes == 3) {
            VolumeAtomic<tgt::Vector3<T> >* vecData = new VolumeAtomic<tgt::Vector3<T> >(dimensions);
            //copy data of the input volumes into the merged data representation
            for (tgt::svec3 i = tgt::svec3(0,0,0); i.z < dimensions.z; ++i.z) {
                //set progress after each slice
                setProgress(static_cast<float>(i.z) / static_cast<float>(dimensions.z));

                for (i.y = 0; i.y < dimensions.y; ++i.y) {
                    for (i.x = 0; i.x < dimensions.x; ++i.x) {
                        tgt::Vector3<T> data = tgt::Vector3<T>(atomicArray[0]->voxel(i), atomicArray[1]->voxel(i), atomicArray[2]->voxel(i));
                        vecData->voxel(i) = data;
                    }
                }
            }
            mergedData = vecData;
        }
        else if (numVolumes == 4) {
            VolumeAtomic<tgt::Vector4<T> >* vecData = new VolumeAtomic<tgt::Vector4<T> >(dimensions);
            //copy data of the input volumes into the merged data representation
            for (tgt::svec3 i = tgt::svec3(0,0,0); i.z < dimensions.z; ++i.z) {
                //set progress after each slice
                setProgress(static_cast<float>(i.z) / static_cast<float>(dimensions.z));

                for (i.y = 0; i.y < dimensions.y; ++i.y) {
                    for (i.x = 0; i.x < dimensions.x; ++i.x) {
                        tgt::Vector4<T> data = tgt::Vector4<T>(atomicArray[0]->voxel(i), atomicArray[1]->voxel(i),
                                                                    atomicArray[2]->voxel(i), atomicArray[3]->voxel(i));
                        vecData->voxel(i) = data;
                    }
                }
            }
            mergedData = vecData;
        }
        else {
            setProgress(0.f);
            LERROR("Invalid number of volumes...");
            delete[] atomicArray;
            return 0;
        }

        setProgress(1.f);
        // free memory
        delete[] atomicArray;

        return mergedData;
    }

    virtual void refresh();

private:

    VolumePort volumeInport_;
    VolumePort volumeInport2_;
    VolumePort volumeInport3_;
    VolumePort volumeInport4_;
    VolumePort volumeOutport_;

    BoolProperty autoRefresh_;          ///< if enabled, the output volume is re-computed as soon as the input changes. Else the refresh-button has to be used.
    ButtonProperty refresh_;            ///< if pressed this clears the outport and invalidates the processor so that the output volume is newly computed

    bool refreshNecessary_;             ///< if set to true, process() will compute the new output volume

    static const std::string loggerCat_;
};

} //namespace

#endif
