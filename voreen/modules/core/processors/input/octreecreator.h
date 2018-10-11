/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_OCTREECREATOR_H
#define VRN_OCTREECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/datastructures/octree/volumeoctreebase.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/stringproperty.h"

namespace voreen {

class OctreeBrickPoolManagerBase;
class VolumeOctreeBase;

struct OctreeCreatorInput {
    bool loadCached;
    std::vector<const VolumeBase*> input;
    OctreeBrickPoolManagerBase* brickPoolManager;
    float homogeneityThreshold;
    size_t brickDim;
    int numThreads;

    OctreeCreatorInput(
        bool loadCached
        , const std::vector<const VolumeBase*>& input
        , OctreeBrickPoolManagerBase* brickPoolManager
        , float homogeneityThreshold
        , size_t brickDim
        , int numThreads
    )
        : loadCached(loadCached)
        , input(input)
        , brickPoolManager(brickPoolManager)
        , homogeneityThreshold(homogeneityThreshold)
        , brickDim(brickDim)
        , numThreads(numThreads)
    {
    }

    OctreeCreatorInput(const OctreeCreatorInput&) = delete;
    OctreeCreatorInput(OctreeCreatorInput&& old)
        : loadCached(old.loadCached)
        , input(old.input)
        , brickPoolManager(std::move(old.brickPoolManager))
        , homogeneityThreshold(old.homogeneityThreshold)
        , brickDim(old.brickDim)
        , numThreads(old.numThreads)
    {
    }
};

struct OctreeCreatorOutput {
    std::unique_ptr<VolumeOctreeBase> octree;
    std::string message;

    OctreeCreatorOutput(
        std::unique_ptr<VolumeOctreeBase> octree
        , const std::string& message
    )
        : octree(std::move(octree))
        , message(message)
    {
    }

    OctreeCreatorOutput(const OctreeCreatorOutput&) = delete;
    OctreeCreatorOutput(OctreeCreatorOutput&& old)
        : octree(std::move(old.octree))
        , message(old.message)
    {
    }
};

class VRN_CORE_API OctreeCreator : public AsyncComputeProcessor<OctreeCreatorInput, OctreeCreatorOutput> {

    friend class VoreenApplication;

public:
    OctreeCreator();
    virtual ~OctreeCreator();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "OctreeCreator";         }
    virtual std::string getCategory() const   { return "Octree";                }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;       }

    virtual bool isReady() const;

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    virtual void setDescriptions() {
        setDescription("Creates an octree volume representation from up to four input volumes. If a multi-channel octree is created, all input channel volumes have to match in dimension, spacing, offset and data type."
                       "<p><strong>Note</strong>: The maximum amount of CPU RAM to be used by the octree as well as the disk storage path of its brick pool are defined globally via application settings.</p>");
    }

    virtual void initialize();
    virtual void deinitialize();

    virtual void adjustPropertiesToInput();

    virtual void setProgressMessage(const std::string& message);

    //void saveOctreeToVVOD();

private:

    std::string getOctreeStoragePath() const;
    std::string getConfigurationHash() const;

    void storeOctreeToCache(const VolumeOctreeBase* octree) const;
    VolumeOctreeBase* restoreOctreeFromCache() const;
    void clearOctree();

    void updatePropertyConfiguration();

    /**
     * Limits the memory used by the octree cache to the specified cache size
     * by deleting octrees in the chronological order of their last access time.
     * Currently called by the VoreenApplication.
     *
     * @param cacheDirectory the cache base directory, usually VoreenApplication::getCachePath()
     * @param cacheSize the maximum cache size in bytes
     * @param keepLatest if true, the latest octree will be kept in the cache even if its size exceeds the limit
     */
    static void limitCacheSize(const std::string& cacheDirectory, const uint64_t maxCacheSize, bool keepLatest);

    /**
     * Deletes the octree part of the cache directory. Currently called by VoreenApplication.
     * @param cacheDirectory the cache base directory, usually VoreenApplication::getCachePath().
     */
    static void deleteCache(const std::string& cacheDirectory);

    VolumePort volumeInport_;
    VolumePort volumeInport2_;
    VolumePort volumeInport3_;
    VolumePort volumeInport4_;
    VolumePort volumeOutport_;
    
    //FileDialogProperty saveOctreeFile_;
    //ButtonProperty saveOctreeButton_;

    IntOptionProperty brickDimensions_;
    IntProperty treeDepth_;
    FloatProperty homogeneityThreshold_;
    BoolProperty useRelativeThreshold_;

    StringOptionProperty brickPoolManager_;
    IntProperty singleBufferMemorySize_;

    IntProperty numThreads_;

    ButtonProperty clearOctree_;

    StringProperty progressMessage_;

    std::string currentConfigurationHash_;

    static const std::string loggerCat_;
};

} //namespace

#endif
