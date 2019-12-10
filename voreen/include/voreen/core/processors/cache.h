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

#ifndef VRN_CACHING_H
#define VRN_CACHING_H

#include <vector>
#include <string>

#include "voreen/core/processors/processor.h"

namespace voreen {

class VRN_CORE_API Cache {
public:
    Cache(Processor* proc);

    void addInport(Port* inport);
    void addAllInports();

    void addOutport(Port* outport);
    void addAllOutports();

    void addProperty(Property* prop);
    void addAllProperties();

    void initialize();
    bool isInitialized() const { return initialized_; }

    std::string getCurrentCacheDir();
    void clearCache();

    std::string getAllInportHashes();
    std::string getPropertyState();
    std::string getPropertyStateHash();

    bool store();
    bool restore();

    bool restoreOutportsFromDir(const std::string& dir);
    bool storeOutportsToDir(const std::string& dir);
protected:
    std::string getInterfaceString();
    bool stringEqualsFileContent(std::string str, std::string fname);
    bool writeLastAccess(std::string dir, int numAccess);
    bool updateLastAccess(std::string dir);
    bool writePropertyState(std::string dir);

    /// category used in logging
    static const std::string loggerCat_;

private:
    Processor* processor_;
    bool initialized_;

    std::vector<std::string> inports_;
    std::vector<std::string> outports_;

    std::vector<std::string> properties_;
};

/// Cleans up the cache dir, will usually be called by VoreenApplication.
class CacheCleaner {
    /// Small helper class
    struct CacheDir {
        std::string dir_;
        int numReads_;
        uint64_t size_;
        uint64_t age_;

        CacheDir(std::string dir, int numReads, uint64_t size, uint64_t age) : dir_(dir), numReads_(numReads), size_(size), age_(age) {}
    };

public:
    CacheCleaner();
    ~CacheCleaner();
    bool initialize(const std::string& dir);
    void deleteUnused();
    void limitCache(size_t maxMB);
    uint64_t getSize() const;

private:
    bool deleteCacheDir(size_t i);

    std::string dir_;
    std::vector<CacheDir> cacheDirs_;
    uint64_t freedSize_;
    int numDeleted_;

    static const std::string loggerCat_;
};

}  // namespace voreen

#endif
