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

#ifndef VRN_FILECOMPARATORS_H
#define VRN_FILECOMPARATORS_H

#include "regressiontestcase.h"

#include "voreen/core/utils/exception.h"
#include "voreen/core/datastructures/geometry/geometry.h"

#include "tgt/vector.h"

#include <string>
#include <vector>
#include <set>

#ifdef VRN_MODULE_DEVIL
#include <IL/il.h>
#include <IL/ilu.h>
#include "modules/devil/devilmodule.h"
#endif

namespace voreen {

class VRN_CORE_API FileComparator {

public:
    virtual void determineFileType(RegressionTestDataset& dataset) const = 0;

    virtual bool supportsFileDiff() const = 0;

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const = 0;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const = 0;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const = 0;

};

//----------------------------------------------------------------------------------------

class VRN_CORE_API BinaryFileComparator : public FileComparator {

public:
    virtual bool supportsFileDiff() const { return false; }

    virtual bool supportsFormat(const RegressionTestDataset& /*dataset*/) const { return true; }

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

};

//----------------------------------------------------------------------------------------

class VRN_CORE_API TextFileComparator : public FileComparator {

public:
    TextFileComparator();

    virtual bool supportsFileDiff() const { return false; }

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const;

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

private:
    std::string loadTextFile(const std::string& filename) const;
    //void saveTextToFile(const std::string& text, const std::string& filename) const;

    std::set<std::string> supportedExtensions_;
};

//----------------------------------------------------------------------------------------

/// Comparator for Voreen Geometry files (.vge)
class VRN_CORE_API VgeFileComparator : public FileComparator {

public:
    VgeFileComparator(float geometryDiffTolerance = 0.f);

    virtual bool supportsFileDiff() const { return false; }

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const;

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

private:
    Geometry* loadGeometry(const std::string& filename) const;

    float geometryDiffTolerance_;
    static const std::string loggerCat_;
};

//----------------------------------------------------------------------------------------

#ifdef VRN_MODULE_DEVIL

class VRN_CORE_API ImageFileComparator : public FileComparator {

public:
    ImageFileComparator(float pixelDiffTolerance = 0.f, int maxErrorPixels = 0,
        int pixelSearchNeighborhood = 0, float diffImageGamma = 1.f, bool diffImageFullAlpha = true);

    virtual bool supportsFileDiff() const { return true; }

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const;

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

    /**
     * Returns the dimensions of the passed image file,
     * or vec2(0) if the dimensions could not be determined.
     */
    static tgt::ivec2 getImageDimensions(const std::string& filename);

private:
    static unsigned int loadImage(const std::string& filename);

    static tgt::ivec2 getImageDimensions(unsigned int image);
    static ILint getImageFormat(unsigned int image);
    static ILint getDataType(unsigned int image);

    float pixelDiffTolerance_;
    int maxErrorPixels_;
    int pixelSearchNeighborhood_;
    float diffImageGamma_;
    bool diffImageFullAlpha_;

    std::set<std::string> supportedExtensions_;

};

#endif // VRN_MODULE_DEVIL

//----------------------------------------------------------------------------------------

class VRN_CORE_API VvdFileComparator : public FileComparator {

public:
    VvdFileComparator(float voxelDiffTolerance = 0.f, int maxErrorVoxels = 0);

    virtual bool supportsFileDiff() const { return false; }

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const;

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

private:
    float voxelDiffTolerance_;
    int maxErrorVoxels_;

    BinaryFileComparator binaryComparator_;

    static const std::string loggerCat_;
};

//----------------------------------------------------------------------------------------

#ifdef VRN_MODULE_HDF5
class VRN_CORE_API HDF5FileComparator : public FileComparator {

public:
    HDF5FileComparator();

    virtual bool supportsFileDiff() const { return false; }

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const;

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

private:

    static const std::string loggerCat_;
};
#endif

//----------------------------------------------------------------------------------------

#ifdef VRN_MODULE_FLOWREEN
class VRN_CORE_API VsdFileComparator : public FileComparator {

public:

    VsdFileComparator();

    virtual bool supportsFileDiff() const { return false; }

    virtual void determineFileType(RegressionTestDataset& dataset) const;

    virtual bool supportsFormat(const RegressionTestDataset& dataset) const;

    virtual bool compare(RegressionTestDataset& refDataset, RegressionTestDataset& outputDataset,
        std::string& report, const RegressionTestCase& testCase) const;

    virtual void generateDiffFile(const RegressionTestDataset& refDataset, const RegressionTestDataset& outputDataset,
        RegressionTestDataset& diffDataset, const RegressionTestCase& testCase) const;

private:

    static const std::string loggerCat_;
};
#endif

} // namespace

#endif
