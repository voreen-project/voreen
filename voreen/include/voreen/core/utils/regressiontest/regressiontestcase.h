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

#ifndef VRN_REGRESSIONTESTCASE_H
#define VRN_REGRESSIONTESTCASE_H

#include "voreen/core/datastructures/datetime.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serialization.h"

#include "tgt/vector.h"

#include <string>
#include <vector>
#include <set>

namespace voreen {

enum RegressionTestResult {
    RegressionTestSuccess,
    RegressionTestFailure,  //< test run was successful, but output data does not match ref data
    RegressionTestIgnored,  //< test run was successful, but its result is ignored
    RegressionTestError,    //< test failed to run (missing ref data dir, seg fault, ...)
    RegressionTestSkipped   //< test was skipped
};

enum FileType {
    ImageFile,
    BinaryFile,
    TextFile,
    FileTypeUnknown
};

struct VRN_CORE_API RegressionTestFile {
    std::string filename_;
    FileType fileType_;
};
VRN_CORE_API std::ostream& operator<<(std::ostream &str, const RegressionTestFile& file);

/// A dataset consists of one or multiple files (e.g. vvd+raw).
struct VRN_CORE_API RegressionTestDataset {
    std::vector<RegressionTestFile> files_;
};


/**
 * Stores the result of a single refFile/outputFile comparison.
 * Note that multiple reference/output files can belong to the same data set (e.g. vvd+raw).
 */
struct VRN_CORE_API FileComparisonResult {
    RegressionTestDataset refDataset_;     ///< reference data set compared
    RegressionTestDataset outputDataset_;  ///< output data set compared
    RegressionTestDataset diffDataset_;    ///< difference created from reference and output data set

    bool match_;                 ///< true, of the reference data set matches the output data set
    std::string report_;         ///< textual comparison report
};

/**
 * Stores configuration options of a test case,
 * which are either loaded from a config file or
 * supplied via command line.
 */
struct VRN_CORE_API RegressionTestCaseConfiguration : public Serializable {
    RegressionTestCaseConfiguration();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    bool enabled_;          ///< if false, test is not run (skipped)
    bool ignored_;          ///< if true, test is run but result is ignored
    int timeout_;           ///< max execution time in seconds

    float pixelDiffTolerance_;
    int maxErrorPixels_;
    int pixelSearchNeighborhood_;

    float voxelDiffTolerance_;
    int maxErrorVoxels_;

    float geometryDiffTolerance_;
};

/// Represents a single test (workspace).
struct VRN_CORE_API RegressionTestCase {
    std::string name_;       ///< name of the test case, usually the test file's basename
    std::string testfile_;   ///< absolute path to the test file
    std::string moduleName_; ///< name of the test file's module
    std::string moduleDir_;  ///< directory name of the test file's module
    std::string group_;      ///< subdirectory of test file in modules/<modulename>/test/
    bool renderingTest_;     ///< if true, test involves OpenGL rendering

    RegressionTestCaseConfiguration config_;
    std::string configfile_; ///< absolute path to config file, if present

    RegressionTestResult result_;
    double duration_;        ///< test runtime in seconds
    std::string errorMsg_;   ///< only set, if result_ == TestError

    std::string call_;       ///< voreentool call used for the test
    int returnCode_;
    std::string htmlLog_;
    std::string consoleLog_;

    /// matching dataset pairs
    std::vector<FileComparisonResult> matchingDatasets_;
    /// corresponding dataset pairs that do not match
    std::vector<FileComparisonResult> mismatchingDatasets_;
    /// reference datasets without corresponding output dataset
    std::vector<RegressionTestDataset> missingRefDatasets_;
    /// output datasets without corresponding output dataset
    std::vector<RegressionTestDataset> unexpectedOutputDatasets_;

    std::string referenceDir_;
    std::string outputDir_;

    /// Returns true, if the test cases have the same module, group and name.
    bool operator==(const RegressionTestCase& t2) const;

    /// Compares by module, then group, then name. Returns true if *this is smaller than t2.
    bool operator<(const RegressionTestCase& t2) const;

};

/// Represents a complete series of RegressionTestCases
struct VRN_CORE_API RegressionTestSuite {
    std::vector<RegressionTestCase> testCases_;
    DateTime date_;

    int numSuccess_;
    int numFailed_;
    int numIgnored_;
    int numError_;
    int numSkipped_;
    double duration_;  //< runtime of all test cases in seconds
};

}

#endif
