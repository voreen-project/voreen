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

#include "voreen/core/utils/regressiontest/regressiontestcase.h"
#include "voreen/core/utils/regressiontest/filecomparators.h"
#include "reportgenerators.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/commandlineparser.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/voreenmodule.h"

#include "tgt/init.h"
#include "tgt/logmanager.h"
#include "tgt/filesystem.h"
#include "tgt/stopwatch.h"

#include <string>
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <algorithm>
#include <ctime>

const std::string APP_BINARY  = "regressiontest";
const std::string APP_NAME    = "Regression-Test";
const std::string APP_DESC    = "Regression test app based on Voreen workspaces.";

const std::string loggerCat_  = "regressiontest.main";

#ifdef WIN32
const std::string VOREENTOOL_BINARY_NAME = "voreentool.exe";
#else
const std::string VOREENTOOL_BINARY_NAME = "voreentool";
#endif

using namespace voreen;
using namespace tgt;

void exitSuccess(const std::string& msg);
void exitError(const std::string& errorMsg, bool printUsage = false);

/**
 * Collects all test cases from the passed paths, reads their config files if present,
 * and returns them as TestSuite.
 *
 * @param testPaths test path to search for test cases
 * @param timeout max execution time of each test case in seconds (written to test case configuration)
 * @param testSuite Test suite containing the found test cases
 * @param modules all modules corresponding to the passed testPaths
 *
 * @throw VoreenException if any of the passed test paths is invalid
 */
void collectTestCases(std::vector<std::string> testPaths, int timeout,
    RegressionTestSuite& testSuite, std::vector<VoreenModule*>& modules);

/**
 * Returns a mapping from modules to a list of their Processors that are to be excluded from the coverage analysis.
 * The ignore lists are retrieved from the coverage-ignore.xml files in the modules' test dirs.
 */
std::map<VoreenModule*, std::vector<Processor*> > collectCoverageIgnoreLists(const std::vector<VoreenModule*>& modules);

/**
 * Runs all test cases of the passed TestSuite and stores the results in it.
 *
 * @param testSuite the TestSuite to execute
 * @param testdataDir absolute path to test data base directory (checked)
 * @param reportDir absolute path to directory where the output data should be copied to (checked)
 * @param keepReferenceData if true, the reference data will be copied to the report dir
 * @param voreentoolPath (checked)
 * @param fileComparators fileComparators to use for comparing output/reference files
 * @param skipRenderingTests if true, only data processing (non OpenGL) tests are executed
 * @param redirectStdOut if true, console output of the testcase will be redirected to file
 */
void runTestSuite(RegressionTestSuite& testSuite, const std::string& testdataDir,
    const std::string& reportDir, bool keepReferenceData,
    const std::string& voreentoolPath,
    const std::vector<FileComparator*>& fileComparators,
    bool skipRenderingTests, bool redirectStdOut, bool glclsharing);

/**
 * Runs a single test case.
 *
 * @param testCase test case to execute
 * @param testdataDir absolute path to test data base directory (checked)
 * @param reportDir absolute path to directory where the output data should be copied to (checked)
 * @param keepReferenceData if true, the reference data will be copied to the report dir
 * @param voreentoolPath (checked)
 * @param fileComparators fileComparators to use for comparing output/reference files
 * @param redirectStdOut if true, console output of the testcase will be redirected to file
 */
void runTestCase(RegressionTestCase& testCase, const std::string& testdataDir,
    const std::string& reportDir, bool keepReferenceData,
    const std::string& voreentoolPath,
    const std::vector<FileComparator*>& fileComparators,
    bool redirectStdOut, bool glclsharing);

/**
 * Determines the test coverage of the passed processors with regard to the passed TestSuite.
 *
 * @return map from each of the passed processors to a vector of test cases that contain this processor
 */
std::map<Processor*, std::vector<RegressionTestCase> > getCoverageMap(const RegressionTestSuite& testSuite, const std::vector<Processor*> processors);


/**
 * Executes a test workspace using the passed voreentool binary and optional working directory/log file paths.
 *
 * @note the passed TestCase must contain a valid testfile_ member
 */
void executeWorkspace(const std::string& binary, RegressionTestCase& testCase, const std::string& workDir = "",
    bool useCaching = false, int timeout = 0, const std::string& htmlLog = "", const std::string& consoleLog = "", bool glclsharing=true);

/**
 * Does a pair-wise comparison of the datasets found in the output and reference dirs of the passed TestCase,
 * using the passed FileComparators.
 *
 * @param testCase testcase whose output and reference dirs are to be compared.
 *  The results are also stored in the passed testcase.
 * @param ignoreFiles list of files to ignore
 * @param fileComparators fileComparators to use
 *
 * @return true, if the output data matches the reference data
 */
bool compareResults(RegressionTestCase& testCase, const std::set<std::string>& ignoreFiles,
    const std::vector<FileComparator*>& fileComparators);

/// Uses the passed FileComparators to determine the type of the files of the passed dataset.
void determineFileType(RegressionTestDataset& dataset, const std::vector<FileComparator*>& fileComparators);

/// Moves all files from the passed srcDir to the passed destDir.
void moveDirectoryContents(const std::string& srcDir, const std::string& destDir);

/// Copies all files from the passed srcDir to the passed destDir.
void copyDirectoryContents(const std::string& srcDir, const std::string& destDir);

int main(int argc, char* argv[]) {

    VoreenApplication vapp(APP_BINARY, APP_NAME, APP_DESC, argc, argv,
        VoreenApplication::ApplicationFeatures(VoreenApplication::APP_ALL &~ VoreenApplication::APP_WIDGETS));

    // prepare command line parser
    CommandLineParser* cmdParser = vapp.getCommandLineParser();
    tgtAssert(cmdParser, "no CommandLineParser");
    int cmdStyle = po::command_line_style::default_style ^ po::command_line_style::allow_guessing;
    cmdParser->setCommandLineStyle(static_cast<po::command_line_style::style_t>(cmdStyle));

    std::vector<std::string> testPaths;
    cmdParser->addMultiOption<std::string>("testcases,t", testPaths, CommandLineParser::RequiredOption,
        "List of test paths containing the test cases to execute. Each test path may either specify a concrete test file "
        "or a directory and has the form <module>[/<testCasePath>] where <testCasePath> is relative to the module's test directory.\n"
        "Examples: \n"
        "- 'base' runs all test cases of the base module\n"
        "- 'base/render' runs all test cases in modules/base/test/render/ \n"
        "- 'base/render/background.vws' runs the single test modules/base/test/render/background.vws \n\n"
        "Note: use 'all' for running all tests\n");

    bool testCoverageMode = false;
    cmdParser->addFlagOption("testcoverage", testCoverageMode, CommandLineParser::MainOption,
        "Determine test coverage instead of actually running the test cases: For each processor "
        "of the specified modules the test cases containing that processor are detected.");

    std::string testdataDir;
    cmdParser->addOption<std::string>("datadir,d", testdataDir, CommandLineParser::RequiredOption,
        "Test data base directory (contains sub directories 'input' and 'reference').");

    std::string voreentoolPath;
    if (FileSystem::fileExists(FileSystem::currentDirectory() + "/" + VOREENTOOL_BINARY_NAME))
        voreentoolPath = FileSystem::cleanupPath(FileSystem::currentDirectory() + "/" + VOREENTOOL_BINARY_NAME);
    else {
#ifdef WIN32
    #ifdef NDEBUG
        voreentoolPath = "$VRN_BASE/bin/Release/" + VOREENTOOL_BINARY_NAME;
    #else
        voreentoolPath = "$VRN_BASE/bin/Debug/" + VOREENTOOL_BINARY_NAME;
    #endif
#else
        voreentoolPath = "$VRN_BASE/bin/" + VOREENTOOL_BINARY_NAME;
#endif
    }
    cmdParser->addOption<std::string>("voreentool", voreentoolPath, CommandLineParser::MainOption,
        "Path to the voreentool binary.",
        FileSystem::cleanupPath(voreentoolPath));

    std::string reportDir;
    cmdParser->addOption<std::string>("reportdir", reportDir, CommandLineParser::MainOption,
        "Directory where the tests' output data and optionally the reference files will be stored. "
        "If not specified, the output data will be discarded.");

    bool keepReferenceData = false;
    cmdParser->addFlagOption("keepReferenceData", keepReferenceData, CommandLineParser::MainOption,
        "Copy reference data to report directory. This may be useful for archiving test results.");

    bool skipRenderingTests = false;
    cmdParser->addFlagOption("skipRenderingTests", skipRenderingTests, CommandLineParser::MainOption,
        "Run data processing tests only.");

    int timeout = 0;
    cmdParser->addOption("timeout", timeout, CommandLineParser::MainOption,
        "Maximum execution time of a single test case in seconds (unix only).");

    float pixelDiffTolerance;
    cmdParser->addOption("pixelDiffTolerance", pixelDiffTolerance, CommandLineParser::MainOption,
        "Maximum allowed per-pixel color distance in image comparisons (RGBA color space). Range: [0.0;1.0]",
        0.02f, "0.02");

    int maxErrorPixels;
    cmdParser->addOption("maxErrorPixels", maxErrorPixels, CommandLineParser::MainOption,
        "Maximum allowed number of pixels with a difference above the tolerance.",
        50);

    int pixelSearchNeighborhood;
    cmdParser->addOption("pixelSearchNeighborhood", pixelSearchNeighborhood, CommandLineParser::MainOption,
        "Radius of the neighborhood around a pixel to search for a match in the reference image.",
        1);

    float diffImageGamma;
    cmdParser->addOption("diffImageGamma", diffImageGamma, CommandLineParser::MainOption,
        "Gamma value to apply to difference images (must be non-negative).",
        0.25f, "0.25");

    bool diffImageFullAlpha;
    cmdParser->addOption("diffImageFullAlpha", diffImageFullAlpha, CommandLineParser::MainOption,
        "Set alpha value of differing pixels to 1.0 in difference images.",
        true, "true");

    float voxelDiffTolerance;
    cmdParser->addOption("voxelDiffTolerance", voxelDiffTolerance, CommandLineParser::MainOption,
        "Maximum allowed per-voxel difference in volume comparisons. Range: [0.0;1.0]",
        1e-06f, "1e-06");

    int maxErrorVoxels;
    cmdParser->addOption("maxErrorVoxels", maxErrorVoxels, CommandLineParser::MainOption,
        "Maximum allowed number of voxels with a difference above the tolerance.",
        0);

    float geometryDiffTolerance;
    cmdParser->addOption("geometryDiffTolerance", geometryDiffTolerance, CommandLineParser::MainOption,
        "Maximum allowed vertex distance in geometry comparisons (must be non-negative)",
        1e-06f, "1e-06");

    std::string htmlFile;
    cmdParser->addOption<std::string>("htmlReport", htmlFile, CommandLineParser::MainOption,
        "Write HTML report to this file.");

    bool splitHtmlReport;
    cmdParser->addOption<bool>("splitHtmlReport", splitHtmlReport, CommandLineParser::AdditionalOption,
        "Split HTML report into multiple files, one for each module.", false, "false");

    std::string junitFile;
    cmdParser->addOption<std::string>("junitReport", junitFile, CommandLineParser::MainOption,
        "Write JUnit XML report to this file.");

    bool redirectStdOut;
    cmdParser->addOption("redirectStdOut", redirectStdOut, CommandLineParser::MainOption,
        "Redirect voreentool console output to separate log file for each test case.",
        true, "true");

    bool noglclsharing = false;
#ifdef VRN_MODULE_OPENCL
    cmdParser->addFlagOption("noglclsharing", noglclsharing, CommandLineParser::MainOption,
        "Disable OpenGL/OpenCL sharing for performance penalty, but improved compatibility in OpenGL tests.");
#endif

    // init application
    try {
        vapp.initialize();
    }
    catch (VoreenException& e) {
        exitError("Failed to initialize application: " + std::string(e.what()));
    }
    std::string vrnBasePath = VoreenApplication::app()->getBasePath();

    //
    // check parameters
    //
    tgtAssert(!testPaths.empty(), "--testcases parameter is empty"); //< required parameter

    if (!testCoverageMode) {
        tgtAssert(!voreentoolPath.empty(), "voreentool path is empty"); //< has default parameter
        voreentoolPath = FileSystem::cleanupPath(strReplaceAll(voreentoolPath, "$VRN_BASE", vrnBasePath));

        if (testdataDir.empty()) {
            exitError("Please specify test data directory (--datadir)");
        }
        testdataDir = strReplaceAll(testdataDir, "$VRN_BASE", vrnBasePath);
        testdataDir = FileSystem::absolutePath(testdataDir);
        if (!FileSystem::dirExists(testdataDir))
            exitError("Test data directory '" + testdataDir + "' does not exist");
        else if (!FileSystem::dirExists(testdataDir + "/reference"))
            exitError("Test data directory '" + testdataDir + "' does not contain 'reference' sub directory");
        else
            LINFO("Test data directory: " << testdataDir);

        // check report dir and clear subdirectories 'output/' and 'reference/'
        if (!reportDir.empty()) {
            reportDir = FileSystem::absolutePath(strReplaceAll(reportDir, "$VRN_BASE", vrnBasePath));
            if (FileSystem::dirExists(reportDir)) {
                LINFO("Report directory: " << reportDir);
                std::string reportOutputDir = FileSystem::cleanupPath(reportDir + "/" + "output");
                if (FileSystem::dirExists(reportOutputDir)) {
                    LDEBUG("Deleting report output directory: " << reportOutputDir);
                    if (!FileSystem::deleteDirectoryRecursive(reportOutputDir))
                        LWARNING("Failed to clear report output directory: " << reportOutputDir);
                }
                std::string reportRefDir = FileSystem::cleanupPath(reportDir + "/" + "reference");
                if (FileSystem::dirExists(reportRefDir)) {
                    LDEBUG("Deleting report reference directory: " << reportRefDir);
                    if (!FileSystem::deleteDirectoryRecursive(reportRefDir))
                        LWARNING("Failed to clear report reference directory: " << reportRefDir);
                }
            }
            else {
                LINFO("Creating report directory: " << reportDir);
                if (!FileSystem::createDirectoryRecursive(reportDir)) {
                    LERROR("Failed to create report directory: " << reportDir);
                    reportDir = "";
                }
            }
        }

    #ifdef WIN32
        if (timeout > 0) {
            LWARNING("Test timeout not supported on Windows");
            timeout = 0;
        }
    #else
        if (timeout < 0)
            exitError("Test timeout must be non-negative: " + itos(timeout));
        else if (timeout == 0)
            LINFO("Running tests without timeout");
        else
            LINFO("Test timeout: " << timeout << " seconds");
    #endif

        if (pixelDiffTolerance < 0.f || pixelDiffTolerance > 1.f)
            exitError("pixelDiffTolerance must be within range [0.0,1.0]");
        if (maxErrorPixels < 0)
            exitError("maxErrorPixels must be non-negative");
        if (pixelSearchNeighborhood < 0)
            exitError("pixelSearchNeighborhood must be non-negative");
        if (diffImageGamma < 0.f)
            exitError("diffImageGamma must be non-negative");

        if (voxelDiffTolerance < 0.f || voxelDiffTolerance > 1.f)
            exitError("voxelDiffTolerance must be within range [0.0,1.0]");
        if (maxErrorVoxels < 0)
            exitError("maxErrorVoxels must be non-negative");

        if (geometryDiffTolerance < 0.f)
            exitError("geometryDiffTolerance must be non-negative");

        if (!htmlFile.empty()) {
            htmlFile = FileSystem::absolutePath(strReplaceAll(htmlFile, "$VRN_BASE", vrnBasePath));
            std::string htmlDir = FileSystem::dirName(htmlFile);
            if (FileSystem::dirExists(htmlDir))
                LINFO("HTML report: " << htmlFile);
            else
                exitError("Directory of HTML report file does not exist: " + htmlDir);
        }
    } // !testCoverageMode

    if (!junitFile.empty()) {
        junitFile = FileSystem::absolutePath(strReplaceAll(junitFile, "$VRN_BASE", vrnBasePath));
        std::string junitDir = FileSystem::dirName(junitFile);
        if (FileSystem::dirExists(junitDir))
            LINFO("JUnit report: " << junitFile);
        else
            exitError("Directory of JUnit report file does not exist: " + junitDir);
    }


    //
    // Setup FileComparators
    //
    std::vector<FileComparator*> fileComparators;
#ifdef VRN_MODULE_DEVIL
    fileComparators.push_back(new ImageFileComparator(pixelDiffTolerance, maxErrorPixels,
        pixelSearchNeighborhood, diffImageGamma, diffImageFullAlpha));
    LINFO("Using ImageFileComparator with "
        << "pixelDiffTolerance: " << pixelDiffTolerance
        << ", maxErrorPixels: " << maxErrorPixels
        << ", pixelSearchNeighorhood: " << pixelSearchNeighborhood
        << ", diffImageGamma: " << diffImageGamma
        << ", diffImageFullAlpha: " << diffImageFullAlpha);
#else
    LWARNING("DevIL module not enabled: images are compared binary");
#endif

#ifdef VRN_MODULE_HDF5
    fileComparators.push_back(new HDF5FileComparator());
    LINFO("Using HDF5FileComparator");
#else
    LWARNING("HDF5 module not enabled: HDF5 files are compared binarily");
#endif
#ifdef VRN_MODULE_FLOWANALYSIS
    fileComparators.push_back(new VsdFileComparator());
    LINFO("Using VsdFileComparator");
#else
    LWARNING("FlowAnalysis module not enabled: Vsd files are compared binarily");
#endif

    fileComparators.push_back(new VvdFileComparator(voxelDiffTolerance, maxErrorVoxels));
    LINFO("Using VvdFileComparator with "
        << "voxelDiffTolerance: " << voxelDiffTolerance
        << ", maxErrorVoxels: " << maxErrorVoxels );

    // needs to inserted before TextFileComparator, since .vge are also text files
    fileComparators.push_back(new VgeFileComparator(geometryDiffTolerance));
    LINFO("Using VgeFileComparator with geometryDiffTolerance: " << geometryDiffTolerance);

    fileComparators.push_back(new TextFileComparator());
    fileComparators.push_back(new BinaryFileComparator()); //< needs to be last, since it accepts all files

    //
    // Collect test cases and tested modules
    //
    RegressionTestSuite testSuite;
    std::vector<VoreenModule*> modules;
    try {
        collectTestCases(testPaths, timeout, testSuite, modules);
    }
    catch (VoreenException& e) {
        exitError(e.what());
    }
    LINFO("Found " << testSuite.testCases_.size() << " test cases.\n");

    //
    // Run tests or determine test coverage
    //
    if (!testCoverageMode) {
        // run tests
        runTestSuite(testSuite, testdataDir, reportDir, keepReferenceData,
            voreentoolPath, fileComparators, skipRenderingTests, redirectStdOut, !noglclsharing);

        // generate file reports
        if (!htmlFile.empty()) {
            try {
                HTMLReportGenerator().generateTestResultReport(htmlFile, testSuite, splitHtmlReport);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
            }
        }
        if (!junitFile.empty()) {
            try {
                JUnitReportGenerator().generateTestResultReport(junitFile, testSuite);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
            }
        }

        // print console report
        LINFO("Success: " << testSuite.numSuccess_ << ", Failure: " << testSuite.numFailed_ << ", Ignored: " << testSuite.numIgnored_ <<
            ", Error: " << testSuite.numError_ << ", Skipped: " << testSuite.numSkipped_ << " (" << testSuite.duration_ << " sec)");

        std::string message;
        message += "---------------------------\n";
        message += "Num Tests: " + itos(testSuite.testCases_.size());
        message += " (" + dtos(testSuite.duration_) + " sec)\n";
        message += "- Success: " + itos(testSuite.numSuccess_) + "\n";
        message += "- Failure: " + itos(testSuite.numFailed_) + "\n";
        message += "- Ignored: " + itos(testSuite.numIgnored_) + "\n";
        message += "- Error:   " + itos(testSuite.numError_) + "\n";
        message += "- Skipped: " + itos(testSuite.numSkipped_) + "\n";
        message += "---------------------------";
        if(testSuite.numFailed_ > 0 || testSuite.numError_ > 0) {
            exitError(message);
        } else {
            exitSuccess(message);
        }
    }
    else { // test coverage mode

        // collect processors of tested modules
        std::vector<Processor*> processors;
        for (size_t i=0; i<modules.size(); i++) {
            VoreenModule* module = modules.at(i);
            std::vector<const Processor*> modProcessors = module->getRegisteredProcessors();
            for (size_t j=0; j<modProcessors.size(); j++)
                processors.push_back(const_cast<Processor*>(modProcessors[j]));
            //processors.insert(processors.end(), module->getProcessors().begin(), module->getProcessors().end());
        }

        // collect coverage ignore lists
        std::map<VoreenModule*, std::vector<Processor*> > ignoreLists = collectCoverageIgnoreLists(modules);

        // determine coverage map for module processors
        std::map<Processor*, std::vector<RegressionTestCase> > coverageMap = getCoverageMap(testSuite, processors);
        tgtAssert(processors.size() == coverageMap.size(), "coverage map size does not match number of processors");

        // generate files reports
        if (!htmlFile.empty()) {
            LWARNING("HTML report not supported in test coverage mode");
        }
        if (!junitFile.empty()) {
            try {
                JUnitReportGenerator().generateTestCoverageReport(junitFile, modules, coverageMap, ignoreLists);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
            }
        }

        exitSuccess(ConsoleReporter().generateTestCoverageReport(modules, coverageMap, ignoreLists));
    }

    tgtAssert(false, "should never get here");
    return 0;
}


//----------------------------------------------------------------------------------------

void exitSuccess(const std::string& msg) {
    //    LINFO(msg);
    VoreenApplication::app()->deinitialize();

    std::cout << std::endl << msg << std::endl << std::endl;
    exit(EXIT_SUCCESS);
}

void exitError(const std::string& errorMsg, bool printUsage) {
    LFATAL("ERROR: " << errorMsg);

    std::string usageString =
        VoreenApplication::app()->getCommandLineParser()->getUsageString(CommandLineParser::AllTypes, false);
    VoreenApplication::app()->deinitialize();

    std::cerr << std::endl << "ERROR: " << errorMsg << std::endl << std::endl;

    if (printUsage)
        std::cerr << usageString << std::endl << std::endl;

    exit(EXIT_FAILURE);
}

void collectTestCases(std::vector<std::string> testPaths, int timeout,
    RegressionTestSuite& testSuite, std::vector<VoreenModule*>& modules) {
    tgtAssert(timeout >= 0, "timeout is negative");
    VoreenApplication* vapp = VoreenApplication::app();
    tgtAssert(vapp, "VoreenApplication not instantiated");
    const std::string vrnBasePath = vapp->getBasePath();

    // replace pseudo-path "all" by all module dirs which have a test/ subdir
    if (std::find(testPaths.begin(), testPaths.end(), "all") != testPaths.end()) {
        // collect directory names of all available modules
        std::set<std::string> moduleDirs;
        for (size_t i=0; i<vapp->getModules().size(); i++)
            moduleDirs.insert(vapp->getModules().at(i)->getDirName());

        testPaths.clear();
        testPaths.insert(testPaths.begin(), moduleDirs.begin(), moduleDirs.end());
    }

    // iterate over test paths and collect test cases
    std::vector<RegressionTestCase> testCases;
    for (size_t i=0; i<testPaths.size(); i++) {
        const std::string testPath = testPaths.at(i);

        std::vector<std::string> testPathParts = FileSystem::splitPath(testPath);
        tgtAssert(testPathParts.size() > 0, "path split is empty");

        // extract test case module
        std::string moduleDir = testPathParts.at(0);
        VoreenModule* module = vapp->getModule(moduleDir);
        if (!module)
            throw VoreenException("Invalid test path '" + testPath + "': module '" + moduleDir + "' not found");
        else if (std::find(modules.begin(), modules.end(), module) == modules.end())
            modules.push_back(module);
        std::string moduleTestDir = module->getModulePath("test");
        testPathParts.erase(testPathParts.begin());
        if (testPathParts.empty() && !FileSystem::dirExists(moduleTestDir)) {
            LDEBUG("Module '" << moduleDir << "' does not have a test/ subdirectory");
            continue;
        }

        // create absolute test path and check whether it exists and whether it references a file or a directory
        std::string testSubPath = strJoin(testPathParts, "/");
        std::string testPathAbs = FileSystem::cleanupPath(moduleTestDir + "/" + testSubPath);
        std::string fileExtension = FileSystem::fileExtension(testPathAbs);
        bool isFile = fileExtension != "";
        if (isFile && !FileSystem::fileExists(testPathAbs))
            throw VoreenException("Invalid test path '" + testPath + "': test file '" + testPathAbs + "' does not exist");
        else if (isFile && fileExtension != "vws" && fileExtension != "vnw")
            throw VoreenException("Invalid test path '" + testPath + "': " + testPathAbs + " is not a valid test file");
        else if (!isFile) {
            if (FileSystem::dirExists(testPathAbs))
                ;
            else if (FileSystem::fileExists(testPathAbs + ".vws")) {
                testPathAbs += ".vws";
                isFile = true;
            }
            else if (FileSystem::fileExists(testPathAbs + ".vnw")) {
                testPathAbs += ".vnw";
                isFile = true;
            }
            else
                throw VoreenException("Invalid test path '" + testPath + "': directory '" + testPathAbs + "' does not exist");
        }

        if (isFile) {
            // test path references a single test file
            tgtAssert(testPathParts.size() > 0, "test path parts empty");
            RegressionTestCase testCase;
            testCase.name_ = FileSystem::baseName(testPathAbs);
            testCase.testfile_ = testPathAbs;
            testCase.moduleName_ = module->getID();
            testCase.moduleDir_ = moduleDir;
            std::vector<std::string> groupElems = FileSystem::splitPath(FileSystem::dirName(testSubPath));
            testCase.group_ = strJoin(groupElems, ".");
            testCase.renderingTest_ = (groupElems.size() > 0 && groupElems.at(0) == "render");
            testCases.push_back(testCase);
        }
        else {
            // test path references a directory => collect all test files in that directory
            std::vector<std::string> testFiles = FileSystem::listFilesRecursive(testPathAbs, true);
            for (size_t j=0; j<testFiles.size(); j++) {
                std::string testFileAbs = FileSystem::cleanupPath(testPathAbs + "/" + testFiles.at(j));
                std::string testFileDir = FileSystem::dirName(testFileAbs);
                std::string fileExtension = FileSystem::fileExtension(testFiles.at(j));
                if (fileExtension != "vws" && fileExtension != "vnw")
                    continue;

                RegressionTestCase testCase;
                testCase.name_ = FileSystem::baseName(testFileAbs);
                testCase.testfile_ = testFileAbs;
                testCase.moduleName_ = module->getID();
                testCase.moduleDir_ = moduleDir;
                std::vector<std::string> groupElems = FileSystem::splitPath(FileSystem::relativePath(testFileDir, moduleTestDir));
                testCase.group_ = strJoin(groupElems, ".");
                testCase.renderingTest_ = (groupElems.size() > 0 && groupElems.at(0) == "render");
                testCases.push_back(testCase);
            }
        }
    }

    // sort test cases and remove duplicates
    std::list<RegressionTestCase> testCasesList;
    testCasesList.insert(testCasesList.begin(), testCases.begin(), testCases.end());
    testCasesList.sort();
    testCasesList.unique();
    testCases.clear();
    testCases.insert(testCases.begin(), testCasesList.begin(), testCasesList.end());

    // configure test cases
    for (size_t i=0; i<testCases.size(); i++) {
        RegressionTestCase& testCase = testCases.at(i);

        testCase.config_.timeout_ = timeout;

        // read config file, if present
        std::string cfgFileAbs = FileSystem::fullBaseName(testCase.testfile_) + ".cfg";
        if (FileSystem::fileExists(cfgFileAbs)) {
            testCase.configfile_ = cfgFileAbs;

            LDEBUG("Reading test config file: " << FileSystem::relativePath(testCase.configfile_, vrnBasePath));
            std::ifstream cfgFile(testCase.configfile_.c_str(), std::ios_base::in);
            if (cfgFile.good()) {
                XmlDeserializer deserializer;
                try {
                    deserializer.read(cfgFile);
                    Deserializer(deserializer).deserialize("TestCaseConfiguration", testCase.config_);
                }
                catch (SerializationException& e) {
                    LWARNING("Failed to deserialize config file '" << testCase.configfile_ << "': " << e.what());
                }
                cfgFile.close();
            }
            else {
                LWARNING("Failed to open config file for reading: " << testCase.configfile_);
            }
        }
    }

    // create test suite
    testSuite.testCases_ = testCases;
    testSuite.date_ = DateTime(time(0));
    testSuite.duration_ = 0.0;
    testSuite.numSuccess_ = 0;
    testSuite.numFailed_ = 0;
    testSuite.numIgnored_ = 0;
    testSuite.numError_ = 0;
    testSuite.numSkipped_ = 0;
}

std::map<VoreenModule*, std::vector<Processor*> > collectCoverageIgnoreLists(const std::vector<VoreenModule*>& modules) {

    tgtAssert(VoreenApplication::app(), "VoreenApplication not instantiated");
    std::string vrnBasePath = VoreenApplication::app()->getBasePath();

    std::map<VoreenModule*, std::vector<Processor*> > ignoreMap;

    for (std::vector<VoreenModule*>::const_iterator it=modules.begin(); it != modules.end(); ++it) {
        VoreenModule* module = *it;
        std::vector<std::string> ignoreProcessorNames;

        // read processor names from file, if present
        std::string ignoreFilePath = module->getModulePath("test/coverage-ignore.xml");
        if (FileSystem::fileExists(ignoreFilePath)) {
            LINFO("Reading coverage ignore file: " << FileSystem::relativePath(ignoreFilePath, vrnBasePath));
            std::ifstream ignoreFile(ignoreFilePath.c_str(), std::ios_base::in);
            if (ignoreFile.good()) {
                XmlDeserializer deserializer;
                try {
                    deserializer.read(ignoreFile);
                    Deserializer(deserializer).deserialize("Processors", ignoreProcessorNames);
                }
                catch (SerializationException& e) {
                    LWARNING("Failed to deserialize coverage ignore file '" << ignoreFilePath << "': " << e.what());
                }
                ignoreFile.close();

                if (ignoreProcessorNames.empty())
                    LWARNING("No processors found in coverage ignore file: " << ignoreFilePath);
            }
            else {
                LWARNING("Failed to open coverage ignore file for reading: " << ignoreFilePath);
            }
        }

        // determine processors from read names
        std::vector<const Processor*> moduleProcessors = module->getRegisteredProcessors();
        std::vector<Processor*> ignoreProcessors;
        for (size_t i=0; i<ignoreProcessorNames.size(); i++) {
            std::string procName = ignoreProcessorNames.at(i);
            bool found = false;
            for (size_t j=0; j<moduleProcessors.size() && !found; j++) {
                if (moduleProcessors.at(j)->getClassName() == procName) {
                    ignoreProcessors.push_back(const_cast<Processor*>(moduleProcessors.at(j)));
                    found = true;
                }
            }
            if (!found) {
                LWARNING("Coverage ignore: processor '" << procName << "' not found in module '" << module->getDirName() << "'");
            }
        }
        ignoreMap.insert(std::make_pair(module, ignoreProcessors));
    }

    return ignoreMap;
}

void runTestSuite(RegressionTestSuite& testSuite, const std::string& testdataDir,
        const std::string& reportDir, bool keepReferenceData,
        const std::string& voreentoolPath,
        const std::vector<FileComparator*>& fileComparators,
        bool skipRenderingTests, bool redirectStdOut, bool glclsharing)
{

    tgtAssert(!voreentoolPath.empty(), "voreentool path is empty");
    tgtAssert(VoreenApplication::app(), "VoreenApplication not instantiated");
    tgtAssert(FileSystem::dirExists(testdataDir), "testdata dir does not exist");  //< checked by caller
    if (reportDir != "")
        tgtAssert(FileSystem::dirExists(reportDir), "report dir does not exist");  //< checked/created by caller

    std::string vrnBasePath = VoreenApplication::app()->getBasePath();

    // check/create output dir
    std::string outputdir = FileSystem::cleanupPath(testdataDir + "/output");
    if (!FileSystem::dirExists(outputdir)) {
        LINFO("Creating output directory: " << outputdir);
        FileSystem::createDirectory(outputdir);
    }
    else
        LDEBUG("Output dir: " << outputdir);

    testSuite.date_ = DateTime(time(0));

    // run test cases
    for (size_t i=0; i<testSuite.testCases_.size(); i++) {
        RegressionTestCase& testCase = testSuite.testCases_.at(i);

        if (testCase.renderingTest_ && skipRenderingTests)
            testCase.config_.enabled_ = false;

        runTestCase(testCase, testdataDir, reportDir, keepReferenceData,
            voreentoolPath, fileComparators, redirectStdOut, glclsharing);
    }

    // collect statistics
    testSuite.numSuccess_ = 0;
    testSuite.numFailed_ = 0;
    testSuite.numIgnored_ = 0;
    testSuite.numError_ = 0;
    testSuite.numSkipped_ = 0;
    testSuite.duration_ = 0.0;
    for (size_t i=0; i<testSuite.testCases_.size(); i++) {
        const RegressionTestCase& testCase = testSuite.testCases_.at(i);
        testSuite.duration_ += testCase.duration_;
        if (testCase.result_ == RegressionTestSuccess)
            testSuite.numSuccess_++;
        else if (testCase.result_ == RegressionTestFailure)
            testSuite.numFailed_++;
        else if (testCase.result_ == RegressionTestIgnored)
            testSuite.numIgnored_++;
        else if (testCase.result_ == RegressionTestError)
            testSuite.numError_++;
        else if (testCase.result_ == RegressionTestSkipped)
            testSuite.numSkipped_++;
    }
    tgtAssert(static_cast<int>(testSuite.testCases_.size()) ==
        testSuite.numSuccess_ + testSuite.numFailed_ + testSuite.numIgnored_ + testSuite.numError_ + testSuite.numSkipped_,
        "invalid statistics");
}

void runTestCase(RegressionTestCase& testCase, const std::string& testdataDir,
        const std::string& reportDir, bool keepReferenceData,
        const std::string& voreentoolPath,
        const std::vector<FileComparator*>& fileComparators,
        bool redirectStdOut, bool glclsharing) {

    tgtAssert(VoreenApplication::app(), "VoreenApplication not instantiated");
    std::string vrnBasePath = VoreenApplication::app()->getBasePath();

    tgtAssert(FileSystem::fileExists(testCase.testfile_), "test file does not exist");
    tgtAssert(testCase.moduleDir_ != "", "test module is empty");

    testCase.result_ = RegressionTestError;
    testCase.duration_ = 0.0;

    // determine test case directories
    std::string baseName = FileSystem::baseName(testCase.testfile_);
    std::string testDir = FileSystem::dirName(testCase.testfile_);
    std::string moduleTestDir = VoreenApplication::app()->getModulePath(testCase.moduleDir_) + "/test";
    std::string testSubDir = testCase.moduleDir_ + "/";
    if (FileSystem::relativePath(testDir, moduleTestDir) != "")
        testSubDir += FileSystem::relativePath(testDir, moduleTestDir) + "/";
    testSubDir += testCase.name_;

    // rendering or data processing test?
    if (testCase.config_.enabled_) {
        if (testCase.renderingTest_)
            LINFO("Running rendering test: " << FileSystem::relativePath(testCase.testfile_, vrnBasePath));
        else
            LINFO("Running data processing test: " << FileSystem::relativePath(testCase.testfile_, vrnBasePath));
    }
    else {
        if (testCase.renderingTest_)
            LINFO("SKIPPED rendering test: " << FileSystem::relativePath(testCase.testfile_, vrnBasePath) << "\n");
        else
            LINFO("SKIPPED data processing test: " << FileSystem::relativePath(testCase.testfile_, vrnBasePath) << "\n");
        testCase.result_ = RegressionTestSkipped;
        return;
    }

    // check reference dir
    std::string refDir = FileSystem::cleanupPath(testdataDir + "/reference/" + testSubDir);
    if (!FileSystem::dirExists(refDir)) {
        testCase.errorMsg_ = "Reference data directory '" + refDir + "' does not exist";
        LERROR("Test ERROR: " << testCase.errorMsg_  << "\n");
        return;
    }
    testCase.referenceDir_ = refDir;

    // check/create/clean output dir
    std::string outputDir = FileSystem::cleanupPath(testdataDir + "/output");
    if (!FileSystem::dirExists(outputDir)) {
        LDEBUG("Creating output directory: " << outputDir);
        if (!FileSystem::createDirectory(outputDir)) {
            testCase.errorMsg_ = "Failed to create output directory: " + outputDir;
            LERROR("Test ERROR: " << testCase.errorMsg_  << "\n");
            return;
        }
    }
    else {
        LDEBUG("Output directory: " << outputDir);
        if (!FileSystem::clearDirectory(outputDir))
            LWARNING("Failed to clear output directory '" + outputDir + "'");
    }
    testCase.outputDir_ = outputDir;

    // run test workspace (write log and console output into output dir)
    std::set<std::string> ignoreFiles;
    testCase.htmlLog_ = baseName + "-log.html";
    std::string htmlLogAbs = FileSystem::cleanupPath(outputDir + "/" + testCase.htmlLog_);
    ignoreFiles.insert(FileSystem::fileName(htmlLogAbs));
    std::string consoleLogAbs;
    if (redirectStdOut) {
        testCase.consoleLog_ = baseName + "-console.txt";
        consoleLogAbs = FileSystem::cleanupPath(outputDir + "/" + testCase.consoleLog_);
        ignoreFiles.insert(FileSystem::fileName(consoleLogAbs));
    }
    bool useCaching = VoreenApplication::app()->useCaching();
    int timeout = testCase.config_.timeout_;
    tgtAssert(timeout >= 0, "invalid timeout value");
    executeWorkspace(voreentoolPath, testCase, testdataDir, useCaching, timeout, htmlLogAbs, consoleLogAbs, glclsharing);
    if (testCase.returnCode_ == 0) {
        LDEBUG("Test run successful (return code: 0)");

        // compare result with reference data
        LDEBUG("Verifying output data ...");
        bool cmpResult = compareResults(testCase, ignoreFiles, fileComparators);
        if (testCase.config_.ignored_) {
            LINFO("Test IGNORED (" << dtos(testCase.duration_) << " sec)\n");
            testCase.result_ = RegressionTestIgnored;
        }
        else if (cmpResult) {
            LINFO("Test SUCCEEDED (" << dtos(testCase.duration_) << " sec)\n");
            testCase.result_ = RegressionTestSuccess;
        }
        else {
            LWARNING("Test FAILED (" << dtos(testCase.duration_) << " sec)");
            LWARNING(ConsoleReporter().generateFileComparisonReport(testCase, false) << "\n");
            testCase.result_ = RegressionTestFailure;
        }
        LDEBUG(ConsoleReporter().generateFileComparisonReport(testCase, true));
    }
    else {
        testCase.errorMsg_ = "Failed to run test. voreentool return code: " + itos(testCase.returnCode_);
        LERROR("Test ERROR: " << testCase.errorMsg_ << " (" << dtos(testCase.duration_) << " sec)\n");
        testCase.result_ = RegressionTestError;
    }

    // copy testfile to output dir
    std::string testFileOutput = FileSystem::cleanupPath(outputDir + "/" + FileSystem::fileName(testCase.testfile_));
    try {
        LDEBUG("Copying test file '" << testCase.testfile_ << "' to output dir: " << outputDir);
        FileSystem::copyFile(testCase.testfile_, testFileOutput);
    }
    catch (tgt::Exception& e) {
        LWARNING("Failed to copy test file '" << testCase.testfile_ << "' to '" << testFileOutput << "': " << e.what());
    }
    // copy configfile (if present) to output dir
    if (testCase.configfile_ != "") {
        std::string configFileOutput = FileSystem::cleanupPath(outputDir + "/" + FileSystem::fileName(testCase.configfile_));
        try {
            LDEBUG("Copying test config file '" << testCase.configfile_ << "' to output dir: " << outputDir);
            FileSystem::copyFile(testCase.configfile_, configFileOutput);
        }
        catch (tgt::Exception& e) {
            LWARNING("Failed to copy test config file '" << testCase.configfile_ << "' to '" << configFileOutput << "': " << e.what());
        }
    }

    // move output and reference data to report dir
    if (!reportDir.empty()) {
        std::string testOutputDir = FileSystem::cleanupPath(reportDir + "/output/" + testSubDir);
        LDEBUG("Moving output data to report directory '" << testOutputDir << "'");

        // create/clear report output dir for test case
        if (FileSystem::dirExists(testOutputDir)) {
            LDEBUG("Clearing report output dir " << testOutputDir);
            FileSystem::clearDirectory(testOutputDir);
        }
        else {
            LDEBUG("Creating report output dir " << testOutputDir);
            if (!FileSystem::createDirectoryRecursive(testOutputDir)) {
                LWARNING("Failed to create report output directory: " << testOutputDir);
                testOutputDir = "";
            }
        }

        // move output data to report dir
        if (testOutputDir != "") {
            try {
                moveDirectoryContents(outputDir, testOutputDir);
                testCase.outputDir_ = testOutputDir;
            }
            catch (VoreenException& e) {
                LWARNING("Failed to move output data to report output directory: " << e.what());
            }
        }

        // copy reference data to report dir
        if (keepReferenceData) {
            // create/clear report reference dir for test case
            std::string testReferenceDir = FileSystem::cleanupPath(reportDir + "/reference/" + testSubDir);
            if (FileSystem::dirExists(testReferenceDir)) {
                LDEBUG("Clearing report reference dir " << testReferenceDir);
                FileSystem::clearDirectory(testReferenceDir);
            }
            else {
                LDEBUG("Creating report reference dir " << testReferenceDir);
                if (!FileSystem::createDirectoryRecursive(testReferenceDir)) {
                    LWARNING("Failed to create report reference directory: " << testReferenceDir);
                    testReferenceDir = "";
                }
            }

            // move reference data to report dir
            if (testReferenceDir != "") {
                try {
                    copyDirectoryContents(refDir, testReferenceDir);
                    testCase.referenceDir_ = testReferenceDir;
                }
                catch (VoreenException& e) {
                    LWARNING("Failed to move reference data to report output directory: " << e.what());
                }
            }
        }
    }
    else {
        testCase.outputDir_ = "";   //< without report directory, output data has to be discarded
    }

}

std::map<Processor*, std::vector<RegressionTestCase> > getCoverageMap(const RegressionTestSuite& testSuite, const std::vector<Processor*> processors) {

    tgtAssert(VoreenApplication::app(), "VoreenApplication not instantiated");
    std::string vrnBasePath = VoreenApplication::app()->getBasePath();

    std::map<Processor*, std::vector<RegressionTestCase> > coverageMap;
    std::map<std::string, Processor*> processorNameMap;

    // create map entries for passed processors
    for (size_t i=0; i<processors.size(); i++) {
        coverageMap.insert(std::make_pair(processors.at(i), std::vector<RegressionTestCase>()));
        processorNameMap.insert(std::make_pair(processors.at(i)->getClassName(), processors.at(i)));
    }

    // iterate over test cases
    std::vector<RegressionTestCase>::const_iterator testcaseIter;
    for (testcaseIter = testSuite.testCases_.begin(); testcaseIter != testSuite.testCases_.end(); ++testcaseIter) {
        const RegressionTestCase& testCase = *testcaseIter;
        if (testCase.testfile_ == "") {
            LWARNING("Test case '" << testCase.name_ << "' has not test file");
            continue;
        }

        // deserialize workspace
        Workspace workspace;
        try {
            LINFO("Loading test workspace: " << FileSystem::relativePath(testCase.testfile_, vrnBasePath));
            workspace.load(testCase.testfile_);
        }
        catch (SerializationException& e) {
            LERROR(e.what());
            continue;
        }

        // add contained processors to coverage map
        tgtAssert(workspace.getProcessorNetwork(), "workspace has no processor network");
        const std::vector<Processor*> testcaseProcs = workspace.getProcessorNetwork()->getProcessors();
        for (size_t i=0; i<testcaseProcs.size(); i++) {
            // ignore processors not contained by the passed list
            if (processorNameMap.find(testcaseProcs.at(i)->getClassName()) == processorNameMap.end())
                continue;
            Processor* proc = processorNameMap[testcaseProcs.at(i)->getClassName()];

            // ignore processors that are not part of the testcase's module
            if (proc->getModuleName() != testCase.moduleName_)
                continue;

            if (coverageMap.find(proc) == coverageMap.end())
                coverageMap.insert(std::make_pair(proc, std::vector<RegressionTestCase>()));
            std::vector<RegressionTestCase>& procTestCases = coverageMap[proc];
            if (std::find(procTestCases.begin(), procTestCases.end(), testCase) == procTestCases.end())
                procTestCases.push_back(testCase);
        }
    }

    return coverageMap;
}

void executeWorkspace(const std::string& binary, RegressionTestCase& testCase, const std::string& workDir, bool useCaching,
        int timeout, const std::string& logFile, const std::string& stdOutFile, bool glclsharing) {
    tgtAssert(!binary.empty(), "binary string empty");
    tgtAssert(!testCase.testfile_.empty(), "testfile string empty");

    testCase.call_ = "";

#ifdef WIN32
    if (timeout > 0)
        LDEBUG("Test timeout not supported on Windows");
#else
    if (timeout > 0)
        testCase.call_ += "timeout " + itos(timeout) + "s ";
#endif

    testCase.call_ += binary;

    testCase.call_ += " -w " + testCase.testfile_;
    if (!workDir.empty())
        testCase.call_ += " --workdir " + workDir;
    testCase.call_ += " --useCaching " + std::string((useCaching ? "true" : "false"));
    if (testCase.renderingTest_) {
        testCase.call_ += " --opengl";
        testCase.call_ += " --trigger-imagesaves";
        if(!glclsharing) {
            testCase.call_ += " --noglclsharing";
        }
    }
    testCase.call_ += " --trigger-volumesaves";
    testCase.call_ += " --trigger-geometrysaves";
    if (!logFile.empty())
        testCase.call_ += " --logFile " + logFile;
    if (!stdOutFile.empty()) {
#ifdef WIN32
        testCase.call_ += " 1>" + stdOutFile + " 2>&1";
#else
        testCase.call_ += " > " + stdOutFile + " 2>&1";
#endif
    }

    // Wait for asynccompute processors to actually process their output before ending network evaluation
#ifdef WIN32
    _putenv_s("VRN_ASYNCCOMPUTEPROCESSOR_WAIT_FOR_RESULT", "1");
#else
    setenv("VRN_ASYNCCOMPUTEPROCESSOR_WAIT_FOR_RESULT", "1", 1);
#endif

    LDEBUG("Calling: " << testCase.call_);

    tgt::Stopwatch stopWatch;
    stopWatch.start();
    testCase.returnCode_ = system(testCase.call_.c_str());
    //testCase.returnCode_ = 1;
    stopWatch.stop();
    testCase.duration_ = static_cast<double>(stopWatch.getRuntime()) / 1000;
}

bool compareResults(RegressionTestCase& testCase, const std::set<std::string>& ignoreFiles,
    const std::vector<FileComparator*>& fileComparators)
{
    tgtAssert(!fileComparators.empty(), "no FileComparators");

    // directories should have been checked by caller
    tgtAssert(FileSystem::dirExists(testCase.outputDir_), "Output directory does not exist");
    tgtAssert(FileSystem::dirExists(testCase.referenceDir_), "Reference directory does not exist");

    std::vector<std::string> errors;

    // retrieve files from directories (sorted)
    const std::vector<std::string> referenceFilesComplete = FileSystem::listFiles(testCase.referenceDir_, true);
    const std::vector<std::string> outputFilesComplete = FileSystem::listFiles(testCase.outputDir_, true);
    // filter ignored files
    std::vector<std::string> referenceFiles;
    for (size_t i=0; i<referenceFilesComplete.size(); i++)
        if (ignoreFiles.find(referenceFilesComplete.at(i)) == ignoreFiles.end())
            referenceFiles.push_back(referenceFilesComplete.at(i));
    size_t numRefFiles = referenceFiles.size();
    std::vector<std::string> outputFiles;
    for (size_t i=0; i<outputFilesComplete.size(); i++)
        if (ignoreFiles.find(outputFilesComplete.at(i)) == ignoreFiles.end())
            outputFiles.push_back(outputFilesComplete.at(i));
    size_t numOutputFiles = outputFiles.size();

    // group reference and output files by basename
    std::map<std::string, std::vector<std::string> > referenceGroups;
    std::map<std::string, std::vector<std::string> > outputGroups;
    for (size_t i=0; i<referenceFiles.size(); i++) {
        std::string baseName = FileSystem::baseName(referenceFiles.at(i));
        tgtAssert(baseName != "", "basename is empty");
        if (referenceGroups.find(baseName) == referenceGroups.end())
            referenceGroups[baseName] = std::vector<std::string>();
        referenceGroups[baseName].push_back(referenceFiles.at(i));
    }
    for (size_t i=0; i<outputFiles.size(); i++) {
        std::string baseName = FileSystem::baseName(outputFiles.at(i));
        tgtAssert(baseName != "", "basename is empty");
        if (outputGroups.find(baseName) == outputGroups.end())
            outputGroups[baseName] = std::vector<std::string>();
        outputGroups[baseName].push_back(outputFiles.at(i));
    }

    // determine file groups that are present in reference and output data
    std::vector<std::string> correspondingGroups;
    std::map<std::string, std::vector<std::string> >::iterator groupIter;
    for (groupIter = referenceGroups.begin(); groupIter != referenceGroups.end(); ++groupIter) {
        std::string baseName = groupIter->first;
        tgtAssert(groupIter->second.size() > 0, "reference file group is empty");
        if (outputGroups.find(baseName) != outputGroups.end())
            correspondingGroups.push_back(baseName);
    }
    std::sort(correspondingGroups.begin(), correspondingGroups.end());

    // run comparison on corresponding file groups
    for (std::vector<std::string>::iterator it = correspondingGroups.begin(); it != correspondingGroups.end(); ++it) {
        std::string groupBasename = *it;
        tgtAssert(referenceGroups.find(groupBasename) != referenceGroups.end(), "file group missing in reference groups");
        tgtAssert(outputGroups.find(groupBasename) != outputGroups.end(), "file group missing in output groups");

        std::vector<std::string> refGroupFiles = referenceGroups[groupBasename];
        std::vector<std::string> outputGroupFiles = outputGroups[groupBasename];
        tgtAssert(refGroupFiles.size() > 0, "ref group is empty");
        tgtAssert(outputGroupFiles.size() > 0, "ref group is empty");

        // create comparison datasets from file groups
        FileComparisonResult comparisonResult;
        for (size_t i=0; i<refGroupFiles.size(); i++) {
            RegressionTestFile file;
            file.filename_ = refGroupFiles.at(i);
            comparisonResult.refDataset_.files_.push_back(file);
        }
        for (size_t i=0; i<outputGroupFiles.size(); i++) {
            RegressionTestFile file;
            file.filename_ = outputGroupFiles.at(i);
            comparisonResult.outputDataset_.files_.push_back(file);
        }

        // check first, if data sets' filenames match
        comparisonResult.match_ = true;
        if (comparisonResult.refDataset_.files_.size() != comparisonResult.outputDataset_.files_.size()) {
            comparisonResult.match_ = false;
            comparisonResult.report_ = "output file count differs from reference file count";
            determineFileType(comparisonResult.refDataset_, fileComparators);
            determineFileType(comparisonResult.outputDataset_, fileComparators);
            testCase.mismatchingDatasets_.push_back(comparisonResult);
        }
        else {
            for (size_t i=0; i<comparisonResult.refDataset_.files_.size() && comparisonResult.match_; i++) {
                std::string refFile = comparisonResult.refDataset_.files_.at(i).filename_;
                std::string outputFile = comparisonResult.outputDataset_.files_.at(i).filename_;
                if (refFile != outputFile) {
                    comparisonResult.match_ = false;
                    comparisonResult.report_ = "output filenames do not match reference filenames";
                    determineFileType(comparisonResult.refDataset_, fileComparators);
                    determineFileType(comparisonResult.outputDataset_, fileComparators);
                    testCase.mismatchingDatasets_.push_back(comparisonResult);
                }
            }
        }
        if (!comparisonResult.match_)
            continue;

        // iterate over FileComparators and use first matching one
        bool compared = false;
        for (size_t i=0; i<fileComparators.size() && !compared; i++) {
            FileComparator* comparator = fileComparators.at(i);
            tgtAssert(comparator, "null pointer");
            if (!comparator->supportsFormat(comparisonResult.refDataset_))
                continue;

            // run comparison
            comparisonResult.match_ = comparator->compare(comparisonResult.refDataset_, comparisonResult.outputDataset_,
                comparisonResult.report_, testCase);
            compared = true;

            // if files do not match, generate diff file
            if (!comparisonResult.match_) {
                if (comparator->supportsFileDiff()) {
                    try {
                        LDEBUG("Generating diff file for ref dataset: " << strJoin(comparisonResult.refDataset_.files_, ", "));
                        comparator->generateDiffFile(comparisonResult.refDataset_, comparisonResult.outputDataset_,
                            comparisonResult.diffDataset_, testCase);
                    }
                    catch (VoreenException& e) {
                        LWARNING("Failed to generate diff file of '" << strJoin(comparisonResult.refDataset_.files_, ", ") << "' and '"
                            << strJoin(comparisonResult.outputDataset_.files_, ", ") << "': " << e.what());
                    }
                }
            }
        }
        if (!compared) {
            LERROR("No compatible FileComparator found for: " << strJoin(comparisonResult.refDataset_.files_, ", "));
            comparisonResult.match_ = false;
        }

        //  store comparison result
        if (comparisonResult.match_)
            testCase.matchingDatasets_.push_back(comparisonResult);
        else
            testCase.mismatchingDatasets_.push_back(comparisonResult);
    }

    // add reference file groups without counterpart in output files
    for (groupIter = referenceGroups.begin(); groupIter != referenceGroups.end(); ++groupIter) {
        std::string groupBaseName = groupIter->first;
        if (std::find(correspondingGroups.begin(), correspondingGroups.end(), groupBaseName) == correspondingGroups.end()) {
            std::vector<std::string> groupFiles = groupIter->second;
            RegressionTestDataset dataSet;
            for (size_t i=0; i<groupFiles.size(); i++) {
                RegressionTestFile file;
                file.filename_ = groupFiles.at(i);
                dataSet.files_.push_back(file);
            }
            determineFileType(dataSet, fileComparators);
            testCase.missingRefDatasets_.push_back(dataSet);
        }
    }

    // add output file groups without counterpart in reference files
    for (groupIter = outputGroups.begin(); groupIter != outputGroups.end(); ++groupIter) {
        std::string groupBaseName = groupIter->first;
        if (std::find(correspondingGroups.begin(), correspondingGroups.end(), groupBaseName) == correspondingGroups.end()) {
            std::vector<std::string> groupFiles = groupIter->second;
            RegressionTestDataset dataSet;
            for (size_t i=0; i<groupFiles.size(); i++) {
                RegressionTestFile file;
                file.filename_ = groupFiles.at(i);
                dataSet.files_.push_back(file);
            }
            determineFileType(dataSet, fileComparators);
            testCase.unexpectedOutputDatasets_.push_back(dataSet);
        }
    }

    // check whether a dataset has been created from all file groups
    tgtAssert(referenceGroups.size() == testCase.matchingDatasets_.size() + testCase.mismatchingDatasets_.size() + testCase.missingRefDatasets_.size(),
        "reference dataset count does not match number of classified datasets");
    tgtAssert(outputGroups.size() == testCase.matchingDatasets_.size() + testCase.mismatchingDatasets_.size() + testCase.unexpectedOutputDatasets_.size(),
        "output dataset count does not match number of classified datasets");

    // check whether all reference and output files have been assigned to a dataset
    size_t numMatchingRefFiles = 0;
    size_t numMatchingOutputFiles = 0;
    for (size_t i=0; i < testCase.matchingDatasets_.size(); i++) {
        numMatchingRefFiles += testCase.matchingDatasets_.at(i).refDataset_.files_.size();
        numMatchingOutputFiles += testCase.matchingDatasets_.at(i).outputDataset_.files_.size();
    }
    size_t numMismatchingRefFiles = 0;
    size_t numMismatchingOutputFiles = 0;
    for (size_t i=0; i < testCase.mismatchingDatasets_.size(); i++) {
        numMismatchingRefFiles += testCase.mismatchingDatasets_.at(i).refDataset_.files_.size();
        numMismatchingOutputFiles += testCase.mismatchingDatasets_.at(i).outputDataset_.files_.size();
    }
    size_t numMissingRefFiles = 0;
    for (size_t i=0; i < testCase.missingRefDatasets_.size(); i++)
        numMissingRefFiles += testCase.missingRefDatasets_.at(i).files_.size();
    size_t numUnexpectedOutputFiles = 0;
    for (size_t i=0; i < testCase.unexpectedOutputDatasets_.size(); i++)
        numUnexpectedOutputFiles += testCase.unexpectedOutputDatasets_.at(i).files_.size();
    tgtAssert(numRefFiles == numMatchingRefFiles + numMismatchingRefFiles + numMissingRefFiles,
        "reference file count does not match number of classified ref files");
    tgtAssert(numOutputFiles == numMatchingOutputFiles + numMismatchingOutputFiles + numUnexpectedOutputFiles,
        "output file count does not match number of classified output files");

    return testCase.mismatchingDatasets_.empty() && testCase.missingRefDatasets_.empty() && testCase.unexpectedOutputDatasets_.empty();
}

void moveDirectoryContents(const std::string& srcDir, const std::string& destDir) {
    if (!FileSystem::dirExists(srcDir))
        throw VoreenException("Src directory does not exist: " + srcDir);
    if (!FileSystem::dirExists(destDir))
        throw VoreenException("Dest directory does not exist: " + destDir);

    std::vector<std::string> files = FileSystem::listFiles(srcDir);
    for (size_t i=0; i<files.size(); i++) {
        std::string srcFile = FileSystem::cleanupPath(srcDir + "/" + files.at(i));
        std::string destFile = FileSystem::cleanupPath(destDir + "/" + files.at(i));

        LDEBUG("Moving file '" << srcFile << "' to '" << destFile << "'");

        if (FileSystem::fileExists(destFile)) {
            LDEBUG("Deleting dest file: " << destFile);
            FileSystem::deleteFile(destFile);
        }

        if (!FileSystem::renameFile(srcFile, destFile, false))
            LWARNING("Failed to move file '" << srcFile << "' to '" << destFile << "'");
    }
}

void copyDirectoryContents(const std::string& srcDir, const std::string& destDir) {
    if (!FileSystem::dirExists(srcDir))
        throw VoreenException("Src directory does not exist: " + srcDir);
    if (!FileSystem::dirExists(destDir))
        throw VoreenException("Dest directory does not exist: " + destDir);

    std::vector<std::string> files = FileSystem::listFiles(srcDir);
    for (size_t i=0; i<files.size(); i++) {
        std::string srcFile = FileSystem::cleanupPath(srcDir + "/" + files.at(i));
        std::string destFile = FileSystem::cleanupPath(destDir + "/" + files.at(i));

        std::string subDir = FileSystem::dirName(files.at(i));
        if (subDir != "") {
            LDEBUG("Creating subdir '" + subDir + "' in destination directory: " + destDir);
            std::string destSubDir = FileSystem::cleanupPath(destDir + "/" + subDir);
            if (!FileSystem::createDirectoryRecursive(destSubDir)) {
                LWARNING("Failed to create directory: " << destSubDir);
                continue;
            }
        }

        LDEBUG("Copying file '" << srcFile << "' to '" << destFile << "'");

        if (FileSystem::fileExists(destFile)) {
            LDEBUG("Deleting dest file: " << destFile);
            FileSystem::deleteFile(destFile);
        }

        try {
            FileSystem::copyFile(srcFile, destFile);
        }
        catch (tgt::Exception& e) {
            LWARNING("Failed to copy file '" << srcFile << "' to '" << destFile << "': " << e.what());
        }
    }
}

void determineFileType(RegressionTestDataset& dataset, const std::vector<FileComparator*>& fileComparators) {
    for (size_t i=0; i<fileComparators.size(); i++) {
        FileComparator* comparator = fileComparators.at(i);
        tgtAssert(comparator, "comparator is null pointer");
        if (comparator->supportsFormat(dataset)) {
            comparator->determineFileType(dataset);
            return;
        }
    }

    LWARNING("No matching FileComparator found for dataset: " << strJoin(dataset.files_, ", "));
}
