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

#include "reportgenerators.h"

#include "tgt/logmanager.h"
#include "tgt/filesystem.h"
#include "tgt/vector.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/serialization/xmlserializationconstants.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/voreenmodule.h"

#include <map>
#include <iostream>

namespace voreen {

using namespace tgt;

std::string ConsoleReporter::generateFileComparisonReport(const RegressionTestCase& testCase, bool includeMatchingFiles) {
    const std::string LIST_SYMB =   "- ";
    const std::string LIST_INDENT = "  ";

    std::ostringstream reportStr;

    size_t numRefFiles = testCase.matchingDatasets_.size() +
        testCase.mismatchingDatasets_.size() +
        testCase.missingRefDatasets_.size();

    if (includeMatchingFiles) {
        reportStr << "\nMatching Datasets: " << testCase.matchingDatasets_.size() << "/" << numRefFiles << "\n";
        for (size_t i=0; i<testCase.matchingDatasets_.size(); i++) {
            FileComparisonResult comparison = testCase.matchingDatasets_.at(i);
            std::vector<RegressionTestFile> refFileGroup = comparison.refDataset_.files_;
            tgtAssert(refFileGroup.size() > 0, "refDataset does not contain any files");
            reportStr << LIST_SYMB << refFileGroup.at(0).filename_ << "\t: " << comparison.report_ << "\n";
            for (size_t j=1; j<refFileGroup.size(); j++)
                reportStr << LIST_INDENT << refFileGroup.at(j).filename_ << "\n";
        }
    }

    if (!testCase.mismatchingDatasets_.empty()) {
        reportStr << "\nMismatching Datasets: " << testCase.mismatchingDatasets_.size() << "/" << numRefFiles << "\n";
        for (size_t i=0; i<testCase.mismatchingDatasets_.size(); i++) {
            FileComparisonResult comparison = testCase.mismatchingDatasets_.at(i);
            std::vector<RegressionTestFile> refFileGroup = comparison.refDataset_.files_;
            tgtAssert(refFileGroup.size() > 0, "refData does not contain any files");
            reportStr << LIST_SYMB << refFileGroup.at(0).filename_ << "\t: " << comparison.report_ << "\n";
            for (size_t j=1; j<refFileGroup.size(); j++)
                reportStr << LIST_INDENT << refFileGroup.at(j).filename_ << "\n";
        }
    }

    if (!testCase.missingRefDatasets_.empty()) {
        reportStr << "\nReference Datasets Missing in Output: " << testCase.missingRefDatasets_.size() << "/" << numRefFiles << "\n";
        for (size_t i=0; i<testCase.missingRefDatasets_.size(); i++) {
            std::vector<RegressionTestFile> fileGroup = testCase.missingRefDatasets_.at(i).files_;
            tgtAssert(fileGroup.size() > 0, "dataset does not contain any files");
            reportStr << LIST_SYMB << fileGroup.at(0).filename_ << "\n";
            for (size_t j=1; j<fileGroup.size(); j++)
                reportStr << LIST_INDENT << fileGroup.at(j).filename_ << "\n";
        }
    }

    if (!testCase.unexpectedOutputDatasets_.empty()) {
        reportStr << "\nUnexpected Output Datasets: " << testCase.unexpectedOutputDatasets_.size() << "\n";

        for (size_t i=0; i<testCase.unexpectedOutputDatasets_.size(); i++) {
            std::vector<RegressionTestFile> fileGroup = testCase.unexpectedOutputDatasets_.at(i).files_;
            tgtAssert(fileGroup.size() > 0, "dataset does not contain any files");
            reportStr << LIST_SYMB << fileGroup.at(0).filename_ << "\n";
            for (size_t j=1; j<fileGroup.size(); j++)
                reportStr << LIST_INDENT << fileGroup.at(j).filename_ << "\n";
        }
    }

    return reportStr.str();
}

std::string ConsoleReporter::generateTestCoverageReport(const std::vector<VoreenModule*>& modules,
    const std::map<Processor*, std::vector<RegressionTestCase> >& coverageMap,
    const std::map<VoreenModule*, std::vector<Processor*> >& ignoreLists)
{
    const std::string LIST_SYMB =   "- ";
    const std::string LIST_INDENT = "  ";

    std::string report;

    size_t numTested = 0;
    size_t numTotal = 0;
    size_t numStableTested = 0;     //< procs with status STABLE/TESTING that are tested
    size_t numStableUntested = 0;   //< procs with status STABLE/TESTING that are not tested
    size_t numUnstableTested = 0;   //< procs with status below TESTING that are tested
    size_t numUnstableUntested = 0; //< procs with status below TESTING that are not tested
    size_t numIgnored = 0;          //< procs listed in the modules coverage-ignore.xml
    // module reports
    for (size_t i=0; i<modules.size(); i++) {
        VoreenModule* module = modules.at(i);
        numTotal += module->getRegisteredProcessors().size();
        std::vector<Processor*> ignoredProcs;
        if (ignoreLists.find(module) != ignoreLists.end())
            ignoredProcs = ignoreLists.find(module)->second;

        // classify each of the module's processors into on of these four classes:
        std::vector<Processor*> stableTested;     //< procs with status STABLE/TESTING that are tested
        std::vector<Processor*> stableUntested;   //< procs with status STABLE/TESTING that are not tested
        std::vector<Processor*> unstableTested;   //< procs with status below TESTING that are tested
        std::vector<Processor*> unstableUntested; //< procs with status below TESTING that are not tested
        std::vector<Processor*> ignored;          //< procs listed in the modules coverage-ignore.xml
        for (size_t j=0; j<module->getRegisteredProcessors().size(); j++) {
            Processor* processor = const_cast<Processor*>(module->getRegisteredProcessors().at(j));
            tgtAssert(processor, "processor is null");
            if (coverageMap.find(processor) == coverageMap.end()) {
                LERRORC("regressiontest.ConsoleReporter", "Processor not found in coverage map: " << processor->getClassName());
                continue;
            }

            // check, if processor is in ignore list
            if (std::find(ignoredProcs.begin(), ignoredProcs.end(), processor) != ignoredProcs.end()) {
                ignored.push_back(processor);
                numIgnored++;
                continue;
            }

            // detemine testing state
            bool procTested = false;
            const std::vector<RegressionTestCase> testCases = coverageMap.find(processor)->second;
            for (size_t k=0; k<testCases.size() && !procTested; k++) {
                if (testCases.at(k).config_.enabled_)
                    procTested = true;
            }

            if (processor->getCodeState() >= Processor::CODE_STATE_TESTING) {
                if (procTested) {
                    stableTested.push_back(processor);
                    numStableTested++;
                }
                else {
                    stableUntested.push_back(processor);
                    numStableUntested++;
                }
            }
            else {
                if (procTested) {
                    unstableTested.push_back(processor);
                    numUnstableTested++;
                }
                else {
                    unstableUntested.push_back(processor);
                    numUnstableUntested++;
                }
            }
        }

        // generate report for this module
        size_t numTestedModule = stableTested.size() + unstableTested.size();
        numTested += numTestedModule;
        report += "\nModule '" + module->getDirName() + "': " + itos(numTestedModule) + "/" + itos(module->getRegisteredProcessors().size())
            + " Processors Tested\n";
        report += std::string(50, '-') + "\n";

        report += "Stable/Testing processors that are tested: " + itos(stableTested.size()) + " \n";
        for (size_t j=0; j<stableTested.size(); j++)
            report += LIST_SYMB + stableTested.at(j)->getClassName() + "\n";

        report += "Stable/Testing processors that are NOT tested: " + itos(stableUntested.size()) + " \n";
        for (size_t j=0; j<stableUntested.size(); j++)
            report += LIST_SYMB + stableUntested.at(j)->getClassName() + "\n";

        report += "Experimental processors that are tested: " + itos(unstableTested.size()) + " \n";
        for (size_t j=0; j<unstableTested.size(); j++)
            report += LIST_SYMB + unstableTested.at(j)->getClassName() + "\n";

        report += "Experimental processors that are not tested: " + itos(unstableUntested.size()) + " \n";
        for (size_t j=0; j<unstableUntested.size(); j++)
            report += LIST_SYMB + unstableUntested.at(j)->getClassName() + "\n";

        report += "Ignored processors: " + itos(ignored.size()) + " \n";
        for (size_t j=0; j<ignored.size(); j++)
            report += LIST_SYMB + ignored.at(j)->getClassName() + "\n";
    }

    // add summary
    std::string summary = "SUMMARY: Num Processors Tested/Total:\t\t " + itos(numTested) + "/" + itos(numTotal) + "\n";
    summary += LIST_SYMB + "Stable/Testing processors that are tested:\t " + itos(numStableTested) + " \n";
    summary += LIST_SYMB + "Stable/Testing processors that are NOT tested: " + itos(numStableUntested) + " \n";
    summary += LIST_SYMB + "Experimental processors that are tested:\t " + itos(numUnstableTested) + " \n";
    summary += LIST_SYMB + "Experimental processors that are not tested:\t " + itos(numUnstableUntested) + " \n";
    summary += LIST_SYMB + "Ignored processors (coverage-ignore.xml):\t " + itos(numIgnored);

    report += "\n" + std::string(56, '-') + "\n" + summary + "\n" + std::string(56, '-');

    return report;
}

//----------------------------------------------------------------------------------------

void HTMLReportGenerator::generateTestResultReport(const std::string& htmlFilename, const RegressionTestSuite& testSuite, bool splitReport) {
    const int IMG_WIDTH = 0; //< set in stylesheet

    // columns widths of overview table
    std::map<std::string, std::string> overviewColWidths;
    overviewColWidths["status"] = "80px";
    overviewColWidths["testname"] = "200px";
    overviewColWidths["time"] = "50px";

    LINFOC("regressiontest.HTMLReportGenerator", "Generating HTML report: " << htmlFilename);
    const std::string documentPath = FileSystem::dirName(htmlFilename);
    const std::string documentBaseName = FileSystem::baseName(htmlFilename);

    // copy stylesheet file to document path
    try {
        std::string cssSrc = VoreenApplication::app()->getBasePath("apps/tests/regressiontest/report.css");
        std::string cssDest = documentPath + "/report.css";
        FileSystem::copyFile(cssSrc, cssDest);
    }
    catch (tgt::Exception& e) {
        throw VoreenException("Failed to copy style sheet: " + std::string(e.what()));
    }

    // group test cases by module
    const std::vector<VoreenModule*>& modules = VoreenApplication::app()->getModules();
    std::map<const VoreenModule*, std::vector<RegressionTestCase> > testCasesMap;
    for (size_t i=0; i<modules.size(); i++) {
        VoreenModule* module = modules.at(i);
        std::vector<RegressionTestCase> moduleTestCases;
        for (size_t t=0; t<testSuite.testCases_.size(); t++) {
            if (testSuite.testCases_.at(t).moduleDir_ == module->getDirName())
                moduleTestCases.push_back(testSuite.testCases_.at(t));
        }
        if (!moduleTestCases.empty())
            testCasesMap[module] = moduleTestCases;
    }

    // create overview document
    TiXmlDocument* document = 0;
    TiXmlElement* htmlNode = 0;
    std::string titleStr = "Voreen Regression Test Report: " + testSuite.date_.toString(false);
    createDocument(htmlFilename, titleStr, document, htmlNode);
    tgtAssert(document && htmlNode, "null pointer returned");

    // body node
    TiXmlElement* body = new TiXmlElement("body");
    htmlNode->LinkEndChild(body);

    // title text
    TiXmlElement* titleDiv = new TiXmlElement("div");
    body->LinkEndChild(titleDiv);
    titleDiv->LinkEndChild(new TiXmlText(titleStr));
    titleDiv->SetAttribute("class", "title");

    // statistics
    TiXmlElement* statisticsElem = new TiXmlElement("div");
    body->LinkEndChild(statisticsElem);
    statisticsElem->SetAttribute("class", "statistics");

    TiXmlElement* spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Test Cases: " + itos(testSuite.testCases_.size())));
    statisticsElem->LinkEndChild(spanElem);
    statisticsElem->LinkEndChild(new TiXmlText(" | "));

    spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Success: " + itos(testSuite.numSuccess_)));
    if (testSuite.numSuccess_ > 0)
        spanElem->SetAttribute("class", "success");
    statisticsElem->LinkEndChild(spanElem);
    statisticsElem->LinkEndChild(new TiXmlText(" | "));

    spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Failed: " + itos(testSuite.numFailed_)));
    spanElem->SetAttribute("class", (testSuite.numFailed_ > 0 ? "failure" : "ignored"));
    statisticsElem->LinkEndChild(spanElem);
    statisticsElem->LinkEndChild(new TiXmlText(" | "));

    spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Error: " + itos(testSuite.numError_)));
    spanElem->SetAttribute("class", (testSuite.numError_ > 0 ? "error" : "ignored"));
    statisticsElem->LinkEndChild(spanElem);
    statisticsElem->LinkEndChild(new TiXmlText(" | "));

    spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Ignored: " + itos(testSuite.numIgnored_)));
    if (testSuite.numIgnored_ == 0)
        spanElem->SetAttribute("class", "ignored");
    statisticsElem->LinkEndChild(spanElem);
    statisticsElem->LinkEndChild(new TiXmlText(" | "));

    spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Skipped: " + itos(testSuite.numSkipped_)));
    if (testSuite.numSkipped_ == 0)
        spanElem->SetAttribute("class", "ignored");
    statisticsElem->LinkEndChild(spanElem);
    statisticsElem->LinkEndChild(new TiXmlText(" | "));

    spanElem = new TiXmlElement("span");
    spanElem->LinkEndChild(new TiXmlText("Duration: " + formatTime((size_t)tgt::iround(testSuite.duration_ * 1000))));
    statisticsElem->LinkEndChild(spanElem);

    // date
    /*TiXmlElement* dateElem = new TiXmlElement("div");
    body->LinkEndChild(dateElem);
    dateElem->SetAttribute("class", "date");
    dateElem->LinkEndChild(new TiXmlText(testSuite.date_.toString(false))); */

    //
    // test overview
    //
    TiXmlElement* curModuleContainer = 0;
    TiXmlElement* curResultTable = 0;
    for (size_t m=0; m<modules.size(); m++) {
        VoreenModule* module = modules.at(m);
        if (testCasesMap.find(module) == testCasesMap.end())
            continue;

        // start new module
        curModuleContainer = new TiXmlElement("div");
        body->LinkEndChild(curModuleContainer);
        curModuleContainer->SetAttribute("class", "overviewModuleContainer");
        /*if (i > 0)
            body->LinkEndChild(new TiXmlElement("hr")); */
        TiXmlElement* moduleHead = new TiXmlElement("div");
        curModuleContainer->LinkEndChild(moduleHead);
        moduleHead->SetAttribute("class", "overviewModuleHead");
        moduleHead->LinkEndChild(new TiXmlText("Module: " + module->getDirName()));

        std::string curGroup = "###";
        for (size_t i=0; i<testCasesMap[module].size(); i++) {
            const RegressionTestCase& testCase = testCasesMap[module].at(i);
            if (curGroup != testCase.group_) {
                // start new group
                TiXmlElement* groupContainer = new TiXmlElement("div");
                curModuleContainer->LinkEndChild(groupContainer);
                groupContainer->SetAttribute("class", "overviewGroupContainer");
                if (testCase.group_ != "") {
                    TiXmlElement* groupHead = new TiXmlElement("div");
                    groupContainer->LinkEndChild(groupHead);
                    groupHead->SetAttribute("class", "overviewGroupHead");
                    groupHead->LinkEndChild(new TiXmlText("Group: " + testCase.group_));
                }
                curGroup = testCase.group_;

                // start group result table
                curResultTable = new TiXmlElement("table");
                groupContainer->LinkEndChild(curResultTable);
                curResultTable->SetAttribute("class", "overviewTable");
                TiXmlElement* tableHead = new TiXmlElement("tr");
                curResultTable->LinkEndChild(tableHead);
                TiXmlElement* th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("status"));
                th->SetAttribute("width", overviewColWidths["status"]);
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("name"));
                th->SetAttribute("width", overviewColWidths["testname"]);
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("secs"));
                th->SetAttribute("width", overviewColWidths["time"]);
            } // group table initialization

            // create row for test case
            TiXmlElement* tr = new TiXmlElement("tr");
            curResultTable->LinkEndChild(tr);

            tgtAssert(curResultTable, "no result table");
            std::string resultStr;
            std::string styleClass;
            if (testCase.result_ == RegressionTestSuccess) {
                resultStr = "SUCCESS";
                styleClass = "success-bg";
            }
            else if (testCase.result_ == RegressionTestFailure) {
                resultStr = "FAILURE";
                styleClass = "failure-bg";
            }
            else if (testCase.result_ == RegressionTestIgnored) {
                resultStr = "IGNORED";
                styleClass = "ignored-bg";
            }
            else if (testCase.result_ == RegressionTestError) {
                resultStr = "ERROR";
                styleClass = "error-bg";
            }
            else if (testCase.result_ == RegressionTestSkipped) {
                resultStr = "SKIPPED";
                styleClass = "skipped-bg";
            }
            else
                tgtAssert(false, "unknown test result");

            tr->SetAttribute("class", styleClass);

            // result column
            TiXmlElement* td = new TiXmlElement("td");
            tr->LinkEndChild(td);
            td->LinkEndChild(new TiXmlText(resultStr));
            td->SetAttribute("align", "center");

            // name column with link
            td = new TiXmlElement("td");
            tr->LinkEndChild(td);
            TiXmlElement* a = new TiXmlElement("a");
            td->LinkEndChild(a);
            a->LinkEndChild(new TiXmlText(testCase.name_));
            // anchor
            const std::string testIDStr = testCase.moduleDir_ + "." + (!testCase.group_.empty() ? testCase.group_ + "." : "") + testCase.name_;
            std::string anchorStr;
            if (splitReport)
                anchorStr =  tgt::FileSystem::cleanupPath(documentBaseName + "-" + testIDStr + ".html");
            else
                anchorStr = "#" + testIDStr;
            a->SetAttribute("href", anchorStr);

            // time column
            td = new TiXmlElement("td");
            tr->LinkEndChild(td);
            td->LinkEndChild(new TiXmlText(dtos(testCase.duration_)));
        } // test case
    } // module

    body->LinkEndChild(new TiXmlElement("hr"));

    // write out root document, if in split mode
    if (splitReport) {
        saveDocument(document, htmlFilename);
        delete document;
        document = 0;
        htmlNode = 0;
        body = 0;
    }


    //
    // test details
    //
    for (size_t m=0; m<modules.size(); m++) {
        VoreenModule* module = modules.at(m);
        if (testCasesMap.find(module) == testCasesMap.end())
            continue;

        // iterate over module test cases
        std::vector<RegressionTestCase> testCases = testCasesMap[module];
        for (size_t t=0; t<testCases.size(); t++) {
            RegressionTestCase testCase = testCases.at(t);
            std::string testCaseRefDirRel;     //< path to test case reference dir, relative to document path
            std::string testCaseOutputDirRel;  //< path to test case output dir, relative to document path
            if (testCase.referenceDir_ != "")
                testCaseRefDirRel = FileSystem::relativePath(testCase.referenceDir_, documentPath);
            if (testCase.outputDir_ != "")
                testCaseOutputDirRel = FileSystem::relativePath(testCase.outputDir_, documentPath);

            // init document for testcase in split mode
            const std::string testCaseIDStr = testCase.moduleDir_ + "." + (!testCase.group_.empty() ? testCase.group_ + "." : "") + testCase.name_;
            const std::string testCaseHTMLFilename = tgt::FileSystem::cleanupPath(documentPath + "/" + documentBaseName + "-" + testCaseIDStr + ".html");
            if (splitReport) {
                std::string titleStr = "Regression Test:" + testCaseIDStr;
                createDocument(testCaseHTMLFilename, titleStr, document, htmlNode);
                tgtAssert(document && htmlNode, "null pointer returned");

                // body node
                body = new TiXmlElement("body");
                htmlNode->LinkEndChild(body);

                // title text
                /*TiXmlElement* titleDiv = new TiXmlElement("div");
                body->LinkEndChild(titleDiv);
                titleDiv->LinkEndChild(new TiXmlText(titleStr));
                titleDiv->SetAttribute("class", "title"); */
            }

            tgtAssert(body, "no body node");

            std::string resultStr;
            std::string styleClass;
            if (testCase.result_ == RegressionTestSuccess) {
                resultStr = "SUCCESS";
                styleClass = "success";
            }
            else if (testCase.result_ == RegressionTestFailure) {
                resultStr = "FAILURE";
                styleClass = "failure";
            }
            else if (testCase.result_ == RegressionTestIgnored) {
                resultStr = "IGNORED";
                styleClass = "ignored";
            }
            else if (testCase.result_ == RegressionTestError) {
                resultStr = "ERROR";
                styleClass = "error";
            }
            else if (testCase.result_ == RegressionTestSkipped) {
                resultStr = "SKIPPED";
                styleClass = "skipped";
            }
            else
                tgtAssert(false, "unknown test result");

            // outer div
            TiXmlElement* detailsContainer = new TiXmlElement("div");
            body->LinkEndChild(detailsContainer);
            detailsContainer->SetAttribute("class", "detailsContainer");

            // anchor
            TiXmlElement* anchor = new TiXmlElement("a");
            detailsContainer->LinkEndChild(anchor);
            std::string anchorStr = testCase.moduleDir_ + ".";
            if (testCase.group_ != "")
                anchorStr += testCase.group_ + ".";
            anchorStr += testCase.name_;
            anchor->SetAttribute("name", anchorStr);

            // title
            TiXmlElement* testTitle = new TiXmlElement("div");
            detailsContainer->LinkEndChild(testTitle);
            testTitle->SetAttribute("class", "detailsTitle " + styleClass + "-bg");
            std::string titleStr = resultStr + ": " + anchorStr;
            if (testCase.result_ != RegressionTestSkipped)
                titleStr += " (" + dtos(testCase.duration_) + " sec)";
            testTitle->LinkEndChild(new TiXmlText(titleStr));

            // details body
            /*TiXmlElement* detailsBody = new TiXmlElement("div");
            detailsContainer->LinkEndChild(detailsBody); */
            TiXmlElement* detailsBody = detailsContainer;

            // if test errored, print error msg
            if (testCase.result_ == RegressionTestError) {
                TiXmlElement* errorElem = new TiXmlElement("div");
                detailsBody->LinkEndChild(errorElem);
                errorElem->SetAttribute("class", "detailsErrorMessage");
                errorElem->LinkEndChild(new TiXmlText(testCase.errorMsg_));
            }

            // matching datasets
            if (!testCase.matchingDatasets_.empty()) {
                tgtAssert(testCaseRefDirRel != "", "TestCase reference dir not set");

                // listing container
                TiXmlElement* matchingFilesContainer = new TiXmlElement("div");
                detailsBody->LinkEndChild(matchingFilesContainer);
                matchingFilesContainer->SetAttribute("class", "detailsListingContainer");

                // title
                TiXmlElement* matchingFilesTitle = new TiXmlElement("div");
                matchingFilesContainer->LinkEndChild(matchingFilesTitle);
                matchingFilesTitle->LinkEndChild(new TiXmlText("Matching Datasets:"));
                matchingFilesTitle->SetAttribute("class", "detailsListingTitle");

                // table head
                TiXmlElement* resultTable = new TiXmlElement("table");
                matchingFilesContainer->LinkEndChild(resultTable);
                resultTable->SetAttribute("class", "detailsListingTable");
                TiXmlElement* tableHead = new TiXmlElement("tr");
                resultTable->LinkEndChild(tableHead);
                TiXmlElement* th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Output"));
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Reference"));
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Comparison Report"));
                //th->SetAttribute("width", "300px");

                // table body
                for (std::vector<FileComparisonResult>::const_iterator compIter = testCase.matchingDatasets_.begin();
                     compIter != testCase.matchingDatasets_.end(); ++compIter)
                {

                    const FileComparisonResult& comparison = *compIter;
                    TiXmlElement* tr = new TiXmlElement("tr");
                    resultTable->LinkEndChild(tr);

                    // output dataset (may not be available)
                    TiXmlElement* td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    for (size_t i=0; i<comparison.outputDataset_.files_.size(); i++) {
                        const RegressionTestFile& outputFile = comparison.outputDataset_.files_.at(i);
                        if (testCaseOutputDirRel != "") {
                            std::string outputFileHref = FileSystem::cleanupPath(testCaseOutputDirRel + "/" + outputFile.filename_);
                            TiXmlNode* innerNode;
                            if (outputFile.fileType_ == ImageFile)
                                innerNode = createImage(outputFileHref, outputFile.filename_, IMG_WIDTH);
                            else
                                innerNode = new TiXmlText(outputFile.filename_);
                            td->LinkEndChild(createLink(outputFileHref, innerNode));
                        }
                        else {
                            td->LinkEndChild(new TiXmlText(outputFile.filename_));
                        }
                        if (i < comparison.outputDataset_.files_.size()-1)
                            td->LinkEndChild(createLineBreak());
                    }

                    // reference dataset (always available)
                    td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    for (size_t i=0; i<comparison.refDataset_.files_.size(); i++) {
                        const RegressionTestFile& refFile = comparison.refDataset_.files_.at(i);
                        std::string refFileHref = FileSystem::cleanupPath(testCaseRefDirRel + "/" + refFile.filename_);
                        TiXmlNode* innerNode;
                        if (refFile.fileType_ == ImageFile)
                            innerNode = createImage(refFileHref, refFile.filename_, IMG_WIDTH);
                        else
                            innerNode = new TiXmlText(refFile.filename_);
                        td->LinkEndChild(createLink(refFileHref, innerNode));

                        if (i < comparison.refDataset_.files_.size()-1)
                            td->LinkEndChild(createLineBreak());
                    }

                    // comparison report
                    td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    td->LinkEndChild(new TiXmlText(comparison.report_));
                }
            }

            // mismatching datasets
            if (!testCase.mismatchingDatasets_.empty()) {
                tgtAssert(testCaseRefDirRel != "", "TestCase reference dir not set");

                // listing container
                TiXmlElement* mismatchingFilesContainer = new TiXmlElement("div");
                detailsBody->LinkEndChild(mismatchingFilesContainer);
                mismatchingFilesContainer->SetAttribute("class", "detailsListingContainer");

                // title
                TiXmlElement* mismatchingFilesTitle = new TiXmlElement("div");
                mismatchingFilesContainer->LinkEndChild(mismatchingFilesTitle);
                mismatchingFilesTitle->LinkEndChild(new TiXmlText("Mismatching Datasets:"));
                mismatchingFilesTitle->SetAttribute("class", "detailsListingTitle failure");

                // table head
                TiXmlElement* resultTable = new TiXmlElement("table");
                mismatchingFilesContainer->LinkEndChild(resultTable);
                resultTable->SetAttribute("class", "detailsListingTable");
                TiXmlElement* tableHead = new TiXmlElement("tr");
                resultTable->LinkEndChild(tableHead);
                TiXmlElement* th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Output"));
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Reference"));
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Difference"));
                th = new TiXmlElement("th");
                tableHead->LinkEndChild(th);
                th->LinkEndChild(new TiXmlText("Comparison Report"));

                // table body
                for (std::vector<FileComparisonResult>::const_iterator compIter = testCase.mismatchingDatasets_.begin();
                     compIter != testCase.mismatchingDatasets_.end(); ++compIter)
                {

                    const FileComparisonResult& comparison = *compIter;
                    TiXmlElement* tr = new TiXmlElement("tr");
                    resultTable->LinkEndChild(tr);

                    // output dataset (may not be available)
                    TiXmlElement* td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    for (size_t i=0; i<comparison.outputDataset_.files_.size(); i++) {
                        const RegressionTestFile& outputFile = comparison.outputDataset_.files_.at(i);
                        if (testCaseOutputDirRel != "") {
                            std::string outputFileHref = FileSystem::cleanupPath(testCaseOutputDirRel + "/" + outputFile.filename_);
                            TiXmlNode* innerNode;
                            if (outputFile.fileType_ == ImageFile)
                                innerNode = createImage(outputFileHref, outputFile.filename_, IMG_WIDTH);
                            else
                                innerNode = new TiXmlText(outputFile.filename_);
                            td->LinkEndChild(createLink(outputFileHref, innerNode));
                        }
                        else {
                            td->LinkEndChild(new TiXmlText(outputFile.filename_));
                        }

                        if (i < comparison.outputDataset_.files_.size()-1)
                            td->LinkEndChild(createLineBreak());
                    }

                    // reference dataset (always available)
                    td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    for (size_t i=0; i<comparison.refDataset_.files_.size(); i++) {
                        const RegressionTestFile& refFile = comparison.refDataset_.files_.at(i);
                        std::string refFileHref = FileSystem::cleanupPath(testCaseRefDirRel + "/" + refFile.filename_);
                        TiXmlNode* innerNode;
                        if (refFile.fileType_ == ImageFile)
                            innerNode = createImage(refFileHref, refFile.filename_, IMG_WIDTH);
                        else
                            innerNode = new TiXmlText(refFile.filename_);
                        td->LinkEndChild(createLink(refFileHref, innerNode));

                        if (i < comparison.refDataset_.files_.size()-1)
                            td->LinkEndChild(createLineBreak());
                    }

                    // difference dataset (may not be available)
                    td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    if (!comparison.diffDataset_.files_.empty() && testCaseOutputDirRel != "") {
                        for (size_t i=0; i<comparison.diffDataset_.files_.size(); i++) {
                            const RegressionTestFile& diffFile = comparison.diffDataset_.files_.at(i);
                            std::string diffFileHref = FileSystem::cleanupPath(testCaseOutputDirRel + "/" + diffFile.filename_);
                            TiXmlNode* innerNode;
                            if (diffFile.fileType_ == ImageFile)
                                innerNode = createImage(diffFileHref, diffFile.filename_, IMG_WIDTH);
                            else
                                innerNode = new TiXmlText(diffFile.filename_);
                            td->LinkEndChild(createLink(diffFileHref, innerNode));

                            if (i < comparison.diffDataset_.files_.size()-1)
                                td->LinkEndChild(createLineBreak());
                        }
                    }
                    else {
                        TiXmlText* text = new TiXmlText("<nobr>&nbsp; not available &nbsp;</nobr>");
                        text->SetCDATA(true);
                        td->LinkEndChild(text);
                        td->SetAttribute("style", "font-style: italic");
                    }

                    // comparison report
                    td = new TiXmlElement("td");
                    tr->LinkEndChild(td);
                    td->LinkEndChild(new TiXmlText(comparison.report_));
                    td->SetAttribute("class", "failure");
                }
            }

            // missing ref files
            if (!testCase.missingRefDatasets_.empty()) {
                tgtAssert(testCaseRefDirRel != "", "TestCase reference dir not set");

                // listing container
                TiXmlElement* missingRefFilesContainer = new TiXmlElement("div");
                detailsBody->LinkEndChild(missingRefFilesContainer);
                missingRefFilesContainer->SetAttribute("class", "detailsListingContainer");

                // title
                TiXmlElement* missingRefFilesTitle = new TiXmlElement("div");
                missingRefFilesContainer->LinkEndChild(missingRefFilesTitle);
                missingRefFilesTitle->LinkEndChild(new TiXmlText("Reference Datasets Missing in Output:"));
                missingRefFilesTitle->SetAttribute("class", "detailsListingTitle failure");

                // table head
                TiXmlElement* resultTable = new TiXmlElement("table");
                resultTable->SetAttribute("class", "detailsListingTable");
                missingRefFilesContainer->LinkEndChild(resultTable);

                // table body
                for (std::vector<RegressionTestDataset>::const_iterator datasetIter = testCase.missingRefDatasets_.begin();
                     datasetIter != testCase.missingRefDatasets_.end(); ++datasetIter)
                {
                    const RegressionTestDataset& refDataset = *datasetIter;
                    TiXmlElement* tr = new TiXmlElement("tr");
                    resultTable->LinkEndChild(tr);

                    // reference files (always available)
                    TiXmlElement* td = new TiXmlElement("td");
                    tr->LinkEndChild(td);

                    for (size_t i=0; i<refDataset.files_.size(); i++) {
                        const RegressionTestFile& refFile = refDataset.files_.at(i);
                        TiXmlNode* innerNode;
                        std::string refFileHref = FileSystem::cleanupPath(testCaseRefDirRel + "/" + refFile.filename_);
                        if (refFile.fileType_ == ImageFile)
                            innerNode = createImage(refFileHref, refFile.filename_, IMG_WIDTH);
                        else
                            innerNode = new TiXmlText(refFile.filename_);
                        td->LinkEndChild(createLink(refFileHref, innerNode));

                        if (i < refDataset.files_.size()-1)
                            td->LinkEndChild(createLineBreak());
                    }
                }
            }

            // unexpected output files
            if (!testCase.unexpectedOutputDatasets_.empty()) {
                tgtAssert(testCaseRefDirRel != "", "TestCase reference dir not set");

                // listing container
                TiXmlElement* unexpectedFilesContainer = new TiXmlElement("div");
                detailsBody->LinkEndChild(unexpectedFilesContainer);
                unexpectedFilesContainer->SetAttribute("class", "detailsListingContainer");

                // title
                TiXmlElement* unexpectedFilesTitle = new TiXmlElement("div");
                unexpectedFilesContainer->LinkEndChild(unexpectedFilesTitle);
                unexpectedFilesTitle->LinkEndChild(new TiXmlText("Unexpected Output Datasets:"));
                unexpectedFilesTitle->SetAttribute("class", "detailsListingTitle failure");

                // table head
                TiXmlElement* resultTable = new TiXmlElement("table");
                unexpectedFilesContainer->LinkEndChild(resultTable);
                resultTable->SetAttribute("class", "detailsListingTable");

                // table body
                for (std::vector<RegressionTestDataset>::const_iterator datasetIter = testCase.unexpectedOutputDatasets_.begin();
                        datasetIter != testCase.unexpectedOutputDatasets_.end(); ++datasetIter)
                {
                    const RegressionTestDataset& outputDataset = *datasetIter;
                    TiXmlElement* tr = new TiXmlElement("tr");
                    resultTable->LinkEndChild(tr);

                    // output dataset (may not be available)
                    TiXmlElement* td = new TiXmlElement("td");
                    tr->LinkEndChild(td);

                    for (size_t i=0; i<outputDataset.files_.size(); i++) {
                        const RegressionTestFile& outputFile = outputDataset.files_.at(i);
                        if (testCaseOutputDirRel != "") {
                            TiXmlNode* innerNode;
                            std::string outputFileHref = FileSystem::cleanupPath(testCaseOutputDirRel + "/" + outputFile.filename_);
                            if (outputFile.fileType_ == ImageFile)
                                innerNode = createImage(outputFileHref, outputFile.filename_, IMG_WIDTH);
                            else
                                innerNode = new TiXmlText(outputFile.filename_);
                            td->LinkEndChild(createLink(outputFileHref, innerNode));
                        }
                        else
                            td->LinkEndChild(new TiXmlText(outputFile.filename_));

                        if (i < outputDataset.files_.size()-1)
                            td->LinkEndChild(createLineBreak());
                    }
                }

            }

            //
            // log files
            //
            TiXmlElement* logContainer = new TiXmlElement("div");
            detailsContainer->LinkEndChild(logContainer);
            logContainer->SetAttribute("class", "detailsLogContainer");
            std::vector<TiXmlNode*> linkNodes;
            // test file
            std::string workspaceHref;
            if (testCaseOutputDirRel != "")
                workspaceHref = testCaseOutputDirRel + "/" + FileSystem::fileName(testCase.testfile_);
            else
                workspaceHref = FileSystem::relativePath(testCase.testfile_, documentPath);
            linkNodes.push_back(createLink(workspaceHref, new TiXmlText(FileSystem::fileName(testCase.testfile_))));
            // test config file
            if (testCase.configfile_ != "") {
                std::string configFileHref;
                if (testCaseOutputDirRel != "")
                    configFileHref = testCaseOutputDirRel + "/" + FileSystem::fileName(testCase.configfile_);
                else
                    configFileHref = FileSystem::relativePath(testCase.configfile_, documentPath);
                linkNodes.push_back(createLink(configFileHref, new TiXmlText(FileSystem::fileName(testCase.configfile_))));
            }
            // html log
            if (testCase.htmlLog_ != "" && testCaseOutputDirRel != "") {
                std::string htmlLogHref = FileSystem::cleanupPath(testCaseOutputDirRel + "/" + testCase.htmlLog_);
                linkNodes.push_back(createLink(htmlLogHref, new TiXmlText(testCase.htmlLog_)));
            }
            // console log
            if (testCase.consoleLog_ != "" && testCaseOutputDirRel != "") {
                std::string consoleLogHref = FileSystem::cleanupPath(testCaseOutputDirRel + "/" + testCase.consoleLog_);
                linkNodes.push_back(createLink(consoleLogHref, new TiXmlText(testCase.consoleLog_)));
            }
            // add links to parent node
            for (size_t i=0; i<linkNodes.size(); i++) {
                logContainer->LinkEndChild(linkNodes.at(i));
                if (i<linkNodes.size()-1)
                    logContainer->LinkEndChild(new TiXmlText(" | "));
            }

            // write out document for current testcase, if in split mode
            if (splitReport) {
                saveDocument(document, testCaseHTMLFilename);
                delete document;
                document = 0;
                htmlNode = 0;
                body = 0;
            }

        } // testcase

    } // module (test detail loop)

    // write out root document, if in single file mode
    if (!splitReport) {
        tgtAssert(document, "no root document");
        saveDocument(document, htmlFilename);
        delete document;
        document = 0;
    }

    tgtAssert(!document, "root document not deleted");
}

void HTMLReportGenerator::createDocument(const std::string& filename, const std::string& title, TiXmlDocument*& document, TiXmlElement*& htmlNode) {

    document = new TiXmlDocument(filename);

    // HTML root node
    htmlNode = new TiXmlElement("html");
    document->LinkEndChild(htmlNode);
    htmlNode->SetAttribute("xmlns", "http://www.w3.org/1999/xhtml");

    // head node
    TiXmlElement* head = new TiXmlElement("head");
    htmlNode->LinkEndChild(head);

    // title node
    TiXmlElement* titleNode = new TiXmlElement("title");
    head->LinkEndChild(titleNode);
    titleNode->LinkEndChild(new TiXmlText(title));

    // encoding
    TiXmlElement* encoding = new TiXmlElement("meta");
    head->LinkEndChild(encoding);
    encoding->SetAttribute("http-equiv", "content-type");
    encoding->SetAttribute("content", "text/html; charset=iso-8859-1");

    // style node
    TiXmlElement* style = new TiXmlElement("link");
    head->LinkEndChild(style);
    style->SetAttribute("rel", "stylesheet");
    style->SetAttribute("type", "text/css");
    style->SetAttribute("href", "report.css");
}

void HTMLReportGenerator::saveDocument(const TiXmlDocument* document, const std::string& filename) {
    tgtAssert(document, "null pointer passsed");

    LDEBUGC("regressiontest.HTMLReportGenerator", "saving document: " << filename);

    // html doc type string
    const std::string DOC_TYPE_STR =
        "<?xml version=\"1.0\"?>\n"
        "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n";

    TiXmlPrinter printer;
    document->Accept(&printer);
    std::string outStr = DOC_TYPE_STR + printer.Str();

    // we do not CDATA tags in an HTML document
    outStr = strReplaceAll(outStr, "<![CDATA[", "");
    outStr = strReplaceAll(outStr, "]]>", "");

    // write output string to file
    std::fstream outFileStream(filename.c_str(), std::ios::out);
    if (!outFileStream.is_open() || outFileStream.bad())
        throw VoreenException("Failed to open HTML file '" + filename + "' for writing");
    outFileStream << outStr;
    bool success = !outFileStream.bad();
    outFileStream.close();

    if (!success)
        throw VoreenException("Failed to write HTML document to file: " + filename);
}

TiXmlElement* HTMLReportGenerator::createLink(const std::string& dest, TiXmlNode* innerNode) {
    TiXmlElement* a = new TiXmlElement("a");
    a->SetAttribute("href", dest);
    a->LinkEndChild(innerNode);
    return a;
}

TiXmlElement* HTMLReportGenerator::createImage(const std::string& src, const std::string& title,
    int width, int height)
{
    TiXmlElement* img = new TiXmlElement("img");
    img->SetAttribute("src", src);
    img->SetAttribute("alt", title);
    img->SetAttribute("title", title);
    if (width > 0)
        img->SetAttribute("width", width);
    if (height > 0)
        img->SetAttribute("height", height);

    return img;
}

TiXmlText* HTMLReportGenerator::createLineBreak() {
    TiXmlText* linebreak = new TiXmlText("<br/>");
    linebreak->SetCDATA(true);
    return linebreak;
}


//----------------------------------------------------------------------------------------

void JUnitReportGenerator::generateTestResultReport(const std::string& filename, const RegressionTestSuite& testSuite) {
    LINFOC("regressiontest.JUnitReportGenerator", "Generating JUnit test result report: " << filename);

    TiXmlDocument document(filename);

    // Insert XML declaration
    document.LinkEndChild(new TiXmlDeclaration(
        XmlSerializationConstants::XMLVERSION,
        XmlSerializationConstants::XMLENCODING,
        XmlSerializationConstants::XMLSTANDALONE));

    // Create XML root node
    TiXmlElement* root = new TiXmlElement("testsuite");
    root->SetAttribute("name", "RegressionTests");
    document.LinkEndChild(root);

    const std::vector<RegressionTestCase>& testCases = testSuite.testCases_;

    // generate testcase elements
    ConsoleReporter textReporter;
    double suiteTime = 0.f;
    for (size_t i=0; i<testCases.size(); i++) {
        const RegressionTestCase& testcase = testCases.at(i);

        TiXmlElement* testcaseElem = new TiXmlElement("testcase");

        std::string classname = testcase.moduleDir_;
        if (!testcase.group_.empty())
            classname += "." + testcase.group_;
        //classname += "." + testcase.name_;
        testcaseElem->SetAttribute("classname", classname);
        testcaseElem->SetAttribute("name", testcase.name_);
        testcaseElem->SetAttribute("group", testcase.group_);
        testcaseElem->SetAttribute("testfile", testcase.testfile_);
        testcaseElem->SetAttribute("module", testcase.moduleDir_);
        testcaseElem->SetAttribute("returncode", testcase.returnCode_);
        testcaseElem->SetDoubleAttribute("time", testcase.duration_);
        suiteTime += testcase.duration_;

        // failure/error elem
        if (testcase.result_ == RegressionTestError) {
            TiXmlElement* errorElem = new TiXmlElement("error");
            errorElem->LinkEndChild(new TiXmlText(testcase.errorMsg_));
            testcaseElem->LinkEndChild(errorElem);
        }
        else if (testcase.result_ == RegressionTestFailure) {
            TiXmlElement* failureElem = new TiXmlElement("failure");
            failureElem->LinkEndChild(new TiXmlText(textReporter.generateFileComparisonReport(testcase, true)));
            testcaseElem->LinkEndChild(failureElem);
        }
        else if (testcase.result_ == RegressionTestIgnored) {
            TiXmlElement* skippedElem = new TiXmlElement("skipped");
            skippedElem->LinkEndChild(new TiXmlText("test result is ignored"));
            testcaseElem->LinkEndChild(skippedElem);
        }
        else if (testcase.result_ == RegressionTestSkipped) {
            TiXmlElement* skippedElem = new TiXmlElement("skipped");
            skippedElem->LinkEndChild(new TiXmlText("test is skipped"));
            testcaseElem->LinkEndChild(skippedElem);
        }
        else if (testcase.result_ == RegressionTestSuccess) {
            /*TiXmlElement* successElem = new TiXmlElement("success");
            successElem->LinkEndChild(new TiXmlText(generateFileComparisonReport(testcase, true)));
            testcaseElem->LinkEndChild(successElem); */
        }
        else {
            tgtAssert(false, "unknown test result"); //< should never get here
        }

        root->LinkEndChild(testcaseElem);
    }
    root->SetDoubleAttribute("time", suiteTime);

    // write XML document to file
    std::fstream outfstr(filename.c_str(), std::ios::out);
    if (!outfstr.is_open() || outfstr.bad())
        throw VoreenException("Failed to open JUnit report file '" + filename + "' for writing");
    TiXmlPrinter printer;
    document.Accept(&printer);
    outfstr << printer.Str();

    bool success = !outfstr.bad();
    outfstr.close();

    if (!success)
        throw VoreenException("Failed to write JUnit report to file '" + filename + "'");
}

void JUnitReportGenerator::generateTestCoverageReport(const std::string& filename, const std::vector<VoreenModule*>& modules,
    const std::map<Processor*, std::vector<RegressionTestCase> >& coverageMap,
    const std::map<VoreenModule*, std::vector<Processor*> >& ignoreLists) {
    LINFOC("regressiontest.JUnitReportGenerator", "Generating JUnit test coverage report: " << filename);

    TiXmlDocument document(filename);

    // Insert XML declaration
    document.LinkEndChild(new TiXmlDeclaration(
        XmlSerializationConstants::XMLVERSION,
        XmlSerializationConstants::XMLENCODING,
        XmlSerializationConstants::XMLSTANDALONE));

    // Create XML root node
    TiXmlElement* root = new TiXmlElement("testsuite");
    root->SetAttribute("name", "RegressionTest - Coverage");
    document.LinkEndChild(root);

    // generate one test element per processor
    for (size_t i=0; i<modules.size(); i++) {
        VoreenModule* module = modules.at(i);
        tgtAssert(module, "module is null");
        std::vector<const Processor*> processors = module->getRegisteredProcessors();
        std::vector<Processor*> ignoredProcs;
        if (ignoreLists.find(module) != ignoreLists.end())
            ignoredProcs = ignoreLists.find(module)->second;
        for (size_t j=0; j<processors.size(); j++) {
            Processor* processor = const_cast<Processor*>(processors.at(j));
            tgtAssert(processor, "processor is null");

            if (coverageMap.find(processor) == coverageMap.end()) {
                LERRORC("regressiontest.JUnitReportGenerator", "Processor not found in coverage map: " << processor->getClassName());
                continue;
            }

            // determine whether processor is tested
            bool procTested = false;
            const std::vector<RegressionTestCase> testCases = coverageMap.find(processor)->second;
            for (size_t k=0; k<testCases.size() && !procTested; k++) {
                if (testCases.at(k).config_.enabled_)
                    procTested = true;
            }

            // processor in ignore list?
            bool procIgnored = std::find(ignoredProcs.begin(), ignoredProcs.end(), processor) != ignoredProcs.end();

            // ignore untested processors with a code state below TESTING
            if (!procTested && (processor->getCodeState() < Processor::CODE_STATE_TESTING))
                continue;

            // generate XML element
            TiXmlElement* testcaseElem = new TiXmlElement("testcase");
            testcaseElem->SetAttribute("classname", module->getDirName());
            testcaseElem->SetAttribute("name", processor->getClassName());
            if (!procTested && !procIgnored) {
                TiXmlElement* failureElem = new TiXmlElement("failure");
                std::string message = "No (enabled) test case covers this processor ";
                message += "(Code State: ";
                if (processor->getCodeState() == Processor::CODE_STATE_STABLE)
                    message += "STABLE";
                else if (processor->getCodeState() == Processor::CODE_STATE_TESTING)
                    message += "TESTING";
                else
                    message += "unknown";
                message += ")";
                failureElem->LinkEndChild(new TiXmlText(message));
                testcaseElem->LinkEndChild(failureElem);
            }
            else if (!procTested && procIgnored) {
                TiXmlElement* skippedElem = new TiXmlElement("skipped");
                skippedElem->LinkEndChild(new TiXmlText("processor is ignored"));
                testcaseElem->LinkEndChild(skippedElem);
            }
            root->LinkEndChild(testcaseElem);
        }
    }

    // write XML document to file
    std::fstream outfstr(filename.c_str(), std::ios::out);
    if (!outfstr.is_open() || outfstr.bad())
        throw VoreenException("Failed to open JUnit report file '" + filename + "' for writing");
    TiXmlPrinter printer;
    document.Accept(&printer);
    outfstr << printer.Str();

    bool success = !outfstr.bad();
    outfstr.close();

    if (!success)
        throw VoreenException("Failed to write JUnit report to file '" + filename + "'");
}

} // namespace
