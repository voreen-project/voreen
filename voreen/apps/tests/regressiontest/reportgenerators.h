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

#ifndef VRN_REPORTGENERATORS_H
#define VRN_REPORTGENERATORS_H

#include "voreen/core/utils/regressiontest/regressiontestcase.h"

#include "voreen/core/utils/exception.h"

#include "tinyxml/tinyxml.h"

#include <string>
#include <vector>
#include <map>

namespace voreen {

class Processor;
class VoreenModule;

class ConsoleReporter {

public:
    /**
     * Generates a report stating the file comparison result (matching/mis-matching files) for a single test case.
     *
     * @param withMatchingFile specifies whether the report should contain matching files
     */
    std::string generateFileComparisonReport(const RegressionTestCase& testCase, bool includeMatchingFiles = true);

    /**
     * Generates a test coverage report for the passed test modules
     * by using the passed coverage map and ignore lists.
     */
    std::string generateTestCoverageReport(const std::vector<VoreenModule*>& modules,
        const std::map<Processor*, std::vector<RegressionTestCase> >& coverageMap,
        const std::map<VoreenModule*, std::vector<Processor*> >& ignoreLists);

};

//----------------------------------------------------------------------------------------

class HTMLReportGenerator {

public:
    void generateTestResultReport(const std::string& filename, const RegressionTestSuite& testSuite, bool splitReport = false);

private:
    /// Creates a document with doctype declaration and HTML header
    void createDocument(const std::string& filename, const std::string& title, TiXmlDocument*& document, TiXmlElement*& htmlNode);

    /// Writes the passed document to the passed filename
    void saveDocument(const TiXmlDocument* document, const std::string& filename);

    TiXmlElement* createLink(const std::string& dest, TiXmlNode* innerNode);
    TiXmlElement* createImage(const std::string& src, const std::string& title,
        int width=0, int height=0);
    TiXmlText* createLineBreak();

};

//----------------------------------------------------------------------------------------

class JUnitReportGenerator {

public:
    void generateTestResultReport(const std::string& filename, const RegressionTestSuite& testSuite);

    void generateTestCoverageReport(const std::string& filename, const std::vector<VoreenModule*>& modules,
        const std::map<Processor*, std::vector<RegressionTestCase> >& coverageMap,
        const std::map<VoreenModule*, std::vector<Processor*> >& ignoreLists);
};

} // namespace

#endif
