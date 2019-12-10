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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <string>
#include <list>
using namespace std;

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/link/propertylink.h"
#include "voreen/core/properties/link/linkevaluatorpython.h"
#include "voreen/core/properties/link/linkevaluatorfactory.h"
using namespace voreen;

typedef void (*TestFunctionPointer)();

int testsNum = 0;
int successNum = 0;
int failureNum = 0;

/*
    This function throws the failureMessage if condition is not fulfilled.
*/
void test(bool condition, std::string failureMessage) {
    if (!condition)
        throw failureMessage;
}

/*
    Uses the test-function to check whether actual and expected are equal and
    appends information to the failure message.
*/
template<class T>
void test(T actual, T expected, std::string failureMessage) {
    std::stringstream s;
    s << failureMessage << "[actual: " << actual << ", expected: " << expected << "]";
    test(actual == expected, s.str());
}

/*
    Runs the given test-function and gives some information on standard output stream.
*/
void runTest(const TestFunctionPointer& testFunction, std::string testName) {
    std::cout << "Testing " << testName << "... ";

    testsNum++;
    try {
        testFunction();
    } catch (std::string failureMessage) {
        std::cout << "[failure]" << std::endl;
        std::cout << "  Reason: " << failureMessage << std::endl;;
        failureNum++;
        return;
    } catch (...) {
        std::cout << "[fatal]" << std::endl;
        std::cout << "  Unknown exception thrown." << std::endl;
        exit(1);
    }

    std::cout << "[success]" << std::endl;
    successNum++;
}

/*
    This function tests scenario 1 described at
    https://dev.voreen.org/wiki/ProSeminarSS2009/PropertyLink;
*/
void testScenario1() {
    IntProp p1("p1", "", 0), p2("p2", "", 0), p3("p3", "", 0);
    PropertyLink l1(&p1, &p2), l2(&p2, &p3);
    p1.set(1);
    test(p2.get(), 1, "Value was not propagated to property p2");
    test(p3.get(), 1, "Value was not propagated to property p3");
}

/*
    This function tests scenario 2 described at
    https://dev.voreen.org/wiki/ProSeminarSS2009/PropertyLink;
*/
void testScenario2() {
    IntProp p1("p1", "", 0), p2("p2", "", 0);
    PropertyLink l1(&p1, &p2), l2(&p2, &p1);
    p1.set(1);
    test(p2.get(), 1, "Value was not propagated to property p2");
    p2.set(2);
    test(p1.get(), 2, "Value was not propagated to property p1");
}

/*
    This function tests scenario 3 described at
    https://dev.voreen.org/wiki/ProSeminarSS2009/PropertyLink;
*/
void testScenario3() {
    IntProp p1("p1", "", 0), p2("p2", "", 0), p3("p3", "", 0);
    PropertyLink l1(&p1, &p2), l2(&p1, &p3);
    p1.set(1);
    test(p2.get(), 1, "Value was not propagated to property p2");
    test(p3.get(), 1, "Value was not propagated to property p3");
}

/*
    This function tests scenario 4 described at
    https://dev.voreen.org/wiki/ProSeminarSS2009/PropertyLink;
*/
void testScenario4() {
    IntProp p1("p1", "", 0), p2("p2", "", 0), p3("p3", "", 0);
    PropertyLink l1(&p1, &p3), l2(&p2, &p3);
    p1.set(1);
    test(p3.get(), 1, "Value from property p1 was not propagated to property p3");
    p2.set(2);
    test(p3.get(), 2, "Value from property p2 was not propagated to property p3");
}

void testPython() {
#ifdef VRN_WITH_PYTHON
    std::list<std::string> liste0 = LinkEvaluatorFactory::ListFunctionNames();
    LinkEvaluatorFactory::RegisterLinkEvaluatorPython("testfunction", "def testfunction(a,b,c): \n    return a+b+c");
    std::list<std::string> liste1 = LinkEvaluatorFactory::ListFunctionNames();

    IntProp p1("p1", "", 0), p2("p2", "", 0);
    LinkEvaluatorBase* lep1 = LinkEvaluatorFactory::CreateLinkEvaluator("testfunction");
    PropertyLink l1(&p1, &p2, lep1);
    p1.set(1);
    test(p2.get(), 1, "Value from property p1 was not propagated to property p2");

    IntProp p3("p3", "", 0), p4("p4", "", 0);
    LinkEvaluatorBase* lep2 = LinkEvaluatorFactory::CreateLinkEvaluator("testfunction");
    PropertyLink l2(&p3, &p4, lep2);
    p3.set(1);
    test(p4.get(), 1, "Value from property p3 was not propagated to property p4");
#endif
}

/*
    This program tests the functionality of property linking in voreen_core.
*/
int main( int argc, char* argv[] ) {
    std::cout << "voreen_core Property Linking test application started..." << std::endl << std::endl;

    runTest(testScenario1, "Scenario1");
    runTest(testScenario2, "Scenario2");
    runTest(testScenario3, "Scenario3");
    runTest(testScenario4, "Scenario4");
    runTest(testPython, "Python");

    std::cout << std::endl << "voreen_core Property Linking test application finished..." << std::endl;
    std::cout << testsNum << " tests run, " << successNum << " successful and " << failureNum << " failed." << std::endl;
}
