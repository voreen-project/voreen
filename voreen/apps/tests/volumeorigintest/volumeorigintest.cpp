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

#include <string>
#include <fstream>
#include <iostream>

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/voreenapplication.h"

using namespace voreen;

typedef void (*TestFunctionPointer)();

int testsNum = 0;
int successNum = 0;
int failureNum = 0;

/**
 * Throws given failureMessage if condition is @c false.
 *
 * @param condition the condition
 * @param failureMessage may thrown failure message
 *
 * @throws std::string if condition is @c false
 */
void test(const bool& condition, const std::string& failureMessage) {
    if (!condition)
        throw failureMessage;
}

/**
 * Throws modified given failure message if @c actual and @c expected are not equal(using @c ==).
 *
 * @param actual actual value for comparison
 * @param expected expected value for comparison
 * @param failureMessage may thrown failure message
 *
 * @throws std::string if @c actual and @c expected are not equal
 */
template<class T>
void test(const T& actual, const T& expected, const std::string& failureMessage) {
    std::stringstream s;
    s << failureMessage << "[actual: " << actual << ", expected: " << expected << "]";
    test(actual == expected, s.str());
}

template<>
void test(const double& actual, const double& expected, const std::string& failureMessage) {
    std::stringstream s;
    s << failureMessage << "[actual: " << actual << ", expected: " << expected << "]";
    test(actual == expected, s.str());
}

template<>
void test(const float& actual, const float& expected, const std::string& failureMessage) {
    std::stringstream s;
    s << failureMessage << "[actual: " << actual << ", expected: " << expected << "]";
    test(actual == expected, s.str());
}

/**
 * Runs the given test function and gives a status report on standard output stream.
 *
 * @note In case of a throw exception except std::string,
 *       the application terminates with error code 1.
 *
 * @param testFunction pointer to test function
 * @param testName displayed test name in status report
 */
void runTest(const TestFunctionPointer& testFunction, const std::string& testName) {
    std::cout << "Testing " << testName << "... ";

    testsNum++;
    try {
        try {
            testFunction();
        }
        catch (const SerializationException& e) {
            test(false, "SerializationException thrown: " + std::string(e.what()));
        }
    }
    catch (const std::string& failureMessage) {
        std::cout << "[failure]" << std::endl;
        std::cout << "  Reason: " << failureMessage << std::endl;;
        failureNum++;
        return;
    }
    catch (...) {
        std::cout << "[fatal]" << std::endl;
        std::cout << "  Unknown exception thrown." << std::endl;
        exit(1);
    }

    std::cout << "[success]" << std::endl;
    successNum++;
}

void testPath() {
    std::string path = "path/to/file.dat";
#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif

    VolumeURL origin(path);

    test(origin.getPath(), pathExpected, "path incorrect");
    test<std::string>(origin.getProtocol(), "", "empty protocol expected");
    test<std::string>(origin.getSearchString(), "", "empty search string expected");
    test<std::string>(origin.getSearchParameter("key1"), "", "search value not empty");
}

void testParseURL() {
    std::string protocol = "dicom";
    std::string path = "path/to/file.dat";
    std::string searchString = "key1=value1&key2=value2";
    std::string url = protocol + "://" + path + "?" + searchString;

#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif
    std::string urlExpected = protocol + "://" + pathExpected + "?" + searchString;

    VolumeURL origin(url);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchParameter("key1"), std::string("value1"), "search string incorrect");
    test(origin.getSearchParameter("key2"), std::string("value2"), "search string incorrect");
}

void testParseURLWinSeparators() {
    std::string protocol = "dicom";
    std::string path = "path\\to\\file.dat";
    std::string searchString = "key1=value1&key2=value2";
    std::string url = protocol + "://" + path + "?" + searchString;

#ifdef WIN32
    std::string pathExpected = path;
#else
    std::string pathExpected = "path/to/file.dat";
#endif
    std::string urlExpected = protocol + "://" + pathExpected + "?" + searchString;

    VolumeURL origin(url);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchParameter("key1"), std::string("value1"), "search string incorrect");
    test(origin.getSearchParameter("key2"), std::string("value2"), "search string incorrect");
}

#ifdef WIN32
void testParseURLWinNetworkPath() {
    std::string protocol = "dicom";
    std::string path = "\\\\path/to\\file.dat";
    std::string searchString = "key1=value1&key2=value2";
    std::string url = protocol + "://" + path + "?" + searchString;

    std::string pathExpected = "\\\\path\\to\\file.dat";
    std::string urlExpected = protocol + "://" + pathExpected + "?" + searchString;

    VolumeURL origin(url);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchParameter("key1"), std::string("value1"), "search string incorrect");
    test(origin.getSearchParameter("key2"), std::string("value2"), "search string incorrect");
}
#endif

void testParseProtocolPath() {
    std::string protocol = "dicom";
    std::string path = "path/to/file.dat";
    std::string url = protocol + "://" + path;

#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif
    std::string urlExpected = protocol + "://" + pathExpected;

    VolumeURL origin(url);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchString(), std::string(""), "search string not empty");
    test(origin.getSearchParameter("key1"), std::string(""), "search value not empty");
}

void testParsePathSearchString() {
    std::string path = "path/to/file.dat";
    std::string searchString = "key1=value1";
    std::string url = path + "?" + searchString;

#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif
    std::string urlExpected = pathExpected + "?" + searchString;

    VolumeURL origin(url);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), std::string(""), "empty protocol exptected");
    test(origin.getSearchString(), searchString, "search string incorrect");
}

void testParseURLSpecialChars() {
    VolumeURL origin("dicom://path/to/file.dat");
    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", "value?a b&c\\d/e=");
    origin.addSearchParameter("key3", "value3");
    std::string url = origin.getURL();

    VolumeURL originParsed(url);

    test(origin == originParsed, "origins not equal");
    test(origin.getSearchParameter("key1"), std::string("value1"), "search string incorrect");
    test(origin.getSearchParameter("key2"), std::string("value?a b&c\\d/e="), "search string incorrect");
    test(origin.getSearchParameter("key3"), std::string("value3"), "search string incorrect");
}

void testConstructURL() {
    std::string protocol = "dicom";
    std::string path = "path/to/file.dat";
    std::string searchString = "key1=value1";

#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif
    std::string urlExpected = protocol + "://" + pathExpected + "?" + searchString;

    VolumeURL origin(protocol, path, searchString);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchParameter("key1"), std::string("value1"), "search string incorrect");
}

void testConstructURLWinSeparators() {
    std::string protocol = "dicom";
    std::string path = "path\\to\\file.dat";
    std::string searchString = "key1=value1&key2=value2";

#ifdef WIN32
    std::string pathExpected = path;
#else
    std::string pathExpected = "path/to/file.dat";
#endif
    std::string urlExpected = protocol + "://" + pathExpected + "?" + searchString;

    VolumeURL origin(protocol, path, searchString);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchParameter("key1"), std::string("value1"), "search string incorrect");
    test(origin.getSearchParameter("key2"), std::string("value2"), "search string incorrect");
}

void testConstructProtocolPath() {
    std::string protocol = "dicom";
    std::string path = "path/to/file.dat";

#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif
    std::string urlExpected = protocol + "://" + pathExpected;

    VolumeURL origin(protocol, path);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), protocol, "protocol incorrect");
    test(origin.getSearchString(), std::string(""), "search string not empty");
}

void testConstructPathSearchString() {
    std::string path = "path/to/file.dat";
    std::string searchString = "key1=value1&key2=value2";

#ifdef WIN32
    std::string pathExpected = "path\\to\\file.dat";
#else
    std::string pathExpected = path;
#endif
    std::string urlExpected = pathExpected + "?" + searchString;

    VolumeURL origin("", path, searchString);

    test(origin.getURL(), urlExpected, "url incorrect");
    test(origin.getPath(), pathExpected, "path incorrect");
    test(origin.getProtocol(), std::string(""), "empty protocol exptected");
    test(origin.getSearchString(), searchString, "search string incorrect");
}

void testEquality() {
    std::string url = "dicom://path/to/file.dat?param1=value1";
    std::string url3 = "dicom://path/to/file.dat";
    std::string url4 = "path/to/file.dat";

    VolumeURL origin(url);
    VolumeURL origin2(url);
    VolumeURL origin3(url3);
    VolumeURL origin4(url4);

    test(origin == origin2, "origins not equal");
    test(!(origin == origin3), "origins falsely equal");
    test(!(origin == origin4), "origins falsely equal");
}

void testSerialization() {
    std::string url = "dicom://path/to/file.dat?param1=value1&param2=value2";
    VolumeURL origin(url);

    std::string baseDir = tgt::FileSystem::absolutePath(VoreenApplication::app()->getProgramPath()); // Forces conversion from absolute in relative path during deserialization.

    XmlSerializer serializer(baseDir);
    serializer.serialize("origin", origin);
    std::stringstream outStream;
    serializer.write(outStream);

    XmlDeserializer deserializer(baseDir);
    std::stringstream inStream;
    inStream << outStream.str();
    deserializer.read(inStream);
    VolumeURL originDeserialized;
    deserializer.deserialize("origin", originDeserialized);

    test(origin == originDeserialized, "origins not equal");
}

void testSerializationSpecialChars() {
    std::string url = "dicom://path/to/file.dat";
    VolumeURL origin(url);
    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", "value?a b&c\\d/e=");
    origin.addSearchParameter("key3", "value3");

    std::string baseDir = tgt::FileSystem::absolutePath(VoreenApplication::app()->getProgramPath()); // Forces conversion from absolute in relative path during deserialization.

    XmlSerializer serializer(baseDir);
    serializer.serialize("origin", origin);
    std::stringstream outStream;
    serializer.write(outStream);

    XmlDeserializer deserializer(baseDir);
    std::stringstream inStream;
    inStream << outStream.str();
    deserializer.read(inStream);
    VolumeURL originDeserialized;
    deserializer.deserialize("origin", originDeserialized);

    test(origin == originDeserialized, "origins not equal");
}

void testGetSearchParameter() {
    VolumeURL origin("path/to/file.dat");
    test(origin.getSearchParameter("key1"), std::string(""), "search parameter not empty");

    origin = VolumeURL("path/to/file.dat");
    origin.addSearchParameter("key1", "value1");
    test(origin.getSearchParameter("key1"), std::string("value1"), "invalid search parameter");
    test(origin.getSearchParameter("KEY1"), std::string(""), "search parameter not empty (case sensitive)");
    test(origin.getSearchParameter("key2"), std::string(""), "search parameter not empty");

    origin = VolumeURL("", "path/to/file.dat");
    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", "value2");
    test(origin.getSearchParameter("key1"), std::string("value1"), "invalid search parameter");
    test(origin.getSearchParameter("key2"), std::string("value2"), "invalid search parameter");
}

void testGetSearchParameterCaseInsensitive() {
    VolumeURL origin("", "path/to/file.dat");

    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", "value2");
    test(origin.getSearchParameter("KEY2", false), std::string("value2"), "invalid search parameter");
}

void testAddSearchParameter() {
    VolumeURL origin ("path/to/file.dat");

    origin.addSearchParameter("key1", "value1");
    test(origin.getSearchParameter("key1"), std::string("value1"), "invalid search parameter");
    test(origin.getSearchParameter("key2"), std::string(""), "search value not empty");

    origin.addSearchParameter("key2", "value2");
    test(origin.getSearchParameter("key1"), std::string("value1"), "invalid search parameter");
    test(origin.getSearchParameter("key2"), std::string("value2"), "invalid search parameter");
}

void testRemoveSearchParameter() {
    VolumeURL origin ("path/to/file.dat");
    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", "value2");

    origin.removeSearchParameter("key0");
    test(origin.getSearchParameter("key1"), std::string("value1"), "invalid search parameter");
    test(origin.getSearchParameter("key2"), std::string("value2"), "invalid search parameter");

    origin.removeSearchParameter("key1");
    test(origin.getSearchParameter("key1"), std::string(""), "search value not empty");
    test(origin.getSearchParameter("key2"), std::string("value2"), "invalid search parameter");

    origin.removeSearchParameter("key2");
    test(origin.getSearchParameter("key1"), std::string(""), "search value not empty");
    test(origin.getSearchParameter("key2"), std::string(""), "search value not empty");
}

void testReplaceSearchParameter() {
    VolumeURL origin ("path/to/file.dat");
    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", "value2");
    origin.addSearchParameter("key3", "value3");

    origin.addSearchParameter("key2", "value2_1");
    origin.addSearchParameter("key1", "value1_1");
    test(origin.getSearchParameter("key1"), std::string("value1_1"), "invalid search parameter");
    test(origin.getSearchParameter("key2"), std::string("value2_1"), "invalid search parameter");
}

void testSearchParameterSpecialChars() {
    VolumeURL origin ("path/to/file.dat");

    std::string valueStr = "value?a b&c\\d/e=";
    origin.addSearchParameter("key1", "value1");
    origin.addSearchParameter("key2", valueStr);
    origin.addSearchParameter("key3", "value3");
    test(origin.getSearchParameter("key1"), std::string("value1"), "invalid search parameter");
    test(origin.getSearchParameter("key2"), valueStr, "invalid search parameter");
    test(origin.getSearchParameter("key3"), std::string("value3"), "invalid search parameter");
}

/**
 * Runs serialization framework tests.
 */
int main(int argc, char** argv) {
    VoreenApplication app("volumeorigintest", "volumeorigintest", "Tests the VolumeURL class", argc, argv);
    app.initialize();
    std::cout << "VolumeURLTest application started..." << std::endl << std::endl;

    // URL parsing
    runTest(testPath, "file path only");
    runTest(testParseURL, "parse complete url");
    runTest(testParseURLWinSeparators, "parse complete url with windows separators");
    runTest(testParseProtocolPath, "parse protocol/path");
    runTest(testParsePathSearchString, "parse path/search string");
    runTest(testParseURLSpecialChars, "parse special chars");

    // URL constructing
    runTest(testConstructURL, "construct complete url");
    runTest(testConstructURLWinSeparators, "construct complete url with windows separators");
#ifdef WIN32
    runTest(testParseURLWinNetworkPath, "construct complete url with windows network path");
#endif
    runTest(testConstructProtocolPath, "construct protocol/path");
    runTest(testConstructPathSearchString, "construct path/search string");

    // serialization
    runTest(testEquality, "equality operator");
    runTest(testSerialization, "serialization");
    runTest(testSerializationSpecialChars, "serialization");

    // search parameter access
    runTest(testGetSearchParameter, "get search parameter");
    runTest(testGetSearchParameterCaseInsensitive, "get search parameter case insensitive");
    runTest(testAddSearchParameter, "add search parameter");
    runTest(testRemoveSearchParameter, "remove search parameter");
    runTest(testReplaceSearchParameter, "replace search parameter");
    runTest(testSearchParameterSpecialChars, "search parameter special chars");

    std::cout << std::endl << "VolumeURLTest application finished..." << std::endl;
    std::cout << std::endl << "---" << std::endl;
    std::cout << testsNum << " tests run, " << successNum << " successful and " << failureNum << " failed." << std::endl;

    app.deinitialize();

    if(successNum == testsNum)
        return 0;
    else
        exit(EXIT_FAILURE);
}
