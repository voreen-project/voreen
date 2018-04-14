/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/voreenapplication.h"
#include "tgt/types.h"

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

/**
 * Tests serialization and deserialization of simple data types and pointers to simple data types.
 */
void testSimpleData() {
    bool b = true;
    char c = 'j';
    signed short ss = -1;
    unsigned short us = -2;
    signed int i = -3;
    unsigned int ui = 4;
    int64_t l = -69534;
    uint64_t ul = 69535;
    float f = 1.1f;
    double d = 1.2;

    //unsigned int* uip = &ui;
    int* ip = new int(-7);

    std::string str = "This is just a <short> string.";
    std::string strn = "This is a string\n with new lines.";
    std::string strr = "This is a string\r with carriage return.";

    std::stringstream stream;

    try {
        XmlSerializer s;
        s.serialize("b", b);
        s.serialize("c", c);
        s.serialize("ss", ss);
        s.serialize("us", us);
        s.serialize("i", i);
        s.serialize("ui", ui);
        s.serialize("l", l);
        s.serialize("ul", ul);
        s.serialize("f", f);
        s.serialize("d", d);
        s.serialize("ip", ip);
        s.serialize("str", str);
        s.serialize("strn", strn);
        s.serialize("strr", strr);
        s.write(stream);
        delete ip;
    }
    catch (...) {
        delete ip;
        throw;
    }

    // Reset all variables to default values...
    b = false;
    c = static_cast<char>(0);
    ss = 0;
    i = 0;
    ui = 0;
    l = 0;
    ul = 0;
    f = 0;
    d = 0;

    ip = 0;

    std::string dstr = "";
    std::string dstrn = "";
    std::string dstrr = "";

    XmlDeserializer ds;
    ds.read(stream);
    ds.deserialize("b", b);
    ds.deserialize("c", c);
    ds.deserialize("ss", ss);
    ds.deserialize("i", i);
    ds.deserialize("ui", ui);
    ds.deserialize("l", l);
    ds.deserialize("ul", ul);
    ds.deserialize("f", f);
    ds.deserialize("d", d);
    ds.deserialize("ip", ip);
    ds.deserialize("str", dstr);
    ds.deserialize("strn", dstrn);
    ds.deserialize("strr", dstrr);

    test(b, true, "b incorrect deserialized");
    test(c, 'j', "c incorrect deserialized");
    test(ss, static_cast<signed short>(-1), "ss incorrect deserialized");
    test(us, static_cast<unsigned short>(-2), "us incorrect deserialized");
    test(i, static_cast<signed int>(-3), "i incorrect deserialized");
    test(ui, static_cast<unsigned int>(4), "ui incorrect deserialized");
    test(l, static_cast<int64_t>(-69534), "l (int64_t) incorrect deserialized");
    test(ul, static_cast<uint64_t>(69535), "l (uint64_t) incorrect deserialized");
    test(f, 1.1f, "f incorrect deserialized");
    test(d, 1.2, "d incorrect deserialized");

    test(ip, "ip still null");

    test(*ip, -7, "incorrect content of ip");

    test(str == dstr, "str incorrect deserialized");
    test(strn == dstrn, "strn incorrect deserialized");
    test(strr == dstrr, "strr incorrect deserialized");

    delete ip;
}

/**
 * Tests serialization and deserialization of tgt types.
 */
void testTgtData() {
    tgt::vec2 v2(1.0f, 2.0f);
    tgt::vec3 v3(3.0f, 4.0f, 5.0f);
    tgt::vec4 v4(6.0f, 7.0f, 8.0f, 9.0f);

    tgt::ivec2 iv2(1, 2);
    tgt::ivec3 iv3(3, 4, 5);
    tgt::ivec4 iv4(6, 7, 8, 9);

    std::stringstream stream;

    XmlSerializer s;
    s.serialize("v2", v2);
    s.serialize("v3", v3);
    s.serialize("v4", v4);
    s.serialize("iv2", iv2);
    s.serialize("iv3", iv3);
    s.serialize("iv4", iv4);
    s.write(stream);

    // Reset all variables to default values...
    v2 = tgt::vec2(0.0f, 0.0f);
    v3 = tgt::vec3(0.0f, 0.0f, 0.0f);
    v4 = tgt::vec4(0.0f, 0.0f, 0.0f, 0.0f);

    iv2 = tgt::ivec2(0, 0);
    iv3 = tgt::ivec3(0, 0, 0);
    iv4 = tgt::ivec4(0, 0, 0, 0);

    XmlDeserializer d;
    d.read(stream);
    d.deserialize("v2", v2);
    d.deserialize("v3", v3);
    d.deserialize("v4", v4);
    d.deserialize("iv2", iv2);
    d.deserialize("iv3", iv3);
    d.deserialize("iv4", iv4);

    test(v2.x, 1.0f, "v2.x incorrect deserialized");
    test(v2.y, 2.0f, "v2.y incorrect deserialized");
    test(v3.x, 3.0f, "v3.x incorrect deserialized");
    test(v3.y, 4.0f, "v3.y incorrect deserialized");
    test(v3.z, 5.0f, "v3.z incorrect deserialized");
    test(v4.x, 6.0f, "v4.x incorrect deserialized");
    test(v4.y, 7.0f, "v4.y incorrect deserialized");
    test(v4.z, 8.0f, "v4.z incorrect deserialized");
    test(v4.w, 9.0f, "v4.w incorrect deserialized");

    test(iv2.x, 1, "iv2.x incorrect deserialized");
    test(iv2.y, 2, "iv2.y incorrect deserialized");
    test(iv3.x, 3, "iv3.x incorrect deserialized");
    test(iv3.y, 4, "iv3.y incorrect deserialized");
    test(iv3.z, 5, "iv3.z incorrect deserialized");
    test(iv4.x, 6, "iv4.x incorrect deserialized");
    test(iv4.y, 7, "iv4.y incorrect deserialized");
    test(iv4.z, 8, "iv4.z incorrect deserialized");
    test(iv4.w, 9, "iv4.w incorrect deserialized");
}

/**
 * Helper function which test serialization and deserialization of a variable and pointer
 * with different serialization and deserialization order.
 */
void testPtrVarOrder(const bool& varFirstAtSerialization, const bool& varFirstAtDeserialization) {
    int i = 1;
    int* ip = &i;

    std::stringstream stream;

    XmlSerializer s;
    if (varFirstAtSerialization) {
        s.serialize("i", i);
        s.serialize("ip", ip);
    }
    else {
        s.serialize("ip", ip);
        s.serialize("i", i);
    }
    s.write(stream);

    // Reset all variables to default values...
    i = 0;
    ip = 0;

    XmlDeserializer d;
    d.read(stream);
    if (varFirstAtDeserialization) {
        d.deserialize("i", i);
        d.deserialize("ip", ip);
    }
    else {
        d.deserialize("ip", ip);
        d.deserialize("i", i);
    }

    test(i, 1, "i incorrect deserialized");
    test(ip != 0, "ip still null");
    test(ip == &i, "ip does not point to i");
}

/**
 * Tests serialization and deserialization process when first a variable and
 * secondly a pointer is serialized and afterwards first the variable and
 * secondly the pointer is deserialized.
 */
void testVarPtrVarPtrOrder() {
    testPtrVarOrder(true, true);
}

/**
 * Tests serialization and deserialization process when first a variable and
 * secondly a pointer is serialized and afterwards first the pointer and
 * secondly the variable is deserialized.
 */
void testVarPtrPtrVarOrder() {
    testPtrVarOrder(true, false);
}

/**
 * Tests serialization and deserialization process when first a pointer and
 * secondly a variable is serialized and afterwards first the variable and
 * secondly the pointer is deserialized.
 */
void testPtrVarVarPtrOrder() {
    testPtrVarOrder(false, true);
}

/**
 * Tests serialization and deserialization process when first a pointer and
 * secondly a variable is serialized and afterwards first the pointer and
 * secondly the variable is deserialized.
 */
void testPtrVarPtrVarOrder() {
    testPtrVarOrder(false, false);
}

/**
 * @c UserDefinedData realizes the @c Serializable interface to implement
 * a de- and serializable user data class.
 *
 * @see Serializable
 */
class UserDefinedData : public Serializable {
public:
    /**
     * Contained data.
     */
    int data;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("data", data);
    }

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("data", data);
    }
};

/**
 * @c UserDefinedDataContainer realizes the @c Serializable interface for de- and serialization.
 *
 * @par
 * This class contains @c UserDefinedData object and a pointer to an @c UserDefinedData object.
 * As you can see this class contains data which itself realizes the @c Serializable
 * interface for de- and serialization.
 *
 * @see UserDefinedData
 * @see Serializable
 */
class UserDefinedDataContainer : public Serializable {
public:
    /**
     * Contained user defined data.
     */
    UserDefinedData data;

    /**
     * Pointer to user defined data.
     */
    UserDefinedData* datap;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("data", data);
        s.serialize("datap", datap);
    }

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("data", data);
        s.deserialize("datap", datap);
    }
};

/**
 * Tests serialization and deserialization of user defined data classes which realizes
 * the @c Serializable interface.
 */
void testUserDefinedData() {
    UserDefinedData d;
    d.data = 1;
    UserDefinedData* dp = &d;
    UserDefinedData* dp2 = new UserDefinedData();
    dp2->data = 2;

    UserDefinedDataContainer dc;
    dc.data.data = 3;
    dc.datap = &dc.data;

    std::stringstream stream;

    try {
        XmlSerializer s;
        s.serialize("d", d);
        s.serialize("dp", dp);
        s.serialize("dp2", dp2);
        s.serialize("dc", dc);
        s.write(stream);
        delete dp2;
    }
    catch (...) {
        delete dp2;
        throw;
    }

    // Reset all variables to default values...
    d.data = 0;
    dp = 0;
    dp2 = 0;
    dc.data.data = 0;
    dc.datap = 0;

    XmlDeserializer de;
    de.read(stream);
    de.deserialize("d", d);
    de.deserialize("dp", dp);
    de.deserialize("dp2", dp2);
    de.deserialize("dc", dc);

    test(d.data, 1, "d incorrect deserialized");

    test(dp != 0, "dp still null");
    test(dp == &d, "dp does not point to d");

    test(dp2 != 0, "dp2 still null");
    test(dp2->data, 2, "dp2 incorrect deserialized");

    delete dp2;

    test(dc.data.data, 3, "dc incorrect deserialized");
    test(dc.datap != 0, "dc.datap still null");
    test(dc.datap == &dc.data, "dc.datap does not point to dc.data");
}

/**
 * Tests serialization and deserialization of STL sequence containrs like @c std::vector.
 */
void testSequenceContainers() {
    std::vector<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);

    std::vector<uint8_t> vec_uint8;
    vec_uint8.push_back(1);
    vec_uint8.push_back(2);
    vec_uint8.push_back(3);

    std::vector<int8_t> vec_int8;
    vec_int8.push_back(1);
    vec_int8.push_back(-2);
    vec_int8.push_back(3);

    std::vector<uint16_t> vec_uint16;
    vec_uint16.push_back(1);
    vec_uint16.push_back(2);
    vec_uint16.push_back(3);

    std::vector<uint64_t> vec_uint64;
    vec_uint64.push_back(1);
    vec_uint64.push_back(2);
    vec_uint64.push_back(1<<30);

    std::vector<int> numbers;
    numbers.push_back(4);
    numbers.push_back(5);

    std::deque<int> deque;
    deque.push_back(1);
    deque.push_back(2);
    deque.push_back(3);

    std::list<int> list;
    list.push_back(1);
    list.push_back(2);
    list.push_back(3);

    std::stringstream stream;

    XmlSerializer s;
    Serializer serializer(s);
    serializer.serialize("v", v);
    serializer.serialize("vec_uint8", vec_uint8);
    serializer.serialize("vec_int8", vec_int8);
    serializer.serialize("vec_uint16", vec_uint16);
    serializer.serialize("vec_uint64", vec_uint64);
    serializer.serialize("numbers", numbers, "number");
    serializer.serialize("deque", deque);
    serializer.serialize("list", list);
    s.write(stream);

    // Reset all variables to default values...
    v.clear();
    test(v.size() == 0, "vector v is not empty");
    numbers.clear();
    test(numbers.size() == 0, "vector numbers is not empty");
    deque.clear();
    test(deque.size() == 0, "deque is not empty");
    list.clear();
    test(list.size() == 0, "list is not empty");

    XmlDeserializer d;
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("v", v);
    deserializer.deserialize("vec_uint8", vec_uint8);
    deserializer.deserialize("vec_int8", vec_int8);
    deserializer.deserialize("vec_uint16", vec_uint16);
    deserializer.deserialize("vec_uint64", vec_uint64);
    deserializer.deserialize("numbers", numbers, "number");
    deserializer.deserialize("deque", deque);
    deserializer.deserialize("list", list);

    test(v.size() == 3, "v: not all vector items deserialized");
    test(v[0], 1, "v: first item incorrect deserialized");
    test(v[1], 2, "v: second item incorrect deserialized");
    test(v[2], 3, "v: third item incorrect deserialized");

    test(vec_uint8.size() == 3, "vec_uint8: not all vector items deserialized");
    test(vec_uint8[0], (uint8_t)1, "vec_uint8: first item incorrect deserialized");
    test(vec_uint8[1], (uint8_t)2, "vec_uint8: second item incorrect deserialized");
    test(vec_uint8[2], (uint8_t)3, "vec_uint8: third item incorrect deserialized");

    test(vec_int8.size() == 3, "vec_int8: not all vector items deserialized");
    test(vec_int8[0], (int8_t)1, "vec_int8: first item incorrect deserialized");
    test(vec_int8[1], (int8_t)-2, "vec_int8: second item incorrect deserialized");
    test(vec_int8[2], (int8_t)3, "vec_int8: third item incorrect deserialized");

    test(vec_uint16.size() == 3, "vec_uint16: not all vector items deserialized");
    test(vec_uint16[0], (uint16_t)1, "vec_uint16: first item incorrect deserialized");
    test(vec_uint16[1], (uint16_t)2, "vec_uint16: second item incorrect deserialized");
    test(vec_uint16[2], (uint16_t)3, "vec_uint16: third item incorrect deserialized");

    test(vec_uint64.size() == 3, "vec_uint64: not all vector items deserialized");
    test(vec_uint64[0], (uint64_t)1, "vec_uint64: first item incorrect deserialized");
    test(vec_uint64[1], (uint64_t)2, "vec_uint64: second item incorrect deserialized");
    test(vec_uint64[2], (uint64_t)(1<<30), "vec_uint64: third item incorrect deserialized");

    test(numbers.size() == 2, "numbers: not all vector items deserialized");
    test(numbers[0], 4, "numbers: first item incorrect deserialized");
    test(numbers[1], 5, "numbers: second item incorrect deserialized");

    test(deque.size() == 3, "deque: not all deque items deserialized");
    test(deque[0], 1, "deque: first item incorrect deserialized");
    test(deque[1], 2, "deque: second item incorrect deserialized");
    test(deque[2], 3, "deque: third item incorrect deserialized");

    test(list.size() == 3, "list: not all list items deserialized");
    test(list.front(), 1, "list: first item incorrect deserialized");
    list.pop_front();
    test(list.front(), 2, "list: second item incorrect deserialized");
    list.pop_front();
    test(list.front(), 3, "list: third item incorrect deserialized");
}

/**
 * Tests serialization and deserialization of STL sequence containers with 'useAttributes' enabled
 */
void testSequenceContainersUseAttributes() {
    std::vector<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);

    std::vector<uint8_t> vec_uint8;
    vec_uint8.push_back(1);
    vec_uint8.push_back(2);
    vec_uint8.push_back(3);

    std::vector<int8_t> vec_int8;
    vec_int8.push_back(1);
    vec_int8.push_back(-2);
    vec_int8.push_back(3);

    std::vector<uint16_t> vec_uint16;
    vec_uint16.push_back(1);
    vec_uint16.push_back(2);
    vec_uint16.push_back(3);

    std::vector<uint64_t> vec_uint64;
    vec_uint64.push_back(1);
    vec_uint64.push_back(2);
    vec_uint64.push_back(1<<30);

    std::vector<int> numbers;
    numbers.push_back(4);
    numbers.push_back(5);

    std::deque<int> deque;
    deque.push_back(1);
    deque.push_back(2);
    deque.push_back(3);

    std::list<int> list;
    list.push_back(1);
    list.push_back(2);
    list.push_back(3);

    std::stringstream stream;

    XmlSerializer s;
    s.setUseAttributes(true);
    Serializer serializer(s);
    serializer.serialize("v", v);
    serializer.serialize("vec_uint8", vec_uint8);
    serializer.serialize("vec_int8", vec_int8);
    serializer.serialize("vec_uint16", vec_uint16);
    serializer.serialize("vec_uint64", vec_uint64);
    serializer.serialize("numbers", numbers, "number");
    serializer.serialize("deque", deque);
    serializer.serialize("list", list);
    s.write(stream);

    // Reset all variables to default values...
    v.clear();
    test(v.size() == 0, "vector v is not empty");
    numbers.clear();
    test(numbers.size() == 0, "vector numbers is not empty");
    deque.clear();
    test(deque.size() == 0, "deque is not empty");
    list.clear();
    test(list.size() == 0, "list is not empty");

    XmlDeserializer d;
    d.setUseAttributes(true);
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("v", v);
    deserializer.deserialize("vec_uint8", vec_uint8);
    deserializer.deserialize("vec_int8", vec_int8);
    deserializer.deserialize("vec_uint16", vec_uint16);
    deserializer.deserialize("vec_uint64", vec_uint64);
    deserializer.deserialize("numbers", numbers, "number");
    deserializer.deserialize("deque", deque);
    deserializer.deserialize("list", list);

    test(v.size() == 3, "v: not all vector items deserialized");
    test(v[0], 1, "v: first item incorrect deserialized");
    test(v[1], 2, "v: second item incorrect deserialized");
    test(v[2], 3, "v: third item incorrect deserialized");

    test(vec_uint8.size() == 3, "vec_uint8: not all vector items deserialized");
    test(vec_uint8[0], (uint8_t)1, "vec_uint8: first item incorrect deserialized");
    test(vec_uint8[1], (uint8_t)2, "vec_uint8: second item incorrect deserialized");
    test(vec_uint8[2], (uint8_t)3, "vec_uint8: third item incorrect deserialized");

    test(vec_int8.size() == 3, "vec_int8: not all vector items deserialized");
    test(vec_int8[0], (int8_t)1, "vec_int8: first item incorrect deserialized");
    test(vec_int8[1], (int8_t)-2, "vec_int8: second item incorrect deserialized");
    test(vec_int8[2], (int8_t)3, "vec_int8: third item incorrect deserialized");

    test(vec_uint16.size() == 3, "vec_uint16: not all vector items deserialized");
    test(vec_uint16[0], (uint16_t)1, "vec_uint16: first item incorrect deserialized");
    test(vec_uint16[1], (uint16_t)2, "vec_uint16: second item incorrect deserialized");
    test(vec_uint16[2], (uint16_t)3, "vec_uint16: third item incorrect deserialized");

    test(vec_uint64.size() == 3, "vec_uint64: not all vector items deserialized");
    test(vec_uint64[0], (uint64_t)1, "vec_uint64: first item incorrect deserialized");
    test(vec_uint64[1], (uint64_t)2, "vec_uint64: second item incorrect deserialized");
    test(vec_uint64[2], (uint64_t)(1<<30), "vec_uint64: third item incorrect deserialized");

    test(numbers.size() == 2, "numbers: not all vector items deserialized");
    test(numbers[0], 4, "numbers: first item incorrect deserialized");
    test(numbers[1], 5, "numbers: second item incorrect deserialized");

    test(deque.size() == 3, "deque: not all deque items deserialized");
    test(deque[0], 1, "deque: first item incorrect deserialized");
    test(deque[1], 2, "deque: second item incorrect deserialized");
    test(deque[2], 3, "deque: third item incorrect deserialized");

    test(list.size() == 3, "list: not all list items deserialized");
    test(list.front(), 1, "list: first item incorrect deserialized");
    list.pop_front();
    test(list.front(), 2, "list: second item incorrect deserialized");
    list.pop_front();
    test(list.front(), 3, "list: third item incorrect deserialized");
}

/**
 * Tests serialization and deserialization of a @c std::set.
 */
void testSet() {
    const int SETELEMENTCOUNT = 5;

    std::set<int> set;
    for (int i = 0; i < SETELEMENTCOUNT; ++i)
        set.insert(i);

    std::stringstream stream;

    XmlSerializer s;
    Serializer serializer(s);
    serializer.serialize("set", set);
    s.write(stream);

    // Reset all variables to default values...
    set.clear();
    test(set.size() == 0, "set is not empty");

    XmlDeserializer d;
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("set", set);

    test((int)set.size() == SETELEMENTCOUNT, "not all set items deserialized");

    bool deserializedValues[SETELEMENTCOUNT];
    for (int i = 0; i < SETELEMENTCOUNT; ++i)
        deserializedValues[i] = false;

    for (std::set<int>::iterator it = set.begin(); it != set.end(); ++it)
        if (*it >= 0 && *it < SETELEMENTCOUNT)
            deserializedValues[*it] = true;

    for (int i = 0; i < SETELEMENTCOUNT; ++i) {
        std::stringstream itemStream;
        itemStream << i;
        test(deserializedValues[i], "int item '" + itemStream.str() + "' not deserialized");
    }
}

/**
 * Tests serialization and deserialization of a @c std::map.
 */
void testMap() {
    std::map<int, std::string> m;
    m[1] = "one";
    m[2] = "two";
    m[3] = "three";

    std::stringstream stream;

    XmlSerializer s;
    Serializer serializer(s);
    serializer.serialize("m", m);
    s.write(stream);

    // Reset all variables to default values...
    m.clear();
    test(m.size() == 0, "map is not empty");

    XmlDeserializer d;
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("m", m);

    test(m.size() == 3, "not all map items deserialized");
    test(m[1] == "one", "first pair incorrect deserialized");
    test(m[2] == "two", "second pair incorrect deserialized");
    test(m[3] == "three", "third pair incorrect deserialized");
}

/**
 * Tests serialization and deserialization of cascaded STL containers
 * like @c std::vector or @c std::map.
 */
void testComplexSTL() {
    std::map<int, std::string> m1, m2, m3;
    m2[1] = "one";
    m3[2] = "two";
    m3[3] = "three";

    std::vector<std::map<int, std::string> > v;
    v.push_back(m1);
    v.push_back(m2);
    v.push_back(m3);

    std::stringstream stream;

    XmlSerializer s;
    Serializer serializer(s);
    serializer.serialize("v", v);
    s.write(stream);

    // Reset all variables to default values...
    v.clear();
    test(v.size() == 0, "vector is not empty");

    XmlDeserializer d;
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("v", v);

    test(v.size() == 3, "not all vector item deserialized");
    test(v[0].size() == 0, "incorrect size of first vector item");
    test(v[1].size() == 1, "incorrect size of second vector item");
    test(v[2].size() == 2, "incorrect size of third vector item");
    test(v[1][1] == "one", "first item of second vector item incorrect deserialized");
    test(v[2][2] == "two", "first item of thrid vector item incorrect deserialized");
    test(v[2][3] == "three", "second item of third vector item incorrect deserialized");
}

/**
 * @c Parent class for polymorphic class de- and serialization tests.
 *
 * @see Child
 * @see Serializable
 */
class Parent : public Serializable {
public:
    /**
     * Contained data.
     */
    int pdata;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("pdata", pdata);
    }

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("pdata", pdata);
    }
};

/**
 * @c Child of @c Parent for polymorphic class de- and serialization tests.
 *
 * @see Parent
 * @see Serializable
 */
class Child : public Parent {
public:
    /**
     * Contained data
     */
    int cdata;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        Parent::serialize(s);
        s.serialize("cdata", cdata);
    }

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        Parent::deserialize(s);
        s.deserialize("cdata", cdata);
    }
};

/**
 * @c Factory for @c Parent and @c Child creation.
 *
 * @see Parent
 * @see Child
 * @see SerializableFactory
 */
class Factory : public SerializableFactory {
public:
    /**
     * @see SerializableFactory::getTypeString
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const {
        if (type == typeid(Parent))
            return "Parent";
        else if (type == typeid(Child))
            return "Child";
        else
            return "";
    }

    /**
     * @see SerializableFactory::createType
     */
    virtual Serializable* createSerializableType(const std::string& typeString) const {
        if (typeString == "Parent")
            return new Parent();
        else if (typeString == "Child")
            return new Child();
        else
            return 0;
    }
};

/**
 * Tests serialization and deserialization of polymorphic classes.
 */
void testPolymorphism() {
    Parent p;
    p.pdata = 1;
    Child c;
    c.pdata = 2;
    c.cdata = 3;
    Parent* pp = &p;
    Parent* pc = &c;

    Parent* child1 = new Child();
    child1->pdata = 4;
    dynamic_cast<Child*>(child1)->cdata = 5;
    Child* child2 = new Child();
    child2->pdata = 6;
    child2->cdata = 7;

    Factory factory;

    std::stringstream stream;

    try {
        XmlSerializer s;
        s.registerFactory(&factory);
        Serializer serializer(s);
        serializer.serialize("p", p);
        serializer.serialize("c", c);
        serializer.serialize("pp", pp);
        serializer.serialize("pc", pc);
        serializer.serialize("child1", child1);
        serializer.serialize("child2", child2);
        s.write(stream);
        delete child1;
        delete child2;
    }
    catch (...) {
        delete child1;
        delete child2;
        throw;
    }

    // Reset all variables to default values...
    p.pdata = 0;
    c.pdata = 0;
    c.cdata = 0;
    pp = 0;
    pc = 0;
    child1 = 0;
    child2 = 0;

    XmlDeserializer d;
    d.registerFactory(&factory);
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("p", p);
    deserializer.deserialize("c", c);
    deserializer.deserialize("pp", pp);
    deserializer.deserialize("pc", pc);
    deserializer.deserialize("child1", child1);
    deserializer.deserialize("child2", child2);

    test(p.pdata, 1, "p.pdata incorrect deserialized");
    test(c.pdata, 2, "c.pdata incorrect deserialized");
    test(c.cdata, 3, "c.cdata incorrect deserialized");
    test(pp != 0, "pp still null");
    test(pp == &p, "pp does not point to p");
    test(pc != 0, "pc still null");
    test(pc == &c, "pc does not point to c");
    test(child1 != 0, "child1 still null");
    test(child1->pdata, 4, "child1.pdata incorrect deserialized");
    test(dynamic_cast<Child*>(child1) != 0, "child1 deserialized without using correct polymorphic type");
    test(dynamic_cast<Child*>(child1)->cdata, 5, "child2.cdata incorrect deserialized");
    test(child2 != 0, "child2 still null");
    test(child2->pdata, 6, "child2.pdata incorrect deserialized");
    test(child2->cdata, 7, "child2.cdata incorrect deserialized");

    delete child1;
    delete child2;
}

/**
 * @c Node class for building a graph with a cycle.
 *
 * @see Serializable
 */
class Node : public Serializable {
public:
    /**
     * Predecessor node.
     */
    Node* predecessor;

    /**
     * Successor node.
     */
    Node* successor;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("predecessor", predecessor);
        s.serialize("successor", successor);
    }

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("predecessor", predecessor);
        s.deserialize("successor", successor);
    }
};

/**
 * Tests serialization and deserialization of graph with a cycle.
 */
void testCycle() {
    Node node1;
    Node node2;
    Node node3;
    node1.successor = &node2;
    node2.successor = &node3;
    node3.successor = &node1;
    node1.predecessor = &node3;
    node2.predecessor = &node1;
    node3.predecessor = &node2;

    std::stringstream stream;

    XmlSerializer s;
    Serializer serializer(s);
    serializer.serialize("tree", node1);
    s.write(stream);

    // ATTENTION: Since every successor and predecessor is deserialized as a pointer reference,
    //            therefore new memory is allocated for node2 and node3. This means that you
    //            cannot check pointer addresses to these objects. But, you can check whether
    //            an equivalent data tree is deserialized.

    // Reset all variables to default values...
    node1.successor = 0;
    node1.predecessor = 0;

    XmlDeserializer d;
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("tree", node1);

    test(node1.successor != 0, "node1.successor still null");
    test(node1.predecessor != 0, "node1.predecessor still null");

    Node* pnode1 = &node1;
    Node* pnode2 = node1.successor;
    Node* pnode3 = node1.predecessor;

    test(pnode2->predecessor == pnode1, "node1 is not the predecessor of node2");
    test(pnode2->successor == pnode3, "node3 is not the successor of node2");
    test(pnode3->predecessor == pnode2, "node2 is not the predecessor of node3");
    test(pnode3->successor == pnode1, "node1 is not the successor of node3");

    delete pnode2;
    delete pnode3;
}

/**
 * @c Abstract class for testing serialization and deserialization
 * of pointers to abstract classes.
 */
class Abstract : public AbstractSerializable {
};

/**
 * @c Specific is a realization of the @c Abstract class, which
 * is not abstract so that instances can be created.
 */
class Specific : public Abstract {
public:
    /**
     * @see IAbstractSerializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("i", i);
    }

    /**
     * @see IAbstractSerializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("i", i);
    }

    /**
     * Data to check serialization and deserialization process.
     */
    int i;
};

class AbstractFactory : public SerializableFactory {
    /**
     * @see SerializableFactory::getTypeString
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const {
        if (type == typeid(Specific))
            return "Specific";
        else
            return "";
    }

    /**
     * @see SerializableFactory::createType
     */
    virtual Serializable* createSerializableType(const std::string& typeString) const {
        if (typeString == "Specific")
            return new Specific();
        else
            return 0;
    }
};

/**
 * Tests serialization and deserialization of pointers to abstract classes.
 */
void testIAbstractSerializable() {
    Abstract* a = new Specific();
    dynamic_cast<Specific*>(a)->i = 1;

    std::stringstream stream;

    AbstractFactory factory;

    XmlSerializer s;
    s.registerFactory(&factory);
    Serializer serializer(s);
    serializer.serialize("Abstract", a);
    s.write(stream);

    // Reset all variables to default values...
    delete a;
    a = 0;

    XmlDeserializer d;
    d.registerFactory(&factory);
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("Abstract", a);

    test(a != 0, "a still null");
    Specific* specific = dynamic_cast<Specific*>(a);
    test(specific != 0, "cast to Specific* not possible");
    test(specific->i, 1, "a incorrect deserialized");

    delete a;
}

/**
 * Tests that direct deseralization of a pointer to an abstract class leads to
 * a @c XMLSerializationMemoryAllocationException.
 */
void testMemoryAllocationException() {
    Abstract* a = new Specific();

    std::stringstream stream;

    XmlSerializer s;
    Serializer serializer(s);
    serializer.serialize("Abstract", a);
    s.write(stream);

    // Reset all variables to default values...
    delete a;
    a = 0;

    XmlDeserializer d;
    d.read(stream);
    try {
        Deserializer deserializer(d);
        deserializer.deserialize("Abstract", a);
        delete a;
        test(false, "No exception raised on abstract class memory allocation try");
    }
    catch (XmlSerializationMemoryAllocationException&) {
    }
}

/**
 * Tests that deserialization of not serialized data attempt
 * leads to an @c XmlSerializationNoSuchDataException.
 */
void testNoSuchDataException() {
    std::stringstream stream;

    XmlSerializer s;
    s.write(stream);

    XmlDeserializer d;
    d.read(stream);
    int i;
    try {
        d.deserialize("NotExistentKey", i);
        test(false, "No exception on deserialization of not existent key");
    }
    catch (XmlSerializationNoSuchDataException&) {
    }
}

/**
 * Tests serialization and deserialization using XML attributes.
 *
 * @attention This is a grey box test which uses implementation details
 *            to keep the test as short as possible. That is why
 *            you have to consider changing the test every time
 *            the implementation details have changed.
 */
void testUseAttributes() {
    int i = 1;
    std::string str = "short string";
    std::string strn = "string\nwith\nnew\nlines";
    tgt::ivec2 vec(2, 3);

    std::vector<int> v;
    v.push_back(4);
    v.push_back(5);

    std::map<int, std::string> m;
    m[6] = "six";
    m[7] = "seven";

    std::vector<tgt::ivec2> vv;
    vv.push_back(tgt::ivec2(8, 9));

    std::stringstream stream;

    XmlSerializer s;
    s.setUseAttributes(true);
    Serializer serializer(s);
    serializer.serialize("i", i);
    serializer.serialize("str", str);
    serializer.serialize("strn", strn);
    serializer.serialize("vec", vec);
    serializer.serialize("v", v);
    serializer.serialize("m", m);
    serializer.serialize("vv", vv);
    s.write(stream);

    // Reset all variables to default values...
    i = 0;
    std::string dstr = "";
    std::string dstrn = "";
    vec = tgt::ivec2(0, 0);
    v.clear();
    m.clear();
    vv.clear();

    test(v.size() == 0, "v not empty");
    test(m.size() == 0, "m not empty");
    test(vv.size() == 0, "vv not empty");

    XmlDeserializer d;
    d.setUseAttributes(true);
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.deserialize("i", i);
    deserializer.deserialize("str", dstr);
    deserializer.deserialize("strn", dstrn);
    deserializer.deserialize("vec", vec);
    deserializer.deserialize("v", v);
    deserializer.deserialize("m", m);
    deserializer.deserialize("vv", vv);

    test(i, 1, "i incorrect deserialized");
    test(str == dstr, "str incorrect deserialized");
    test(strn == dstrn, "strn incorrect deserialized");
    test(vec.x, 2, "vec.x incorrect deserialized");
    test(vec.y, 3, "vec.y incorrect deserialized");
    test(v.size() == 2, "v: not all items deserialized");
    test(v[0], 4, "v: first item incorrect deserialized");
    test(v[1], 5, "v: second item incorrect deserialized");
    test(m.size() == 2, "m: not all items deserialized");
    test(m[6] == "six", "m: first item incorrect deserialized");
    test(m[7] == "seven", "m: second item incorrect deserialized");
    test(vv.size() == 1, "vv: not all items deserialized");
    test(vv[0].x, 8, "vv: first item x component incorrect deserialized");
    test(vv[0].y, 9, "vv: first item y component incorrect deserialized");
}

/**
 * Tests serialization and deserialization using pointer content serialization
 * mixed up with using XML attributes.
 *
 * @attention This is a grey box test which uses implementation details
 *            to keep the test as short as possible. That is why
 *            you have to consider changing the test every time
 *            the implementation details have changed.
 */
void testUsePointerContentSerialization() {
    int i = 1;
    int* ip = &i;
    std::vector<int*> v;
    v.push_back(new int(2));
    std::map<std::string, int*> m;
    m["three"] = new int(3);

    std::stringstream stream;

    XmlSerializer s;
    s.setUsePointerContentSerialization(true);
    Serializer serializer(s);
    serializer.serialize("i", i);
    serializer.serialize("ip", ip);
    serializer.serialize("v", v);
    serializer.serialize("m", m, "number");
    s.setUseAttributes(true);
    serializer.serialize("ai", i);
    serializer.serialize("aip", ip);
    serializer.serialize("av", v);
    serializer.serialize("am", m, "number");
    s.write(stream);

    // Reset all variables to default values...
    i = 0;
    ip = new int(0);
    *v[0] = 0;
    *m["three"] = 0;

    int ai = 0;
    int* aip = new int(0);
    std::vector<int*> av;
    av.push_back(new int(0));
    std::map<std::string, int*> am;
    am["three"] = new int(0);

    XmlDeserializer d;
    d.read(stream);
    Deserializer deserializer(d);
    deserializer.setUsePointerContentSerialization(true);
    deserializer.deserialize("i", i);
    deserializer.deserialize("ip", ip);
    deserializer.deserialize("v", v);
    deserializer.deserialize("m", m, "number");
    d.setUseAttributes(true);
    deserializer.deserialize("ai", ai);
    deserializer.deserialize("aip", aip);
    deserializer.deserialize("av", av);
    deserializer.deserialize("am", am, "number");

    test(i, 1, "i incorrect deserialized");
    test(ip != &i, "ip is pointing to adress of i");
    test(*ip, 1, "ip incorrect deserialized");
    test(v.size() == 1, "v: incorrect size");
    test(*v[0], 2, "v: first item incorrect deserialized");
    test(m.size() == 1, "m: incorrect size");
    test(*m["three"], 3, "m: first item incorrect deserialized");

    test(ai, 1, "ai incorrect deserialized");
    test(aip != &ai, "aip is pointing to adress of ai");
    test(*aip, 1, "aip incorrect deserialized");
    test(av.size() == 1, "av: incorrect size");
    test(*av[0], 2, "av: first item incorrect deserialized");
    test(am.size() == 1, "am: incorrect size");
    test(*am["three"], 3, "am: first item incorrect deserialized");

    delete ip;
    delete v[0];
    delete m["three"];
    delete aip;
    delete av[0];
    delete am["three"];
}

/**
 * Tests serialization and deserialization of @c AbstractSerializable null pointers.
 */
void testBugAbstractSerializableNullPointerSerialization() {
    Abstract* a = 0;

    std::stringstream stream;

    AbstractFactory factory;

    XmlSerializer s;
    s.setUseAttributes(true);
    s.registerFactory(&factory);
    Serializer serializer(s);
    serializer.serialize("Abstract", a);
    s.write(stream);

    XmlDeserializer d;
    d.setUseAttributes(true);
    d.registerFactory(&factory);
    d.read(stream);
    try {
        Deserializer deserializer(d);
        deserializer.deserialize("Abstract", a);
    }
    catch (XmlSerializationMemoryAllocationException&) {
        test(false, "bug occured, since memory allocation exception is thrown for 0 pointer");
    }
}

/**
 * Runs serialization framework tests.
 */
int main(int argc, char** argv) {
    VoreenApplication app("serializertest", "serializertest", "Tests the XMLSerializer", argc, argv);
    app.initialize();
    std::cout << "voreen_core Serializer test application started..." << std::endl << std::endl;

    // Testing simple data serialization...
    runTest(testSimpleData, "simple data serialization");
    runTest(testTgtData, "tgt data serialization");

    // Testing different ordered pointer/variable serialization...
    runTest(testVarPtrVarPtrOrder, "var/ptr serialization var/ptr deserialization");
    runTest(testVarPtrPtrVarOrder, "var/ptr serialization ptr/var deserialization");
    runTest(testPtrVarVarPtrOrder, "ptr/var serialization var/ptr deserialization");
    runTest(testPtrVarPtrVarOrder, "ptr/var serialization ptr/var deserialization");

    // Testing user defined data serialization...
    runTest(testUserDefinedData, "user defined data serialization");

    // Testing STL sequence containers, like std::vector, serialization...
    runTest(testSequenceContainers, "STL sequence containers, like std::vector, serialization");

    // Testing STL sequence containers, like std::vector, serialization...
    runTest(testSequenceContainersUseAttributes, "STL sequence containers with 'useAttributes'");

    // Testing std::map serialization...
    runTest(testMap, "std::map serialization");

    // Testing std::set serialization...
    runTest(testSet, "std::set serialization");

    // Testing complex STL container serialization...
    runTest(testComplexSTL, "complex STL container serialization");

    // Testing polymorphic user definded data serialization...
    runTest(testPolymorphism, "polymorphic user definded data serialization");

    // Testing cycle serialization...
    runTest(testCycle, "cycle serialization");

    // Testing serialization using AbstractSerializable...
    runTest(testIAbstractSerializable, "serialization using AbstractSerializable");

    // Testing memory allocation exception...
    runTest(testMemoryAllocationException, "exception when deserializing pointer to abstract class");

    // Testing no such data exception...
    runTest(testNoSuchDataException, "no such data exception when deserializing not existent key");

    // Testing use of XML attributes for serialization...
    runTest(testUseAttributes, "serialization using XML attributes");

    // Testing use of pointer content serialization...
    runTest(testUsePointerContentSerialization, "serialization using pointer content serialiaztion");

    // Testing functionality of bug fixes...
    runTest(testBugAbstractSerializableNullPointerSerialization, "bug fix (serialization of AbstractSerializable 0 pointer)");

    std::cout << std::endl << "voreen_core Serializer test application finished..." << std::endl;
    std::cout << std::endl << "---" << std::endl;
    std::cout << testsNum << " tests run, " << successNum << " successful and " << failureNum << " failed." << std::endl;

    app.deinitialize();

    if(successNum == testsNum)
        return 0;
    else
        exit(EXIT_FAILURE);
}
