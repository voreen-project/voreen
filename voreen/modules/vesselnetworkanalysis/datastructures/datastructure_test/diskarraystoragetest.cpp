/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include <cstdio>
#include <string>
#define TGT_ASSERT_H
#define tgtAssert(cond, text) if(!(cond)) { std::cerr << "Assertion " << #cond << " failed: " << text << std::endl; exit(1); }
#define TGT_FILESYSTEM_H
namespace tgt{
    struct FileSystem {
        static void deleteFile(const std::string& name) {
            std::remove(name.c_str());
        };
    };
}
#include "tgt/vector.h"
#include "diskarraystorage.h"

#include <random>

template<typename T>
bool testBuilder(int maxElements, int numVectors, std::function<T()> random, std::default_random_engine& generator) {
    DiskArrayStorage<T> storage("test.tmp");
    std::vector<std::vector<T>> vectors;
    std::vector<DiskArray<T>> arrays;

    std::uniform_int_distribution<int> dist(0,maxElements);

    for(int i=0; i<numVectors; ++i) {
        std::vector<T> vec;
        auto builder = storage.build();
        for(int j=0; j<dist(generator); ++j) {
            auto v = random();
            vec.push_back(v);
            builder.push(v);
        }

        arrays.push_back(std::move(builder).finalize());
        vectors.push_back(std::move(vec));
    }

    tgtAssert(vectors.size() == arrays.size(), "size does not match");

    for (size_t i = 0; i < vectors.size(); i++) {
        const auto& vec = vectors[i];
        const auto& arr = arrays[i];
        if(vec.size() != arr.size()) {
            std::cerr << "Fail!" << std::endl;
            std::cerr << "Expected: Array size " << vec.size() << std::endl;
            std::cerr << "Found   : Array size " << arr.size() << std::endl;
            assert(false);
            return false;
        }
        for(int j=0; j<vec.size(); ++j) {
            if(vec[j] != arr[j]) {
                std::cerr << "Fail at pos " << j << " in array " << i << "!" << std::endl;
                std::cerr << "Expected: " << +vec[j] << std::endl;
                std::cerr << "Found   : " << +arr[i] << std::endl;
                assert(false);
                return false;
            }
        }
    }
    return true;
}

template<typename T>
bool testStore(int maxElements, int numVectors, std::function<T()> random, std::default_random_engine& generator) {
    DiskArrayStorage<T> storage("test.tmp");
    std::vector<std::vector<T>> vectors;
    std::vector<DiskArray<T>> arrays;

    std::uniform_int_distribution<int> dist(0,maxElements);

    for(int i=0; i<numVectors; ++i) {
        std::vector<T> vec;
        for(int j=0; j<dist(generator); ++j) {
            vec.push_back(random());
        }
        arrays.push_back(storage.store(vec));
        vectors.push_back(std::move(vec));
    }

    tgtAssert(vectors.size() == arrays.size(), "size does not match");

    for (size_t i = 0; i < vectors.size(); i++) {
        const auto& vec = vectors[i];
        const auto& arr = arrays[i];
        if(vec.size() != arr.size()) {
            std::cerr << "Fail!" << std::endl;
            std::cerr << "Expected: Array size " << vec.size() << std::endl;
            std::cerr << "Found   : Array size " << arr.size() << std::endl;
            assert(false);
            return false;
        }
        for(int j=0; j<vec.size(); ++j) {
            if(vec[j] != arr[j]) {
                std::cerr << "Fail at pos " << j << " in array " << i << "!" << std::endl;
                std::cerr << "Expected: " << +vec[j] << std::endl;
                std::cerr << "Found   : " << +arr[i] << std::endl;
                assert(false);
                return false;
            }
        }
    }
    return true;
}

template<typename T>
bool test(int maxElements, int numVectors, std::function<T()> random, std::default_random_engine& generator) {
    bool res = true;
    res &= testBuilder(maxElements, numVectors, random, generator);
    res &= testStore(maxElements, numVectors, random, generator);
    return res;
}

template<typename T>
std::function<T()> randomInt(std::default_random_engine& generator) {
    return [&] () {
        return std::uniform_int_distribution<T>(
                std::numeric_limits<T>::min(),
                std::numeric_limits<T>::max()
                )(generator);
    };
}

template<typename T>
std::function<T()> randomFloat(std::default_random_engine& generator) {
    return [&] () {
        return std::uniform_real_distribution<T>(
                std::numeric_limits<T>::min(),
                std::numeric_limits<T>::max()
                )(generator);
    };
}

int main() {
    std::default_random_engine generator;
    generator.seed(0xdeadbeef);
#ifndef _MSC_VER
    test<uint8_t>(1000, 1000, randomInt<uint8_t>(generator), generator);
#endif
    test<uint16_t>(1000, 1000, randomInt<uint16_t>(generator), generator);
    test<uint32_t>(1000, 1000, randomInt<uint32_t>(generator), generator);
    test<uint64_t>(1000, 1000, randomInt<uint64_t>(generator), generator);
#ifndef _MSC_VER
    test<int8_t>(1000, 1000, randomInt<int8_t>(generator), generator);
#endif
    test<int16_t>(1000, 1000, randomInt<int16_t>(generator), generator);
    test<int32_t>(1000, 1000, randomInt<int32_t>(generator), generator);
    test<int64_t>(1000, 1000, randomInt<int64_t>(generator), generator);
    test<float>(1000, 1000, randomFloat<float>(generator), generator);
    test<double>(1000, 1000, randomFloat<double>(generator), generator);

    std::cout << "Test successful!" << std::endl;
}
