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

#include "voreen/core/utils/stringutils.h"

#include <algorithm>
#include <cstdio>
#include <iomanip>

#include "tgt/vector.h" // for tgt::round

namespace voreen {

std::string itos(int i, int stringLength /*= -1*/, char fillChar /*= '0'*/) {
    std::stringstream s;
    if(stringLength > 0)
        s << std::setw(stringLength) << std::setfill(fillChar) << i;
    else
        s << i;
    return s.str();
}

std::string itos(uint32_t i, int stringLength /*= -1*/, char fillChar /*= '0'*/) {
    std::stringstream s;
    if(stringLength > 0)
        s << std::setw(stringLength) << std::setfill(fillChar) << i;
    else
        s << i;
    return s.str();
}

std::string itos(uint64_t i, int stringLength /*= -1*/, char fillChar /*= '0'*/) {
    std::stringstream s;
    if(stringLength > 0)
        s << std::setw(stringLength) << std::setfill(fillChar) << i;
    else
        s << i;
    return s.str();
}

#ifdef __APPLE__
std::string itos(long unsigned int i, int stringLength /*= -1*/, char fillChar /*= '0'*/) {
    std::stringstream s;
    if(stringLength > 0)
        s << std::setw(stringLength) << std::setfill(fillChar) << i;
    else
        s << i;
    return s.str();
}
#endif

std::string ftos(float f, int precision /* =-1*/) {
    char buffer[50];
    if (precision > -1)
        sprintf(buffer, "%.*f", precision, f);
    else
        sprintf(buffer, "%f", f);
    return std::string(buffer);
}

std::string dtos(double d, int precision /* = -1*/) {
    char buffer[50];
    if (precision > -1)
        sprintf(buffer, "%.*f", precision, d);
    else
        sprintf(buffer, "%f", d);
    return std::string(buffer);
}

#if (!defined(_MSC_VER) || (_MSC_VER < 1600)) && !defined(__APPLE__)

int stoi(const std::string& s) {
    int result = 0;
    std::istringstream isst;
    isst.str(s);
    isst >> result;
    return result;
}

float stof(const std::string& s) {
    float result = 0.f;
    std::istringstream isst;
    isst.str(s);
    isst >> result;
    return result;
}

double stod(const std::string& s) {
    double result = 0.0;
    std::istringstream isst;
    isst.str(s);
    isst >> result;
    return result;
}

#endif

char* strToChr(const std::string& s){
    char* result = new char[s.size()+1];
#ifdef WIN32
#pragma warning(disable:4996)
#endif
    s.copy(result, s.size());
#ifdef WIN32
#pragma warning(default:4996)
#endif
    result[s.size()] = '\0';
    return result;
}


std::string strReplaceAll(const std::string& str, const std::string& from, const std::string& to) {
    std::string strConv(str);
    size_t lookHere = 0;
    size_t foundHere;
    while((foundHere = strConv.find(from, lookHere)) != std::string::npos) {
        strConv.replace(foundHere, from.size(), to);
        lookHere = foundHere + to.size();
    }
    return strConv;
}

std::string strReplaceFirst(const std::string& str, const std::string& from, const std::string& to) {
    std::string strConv(str);
    size_t foundHere;
    if((foundHere = strConv.find(from)) != std::string::npos) {
        strConv.replace(foundHere, from.size(), to);
    }
    return strConv;
}

std::string strReplaceLast(const std::string& str, const std::string& from, const std::string& to) {
    std::string strConv(str);
    size_t foundHere;
    if((foundHere = strConv.rfind(from)) != std::string::npos) {
        strConv.replace(foundHere, from.size(), to);
    }
    return strConv;
}

std::string trim(std::string str, const std::string& charlist) {
    //trim from beginning:
    while(!str.empty()) {
        char curChar = str[0];

        if(charlist.find(curChar) != std::string::npos)
            str.erase(0, 1); //delete first character
        else
            break;
    }

    //trim from end:
    while(!str.empty()) {
        char curChar = str[str.size()-1];

        if(charlist.find(curChar) != std::string::npos)
            str.erase(str.size()-1, 1); //delete last character
        else
            break;
    }

    return str;
}

std::string toLower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

std::string toUpper(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

std::vector<std::string> strSplit(const std::string& str, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(str);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> strSplit(const std::string& str, const std::string& delim) {
    std::vector<std::string> result;
    if (delim.length() == 0 || str.length() == 0) {
        result.push_back(str);
        return result;
    }

    size_t pos = 0;
    while (pos != std::string::npos) {
        size_t last = pos;
        pos = str.find(delim, last);
        result.push_back(str.substr(last, pos-last));
        if (pos != std::string::npos)
            pos += delim.size();
    }
    return result;
}

bool endsWith(const std::string& input, const std::string& ending) {
    size_t li = input.length();
    size_t le = ending.length();
    if(li >= le) {
        for(size_t i=0; i<le; i++) {
            if(ending[i] != input[li - le + i])
                return false;
        }
        return true;
    }
    else
        return false;
}

bool startsWith(const std::string& input, const std::string& start) {
    size_t li = input.length();
    size_t ls = start.length();
    if(li >= ls) {
        for (size_t i=0; i<ls; i++) {
            if (start[i] != input[i])
                return false;
        }
        return true;
    }
    else
        return false;
}

std::string formatMemorySize(uint64_t bytes) {
    // calculate GB/MB/kB with 0.1f precision
    float gb = tgt::round(static_cast<float>(bytes)/107374182.4f) / 10.f;
    float mb = tgt::round(static_cast<float>(bytes)/104857.6f) / 10.f;
    float kb = tgt::round(static_cast<float>(bytes)/102.4f) / 10.f;

    // select appropriate unit
    std::stringstream out;
    if (gb >= 1.0f)
        out << gb << " GB";
    else if (mb >= 1.0f)
        out << mb << " MB";
    else if (kb >= 1.0f)
        out << kb << " kB";
    else
        out << bytes << " bytes";

    return out.str();
}

std::string formatSpatialLength(float base_mm) {
    // calculate km/m/cm/mm/micron/nm with 0.1f precision
    float km  = tgt::round(base_mm/100000.f) / 10.f;
    float m   = tgt::round(base_mm/100.f) / 10.f;
    float cm  = tgt::round(base_mm) / 10.f;
    float mm  = tgt::round(base_mm*10.f) / 10.f;
    float mu  = tgt::round(base_mm*10000.f) / 10.f;
    float nm  = tgt::round(base_mm*10000000.f) / 10.f;

    // select appropriate unit
    std::stringstream out;
    if (km >= 1.0f)
        out << km << " km";
    else if (m >= 1.0f)
        out << m << " m";
    else if (cm >= 1.0f)
        out << cm << " cm";
    else if (mm >= 1.0f)
        out << mm << " mm";
    else if (mu >= 1.0f)
        out << mu << " µm";
    else
        out << nm << " nm";

    return out.str();
}

VRN_CORE_API std::string formatSpatialLength(const tgt::vec3& position) {
    // calculate km/m/cm/mm/micron/nm with 0.1f precision
    tgt::vec3 posKM  = tgt::round(position/100000.f) / 10.f;
    tgt::vec3 posM   = tgt::round(position/100.f) / 10.f;
    tgt::vec3 posCM  = tgt::round(position) / 10.f;
    tgt::vec3 posMM  = tgt::round(position*10.f) / 10.f;
    tgt::vec3 posMU = tgt::round(position*10000.f) / 10.f;
    tgt::vec3 posNM  = tgt::round(position*10000000.f) / 10.f;

    // select appropriate unit
    std::stringstream out;
    if (posKM.x >= 1.0f || posKM.y >= 1.0f || posKM.z > 1.0f)
        out << posKM << " km";
    else if (posM.x >= 1.0f || posM.y >= 1.0f || posM.z > 1.0f)
        out << posM << " m";
    else if (posCM.x >= 1.0f || posCM.y >= 1.0f || posCM.z > 1.0f)
        out << posCM << " cm";
    else if (posMM.x >= 1.0f || posMM.y >= 1.0f || posMM.z > 1.0f)
        out << posMM << " mm";
    else if (posMU.x >= 1.0f || posMU.y >= 1.0f || posMU.z > 1.0f)
        out << posMU << " µm";
    else
        out << posNM << " nm";

    return out.str();
}

VRN_CORE_API std::string formatTime(const size_t time) {
    // calculate min/sec/ms
    size_t min  = time / 60000; // 1000ms*60sec
    size_t sec  = (time - min*60000)/1000;
    size_t ms   = time % 1000;

    // select appropriate unit
    std::stringstream out;
    if (min >= 1) {
        out << min << ":" << itos(sec,2) << "." << itos(ms,3) << " min";
    }
    else if (sec >= 1)
        out << sec << "." << itos(ms,3) << " sec";
    else
        out << ms << " ms";

    return out.str();
}

VRN_CORE_API std::string convertNewlinesWindowsToUnix(const std::string& str) {
    return strReplaceAll(str, "\r\n", "\n");
}

VRN_CORE_API std::string convertNewlinesUnixToWindows(const std::string& str) {
    return strReplaceAll(str, "\n", "\r\n");
}

std::string getFormatFromBaseTypeAndChannels(const std::string baseType, size_t numChannels) {
    if(numChannels == 1) {
        return baseType;
    } else {
        return "Vector" + std::to_string(numChannels) + "(" + baseType + ")";
    }
}

} // namespace
