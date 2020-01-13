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

#ifndef VRN_BRUKERVOLUMEREADER_H
#define VRN_BRUKERVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"
#include "voreen/core/utils/exception.h"

namespace voreen {

class VRN_CORE_API JCampParserException : public VoreenException {
public:
    JCampParserException(const std::string& what = "") : VoreenException(what) {}
    virtual ~JCampParserException() throw() {}
};

struct VRN_CORE_API JCampToken {
    enum Type {
        OPENING_BRACKET,
        CLOSING_BRACKET,
        COMMA,
        NUMBER,
        STRING,
        ENUM
    };

    JCampToken(Type type, std::string content) : type_(type), content_(content) {}
    std::string toString() const;

    Type getType() const { return type_; }
    bool isInteger() const;
    int getInteger() const;
    float getFloat() const;
    std::string getString() const { return content_; }
    std::string getEnum() const {
        if(type_ != ENUM)
            throw JCampParserException("Expected enum!");

        return content_;
    }

private:
    Type type_;
    std::string content_;
};

class VRN_CORE_API JCampElement {
public:
    void addLine(std::string line);
    bool parse();

    std::string getName() const { return name_; }
    size_t getNumDimensions() const { return dims_.size(); }

    void verifyNumElements(size_t expected) {
        if(data_.size() != expected)
            throw JCampParserException("Element " + getName() + ": NumElements mismatch!");
    }

    void verifyNumDimensions(size_t expected) {
        if(getNumDimensions() != expected)
            throw JCampParserException("Element " + getName() + ": Dimensions mismatch!");
    }

    void verifySize(size_t dimension, size_t expectedSize) {
        if(getSize(dimension) != expectedSize)
            throw JCampParserException("Element " + getName() + ": Dimension size mismatch!");
    }

    size_t getSize(size_t dimension) const {
        if(dimension > getNumDimensions())
            throw JCampParserException("Element " + getName() + ": Dimension index out of range!");

        return dims_[dimension];
    }

    JCampToken getData(size_t i=0) const {
        if(i<data_.size())
            return data_[i];
        else
            throw JCampParserException("Element " + getName() + ": Data index out of range!");
    }

    JCampToken getData(size_t x, size_t y) const {
        return getData(getSize(1) * x + y);
    }

    bool isPrivate() { return private_; }

protected:
    std::vector<JCampToken> tokenize(std::string data);
    std::vector<std::string> lines_;

    std::string name_;
    bool private_;
    std::vector<JCampToken> data_;

    std::vector<int> dims_;

    static const std::string loggerCat_;
};

class VRN_CORE_API JCampParser {
public:
    bool parseFile(std::string filename);

    bool isCommentLine(const std::string& line);
    bool isParameterLine(const std::string& line);

    JCampElement getElement(const std::string& key) {
        if(elements_.find(key) != elements_.end())
            return elements_[key];
        else
            throw JCampParserException("Element " + key + " not found!");
    }

    JCampElement getElement(const std::string& key) const {
        std::map<std::string, JCampElement>::const_iterator it = elements_.find(key);
        if(it != elements_.end())
            return (it->second);
            //return elements_[key];
        else
            throw JCampParserException("Element " + key + " not found!");
    }

    ///Check if an element named key exists, whether it is one number and return it
    float getFloat(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumElements(1);

        return e.getData(0).getFloat();
    }

    ///Check if an element named key exists, whether it is one number and return it
    int getInteger(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumElements(1);

        return e.getData(0).getInteger();
    }

    ///Check if an element named key exists, whether it is one enum and return it
    std::string getEnum(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumElements(1);

        return e.getData(0).getEnum();
    }

    ///Check if an element named key exists, whether it is one vec2 and return it
    tgt::vec2 getVec2(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumDimensions(1);
        e.verifySize(0, 2);

        return tgt::vec2(e.getData(0).getFloat(), e.getData(1).getFloat());
    }

    ///Check if an element named key exists, whether it is one ivec2 and return it
    tgt::ivec2 getIVec2(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumDimensions(1);
        e.verifySize(0, 2);

        return tgt::ivec2(e.getData(0).getInteger(), e.getData(1).getInteger());
    }

    ///Check if an element named key exists, whether it is one vec3 and return it
    tgt::vec3 getVec3(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumDimensions(1);
        e.verifySize(0, 3);

        return tgt::vec3(e.getData(0).getFloat(), e.getData(1).getFloat(), e.getData(3).getFloat());
    }

    ///Check if an element named key exists, whether it is one ivec3 and return it
    tgt::ivec3 getIVec3(const std::string& key) {
        JCampElement e = getElement(key);
        e.verifyNumDimensions(1);
        e.verifySize(0, 3);

        return tgt::ivec3(e.getData(0).getInteger(), e.getData(1).getInteger(), e.getData(2).getInteger());
    }

    std::vector<float> getFloatArray(const std::string& key) {
        std::vector<float> ret;

        JCampElement e = getElement(key);
        e.verifyNumDimensions(1);
        size_t numEntries = e.getSize(0);

        ret.reserve(numEntries);
        for (size_t i=0; i < numEntries; ++i) {
            JCampToken t = e.getData(static_cast<int>(i));

            ret.push_back(t.getFloat());
        }

        return ret;
    }

    bool containsElement(const std::string& key) const { return elements_.find(key) != elements_.end(); }

protected:
    std::map<std::string, JCampElement> elements_;
    static const std::string loggerCat_;
};

/**
 * Reader for <tt>2dseq</tt> volume files (Bruker ParaVision format).
 * TODO: This reader is still incomplete and largely untested.
 *
 * For some hints about the file format see
 *  http://imaging.mrc-cbu.cam.ac.uk/imaging/FormatBruker
 *  http://www.cabiatl.com/mricro/mricro/bru2anz/index.html
 *  http://www.orfeo-toolbox.org/doxygen/itkBruker2DSEQImageIO_8cxx-source.html
 *
 * Hooray, found some official docs: http://filer.case.edu/vxs33/pvman/D/Docs/
 *                                   http://filer.case.edu/vxs33/pvman/A/Docs/
 */
class VRN_CORE_API BrukerVolumeReader : public VolumeReader {
    struct FrameGroup {
        int len_;
        std::string groupId_;
        std::string groupComment_;
        int valsStart_;
        int valsCnt_;
    };

    struct Slice {
        Slice(int id) : id_(id) {}
        int id_;
        std::vector<int> groupIds_;
    };

    struct SliceCollection {
        SliceCollection(int numSlices) { slices_.insert(slices_.begin(), numSlices, Slice(-1)); }

        bool hasMatchingGroupIds(const std::vector<int>& groupIds, size_t sliceDim) {
            if(slices_.front().groupIds_.empty())
                return true;

            if(slices_.front().groupIds_.size() == groupIds.size()) {
                for(size_t i=0; i<groupIds.size(); i++) {
                    if(i == sliceDim)
                        continue;
                    if(slices_.front().groupIds_[i] != groupIds[i])
                        return false;
                }
                return true;
            }
            else
                return false;
        }

        void setSlice(size_t i, Slice s) {
            if(i < slices_.size()) {
                if(slices_[i].id_ != -1)
                    std::cout << "slice already set!";
                else
                    slices_[i] = s;
            }
            else
                std::cout << "slice index out of range!";
        }

        Slice& getSlice(size_t i) { return slices_[i]; }
        size_t getNumSlices() { return slices_.size(); }
    protected:
        tgt::ivec2 sliceResolution_;
        std::vector<Slice> slices_;
    };
public:
    BrukerVolumeReader(ProgressBar* progress = 0);

    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "BrukerVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Bruker ParaVision format"; }

    /**
     * Loads a single volume from the passed origin.
     *
     * \see VolumeReader::read
     **/
    virtual VolumeBase* read(const VolumeURL& origin);

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * \see VolumeReader::read
     **/
    virtual VolumeList* read(const std::string& url);

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * \param   url         url to load volume from
     * \param   volumeId    id to select the volume, if -1 all volumes will be selected
     **/
    virtual VolumeList* read(const std::string& url, int volumeId);

    std::vector<VolumeURL> listVolumes(const std::string& url) const;
protected:
    std::vector<SliceCollection> listSliceCollections(const std::string& url) const;
    bool mapTokenToMetaData(const std::string& voreenKey, Volume* vh, const JCampToken& t);
    bool mapMetaData(const std::string& key, const std::string& voreenKey, Volume* vh, const JCampParser& parser, size_t id, size_t numFrames);

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_BRUKERVOLUMEREADER_H
