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

#ifndef VRN_TRANSFUNC1DKEYS_H
#define VRN_TRANSFUNC1DKEYS_H

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "utils/transfuncmappingkey.h"

#include "tgt/vector.h"

#include <queue>

#include "voreen/core/utils/glsl.h"

namespace voreen {

/**
 * One dimensional, piece-wise linear transfer function based on key values.
 *
 * Internally, it is represented by a one-dimensional RGBA texture of type GL_UNSIGNED_BYTE.
 */
class VRN_CORE_API TransFunc1DKeys : public TransFunc1D {
public:
    /**
     * Constructor
     *
     * @param width desired width of the transfer function
     */
    TransFunc1DKeys(int width = 1024);

    /**
     * Copy constructor.
     *
     * @param tf reference to the TransFuncIntensity instance, which is copied
     */
    TransFunc1DKeys(const TransFunc1DKeys& tf);

    /**
     * Destructor - deletes the keys of the transfer function
     */
    virtual ~TransFunc1DKeys();

    virtual std::string getClassName() const { return "TransFunc1DKeys";     }
    virtual TransFunc1DKeys* create() const  { return new TransFunc1DKeys(); }
    virtual TransFunc1DKeys* clone() const;
    virtual void setMemberValuesFrom(const TransFuncBase* transfunc);
    virtual bool compareTo(const TransFuncBase& tf) const;
    virtual void reset();

    //--------------------------------------
    //  handle texture
    //--------------------------------------
protected:
    /** @override */
    virtual tgt::Vector4<GLfloat> getMappingForValueFloat(float x);
    /** @override */
    virtual tgt::Vector4<GLubyte> getMappingForValueUByte(float x);

    //--------------------------------------
    //  handle keys
    //--------------------------------------
public:
    /**
     * Returns the number of keys in this transfer function.
     *
     * @return the number of keys
     */
    int getNumKeys() const;

    /**
     * Returns the key at i-th position. Keys are sorted by their intensities in ascending order.
     * If a value outside of [0, getNumKeys()] is passed, it will be clamped to the appropriate values.
     *
     * @param i the i-th key will be returned
     * @return the pointer to the appropriate key
     */
    const TransFuncMappingKey* getKey(int i) const;

    TransFuncMappingKey* getKey(int i);

    /**
     * Returns all keys. Keys are sorted by their intensities in ascending order.
     *
     * @return a vector containing all the keys
     */
    const std::vector<TransFuncMappingKey*> getKeys() const;

    /**
    * Replaces the current keys by the passed ones.
    *
    * @param keys the new keys
    */
    void setKeys(std::vector<TransFuncMappingKey*> keys);

    /**
     * Adds a key to the property mapping function.
     * It will be automatically inserted into the correct position.
     *
     * @param key the key to be added
     */
    void addKey(TransFuncMappingKey* key);

    /**
     * Updates a key within the property mapping function.
     * Call this method when intensity of a key is changed.
     *
     * @param key the key to be updated
     */
    void updateKey(TransFuncMappingKey* key);

    /**
     * Removes the given key from the transfer function.
     * Also deletes the passed key.
     *
     * @param key the key that will be removed.
     */
    void removeKey(TransFuncMappingKey* key);

    /**
     * Removes all keys from the transfer function.
     */
    void clearKeys();

    /**
     * This method returns whether the mapping function is empty.
     * i.e., it contains no mapping keys.
     *
     * @return Is the mapping function empty?
     */
    bool isEmpty() const;
protected:
    /**
     * Calculates the average for the segment [segStart,segEnd).
     *
     * @param segStart the start of the segment
     * @param segEnd the end of the segment; the value itself is not part of the calculation
     * @return a vector containing the average value for the segment
     */
    tgt::vec4 getMeanValue(float segStart, float segEnd) const;

    /**
     * This method generates keys out of the given data.
     * The given data has to have the shape RGBA, otherwise the creation will be unsuccessful.
     * The method extrapolates the extrema of the colorchannels and puts a key in those places,
     * where the difference between neighboring entries is not linear.
     *
     * @param data an array of width*4 entries of bytes
     */
    void generateKeys(unsigned char* data);
    //--------------------------------------
    //  function helper
    //--------------------------------------
public:
    /**
     * Creates a default function.
     * Generates two keys:
     * One at intensity 0 with the color (0,0,0,0) ;
     * another one at intensity 1 with the color (255,255,255,255)
     */
    void setToStandardFunc();

    /// Returns true if the TransferFunction is the default function (@see setToStandardFunc())
    bool isStandardFunc() const;

    /*
     * Transforms this TF into a ramp.
     * (i.e., first key is at 0 (normalized) and alpha=0; last key is at 1 and has alpha=1)
     */
    void makeRamp();

    /**
     * Inverts the order of all keys.
     */
    void invertKeys();

    //--------------------------------------
    //  load and save
    //--------------------------------------
public:
    virtual const std::vector<std::string> getLoadFileFormats() const;
    virtual const std::vector<std::string> getSaveFileFormats() const;

        /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * The central entry point for loading a transfer function. The file extension is extracted
     * and based on that, explicit, private or protected load-procedures are called.
     * If there is no extension, loading will be unsuccessful.
     *
     * Currently supported extensions include:
     * tfi , lut , table , plist , jpg , png
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    virtual bool load(const std::string& filename);

    /**
     * Saves the transfer function to a file. Any data in the file will be overwritten.
     * The supported extensions include:
     * tfi, lut, png
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true, if the operation was successfull, false otherwise
     */
    virtual bool save(const std::string& filename) const;
protected:
    /**
     * Saves transfer function to a XML file. The extension of the file is tfi.
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true, if the operation was successful, false otherwise
     */
    bool saveTfi(const std::string& filename) const;

    /**
     * Saves transfer function as image using the DevIL library.
     * Any data in the file will be overwritten.
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true if the operation was successful, false otherwise
     */
    bool saveImage(const std::string& filename) const;

    /**
     * Saves the transfer function as a regular text file with 255 rows of R G B values.
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true if the operation was successful, false otherwise
     */
    bool saveLUT(const std::string& filename) const;

        /**
     * Loads a transfer function from an file with ending tfi.
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    bool loadTfi(const std::string& filename);

    /**
     * Loads a transfer function out of an ordinary image file. For this method, DevIL is required.
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    bool loadImage(const std::string& filename);

    /**
     * Loads a transfer function from a text file.
     * This format is used by Klaus Engel in his preintegration volume renderer.
     * Its just 256 rows of 4 entries each. -> RGBA
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    bool loadTextTable(const std::string& filename);

    /**
     * Loads a transfer function from an Osirix CLUT file.
     * Osirix CLUT doesn't support an alpha channel. So all values will be set opaque.
     *
     * @param filename the name of the CLUT file
     * @return true if loading succeeds, false otherwise
     */
    bool loadOsirixCLUT(const std::string& filename);

    /**
     * Loads a transfer function from a LUT used by ImageJ (http://rsbweb.nih.gov/ij/)
     * Those Lookup-Tables might come in three different kinds:
     * i)   binary LUT's coming from the National Institutes of Health (NIH)
     *      those files include a header containing additional information about the table
     * ii)  binary LUT's saved in a raw format, lacking additional data.
     * iii) a LUT in a simple text format, also with a 'missing' header
     *
     * @param filename the name of the file, which should be opened
     * @return true if loading succeeds, false otherwise
     */
    bool loadImageJ(const std::string& filename);

    /**
     * Opens a binary LUT-File. There a two possibilities for a binary file; raw and NIH
     * the NIH type includes additional data about version information and the number of colors.
     * No alpha-channel information is included in a LUT file. All values will be opaque.
     *
     * @param fileStream the already opened file stream used to extract the data
     * @param raw should the file be treated as raw data?
     * @return 256 if the load was successful, 0 otherwise
     */
    int openImageJBinary(std::ifstream& fileStream, bool raw);

    /**
     * Opens a LUT file containing textual information. Currently, two types are supported;
     * i)  256 rows , 3 columns  ; one column for each color
     * ii) 256 rows , 4 columns  ; the first column is an index from 0 to 255. The others like i)
     * Each entry should be seperated by a tabstop.
     * In both cases rows beginning with a non-integer character will be ignored.
     * No alpha-channel information is included in a the file. All values will be opaque.
     *
     * @param fileStream the already opened file stream used to extract the data
     * @return 256 if the load was successful, 0 otherwise
     */
    int openImageJText(std::ifstream& fileStream);

protected:
    std::vector<TransFuncMappingKey*> keys_; ///< internal representation of the transfer function as a set of keys

    static const std::string loggerCat_; ///< logger category
};

///------------------------------------------------------------------------------------------------

/**
 * Metadata encapsulating TransFunc1DKeys.
 */
class VRN_CORE_API TransFunc1DKeysMetaData : public TransFuncMetaDataGeneric<TransFunc1DKeys> {
public:
    TransFunc1DKeysMetaData();
    TransFunc1DKeysMetaData(TransFunc1DKeys* transfunc);
    TransFunc1DKeysMetaData(const std::vector<TransFunc1DKeys*>& transfunc);
    virtual MetaDataBase* clone() const;
    virtual MetaDataBase* create() const     { return new TransFunc1DKeysMetaData(); }
    virtual std::string getClassName() const { return "TransFunc1DKeysMetaData";     }
    virtual std::string toString() const {return "TransFunc1DKeysMetaData";}
};


} // namespace voreen

#endif // VRN_TRANSFUNC1DKEYS_H
