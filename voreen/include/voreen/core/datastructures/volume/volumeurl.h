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

#ifndef VRN_VOLUMEURL_H
#define VRN_VOLUMEURL_H

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

/**
 * A VolumeURL encapsulates a URL that specifies the location of a \em single volume.
 *
 * The structure of a origin URL is as follows:
 * \verbatim
 *     protocol://filepath?key1=value&key2=value2...
 * \endverbatim
 * where only the filepath component is obligatory. The optional protocol string
 * specifies the data type of the referenced volume. The search string consisting
 * of key/value pairs may be used to encode additional information necessary for
 * distinctly identifying the referenced volume within a container file.
 * Some examples for valid origin URLs are:
 * - path/to/myvolume.dat
 * - dat://path/to/myvolume.dat
 * - dicom://path/to/DICOMDIR?SeriesInstanceUID=1.3.12.2
 *
 * The VolumeURL's MetaDataContainer may be used to provide optional information
 * that can be presented in a user interface. The MetaDataContainer is not persisted.
 *
 */
struct VRN_CORE_API VolumeURL : public Serializable {

    VolumeURL();

    /**
     * Constructs the origin from the passed URL.
     *
     * @note Do not use this constructor, if the URL string contains search values with unescaped special chars (?&=\<space>).
     *  Use addSearchParameter() instead in this case.
     */
    VolumeURL(const std::string& URL);

    /**
     * Constructs the origin from the specified protocol string, filepath and optional search string.
     *
     * @note Do not use this constructor, if the search values contain unescaped special chars (?&=\<space>).
     *  Use addSearchParameter() instead in this case.
     */
    VolumeURL(const std::string& protocol, const std::string& filepath, const std::string& searchString = "");

    virtual ~VolumeURL();

    VolumeURL& operator=(const VolumeURL& rhs);
    VolumeURL(const VolumeURL& rhs);
    bool operator==(const VolumeURL& rhs) const;

    /// Returns the complete URL where volume is loaded from.
    std::string getURL() const;

    /// Returns the protocol portion of the URL, which specifies the data format. May be empty.
    std::string getProtocol() const;

    /// Returns the path portion of the URL, without the protocol specifier and the trailing search string.
    std::string getPath() const;

    /// Returns the file name component of the URL.
    std::string getFilename() const;

    /// Returns the search string portion of the URL. May be empty.
    std::string getSearchString() const;

    /**
     * Appends the given search parameter to the URL in the form: "key=value"
     * If the key is already present, its value is overridden.
     */
    void addSearchParameter(const std::string& key, const std::string& value);

    /**
     * Removes the passed search key from the URL. If no parameter with this key exists,
     * the call is ignored.
     */
    void removeSearchParameter(const std::string& key);

    /**
     * Returns the value corresponding to the passed key in the URL's search string,
     * or an empty string, if the key is not found.
     *
     * @param key name of the search string attribute to extract
     * @param caseSensitive if true, the key name is compared case-sensitively
     */
    std::string getSearchParameter(const std::string& key, bool caseSensitive = true) const;

    /**
     * Return the VolumeURL's MetaDataContainer, which may be used to store
     * additional information about the referenced volume that is not required for
     * distinctly identifying it.
     *
     * @note The MetaDataContainer is not serialized.
     */
    MetaDataContainer& getMetaDataContainer();

    /// @overload
    const MetaDataContainer& getMetaDataContainer() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    static std::string convertURLToRelativePath(const std::string& url, const std::string& basePath);
    static std::string convertURLToAbsolutePath(const std::string& url, const std::string& basePath);

private:
    static std::string constructURL(const std::string& protocol, const std::string& path, const std::map<std::string, std::string>& searchParameters);
    static void parseURL(const std::string& url, std::string& protocol, std::string& path, std::map<std::string, std::string>& searchParameters);

    static std::string constructSearchString(const std::map<std::string, std::string>& searchParameterMap);
    static std::map<std::string, std::string> parseSearchString(std::string searchString);

    static std::string escapeString(const std::string& str);
    static std::string unescapeString(const std::string& str);

    std::string protocol_;  //< protocol portion of the URL (may be empty)
    std::string path_;      //< path portion of the URL
    std::map<std::string, std::string> searchParameterMap_; //< search parameters as key/value pairs

    /// May contain additional meta information about the volume (not serialized).
    MetaDataContainer metaDataContainer_;

    static const std::string loggerCat_;
};



}   // namespace

#endif
