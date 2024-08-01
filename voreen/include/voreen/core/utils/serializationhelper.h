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

#ifndef VRN_SERIALIZATIONHELPER_H
#define VRN_SERIALIZATIONHELPER_H

#include "voreen/core/utils/exception.h"
#include "tgt/filesystem.h"

#include <string>
#include <iostream>
#include <fstream>

#ifdef __unix__
#include <unistd.h>
#endif

namespace voreen {

class SerializationHelper {
public:
    template<typename T>
    static void load(const std::string& filename, const std::string& type, T& obj) {
            // open file for reading
            std::fstream fileStream(filename.c_str(), std::ios_base::in);
            if (fileStream.fail()) {
                throw SerializationException("Failed to open " + type + " file '" + tgt::FileSystem::absolutePath(filename) + "' for reading.");
            }

            // read data stream into deserializer
            XmlDeserializer d(FileSys.dirName(filename));
            d.setUseAttributes(true);
            try {
                d.read(fileStream);
            }
            catch (SerializationException& e) {
                throw SerializationException("Failed to read serialization data stream from " + type + " file '"
                        + filename + "': " + e.what());
            }
            catch (...) {
                throw SerializationException("Failed to read serialization data stream from " + type + " file '"
                        + filename + "' (unknown exception).");
            }

            // deserialize " + type + " from data stream
            try {
                d.deserialize(type, obj);
            }
            catch (std::exception& e) {
                throw SerializationException("Deserialization from " + type + " file '" + filename + "' failed: " + e.what());
            }
            catch (...) {
                throw SerializationException("Deserialization from " + type + " file '" + filename + "' failed (unknown exception).");
            }
        }

    template<typename T>
    static void save(const std::string& filename, const std::string& type, const T& obj) {
        // serialize
        XmlSerializer s(FileSys.dirName(filename));
        s.setUseAttributes(true);
        s.serialize(type, obj);

        // write serialization data to temporary string stream
        std::ostringstream textStream;

        try {
            s.write(textStream);
            if (textStream.fail())
                throw SerializationException("Failed to write serialization data to string stream.");
        }
        catch (std::exception& e) {
            throw SerializationException("Failed to write serialization data to string stream: " + std::string(e.what()));
        }
        catch (...) {
            throw SerializationException("Failed to write serialization data to string stream (unknown exception).");
        }

        // Now we have a valid StringStream containing the serialization data.
        // => Open output file and write it to the file.
        // For added data security we write to a temporary file and afterwards move it into place
        // (which should be an atomic operation).
        const std::string tmpfilename = filename + ".tmp";
        std::fstream fileStream(tmpfilename.c_str(), std::ios_base::out);
        if (fileStream.fail())
            throw SerializationException("Failed to open file '" + tmpfilename + "' for writing.");

        try {
            fileStream << textStream.str();
        }
        catch (std::exception& e) {
            throw SerializationException("Failed to write serialization data stream to file '"
                    + tmpfilename + "': " + std::string(e.what()));
        }
        catch (...) {
            throw SerializationException("Failed to write serialization data stream to file '"
                    + tmpfilename + "' (unknown exception).");
        }
        fileStream.close();

        // Finally move the temporary file into place. It is important that this happens in-place,
        // without deleting the old file first.
        bool success;
#ifdef WIN32
        // rename() does not replace existing files on Windows, so we have to use this
        success = (MoveFileEx(tmpfilename.c_str(), filename.c_str(),
                    MOVEFILE_REPLACE_EXISTING | MOVEFILE_COPY_ALLOWED) != 0);

#else
        // atomic replace
        success = (rename(tmpfilename.c_str(), filename.c_str()) == 0);
#endif
        if (!success) {
#ifdef WIN32
            _unlink(tmpfilename.c_str()); // ignore failure here
#else
            unlink(tmpfilename.c_str()); // ignore failure here
#endif
            throw SerializationException("Failed to rename temporary file '" + tmpfilename + "' to '"
                    + filename + "'");
        }
    }
};

} // namespace

#endif // VRN_STRINGUTILS_H
