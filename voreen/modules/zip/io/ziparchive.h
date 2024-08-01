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

#ifndef VRN_ZIPARCHIVE_H
#define VRN_ZIPARCHIVE_H

#include "voreen/core/voreencoreapi.h"
#include "tgt/filesystem.h"
#include "tgt/types.h"

#include <cstring>
#include <map>

namespace voreen {

/**
 * Class for reading and writing zip files.
 * This reader is fairly simple: it can only read and write unencrypted zip files
 * using the deflate/inflate methods and depends on zlib (version 1.2.3).
 *
 * Zip64 file format is not supported either yet and it is based on the appnote.txt
 * file by PKWARE version 6.3.2 (http://www.pkware.com/documents/casestudies/APPNOTE.TXT)
 *
 * @author  Dirk Feldmann, November 2009
 */
class VRN_CORE_API ZipArchive {
public:
    enum ArchiveTarget { TARGET_DISK, TARGET_MEMORY };

public:
    ZipArchive(const std::string& archiveName, const bool autoOpen = true);

    /**
     * Dtor automatically calling <code>close()</code> and freein internal
     * resources.
     */
    ~ZipArchive();

    /**
     * Adds the file given in fileName to the archive. The file will be stored
     * in the archive when <code>save()</code> is called. Until then, the file
     * is marked as "new in archive" and any attempt to extract it by calling
     * <code>extractFile()</code> will fail.
     *
     * The file will be stored in the archive directly below its "root" and no
     * directory will be included, unless the parameter internalDirectory contains
     * a valid unix-like directory without leading and trailing slashed, e.g.
     * "data/foo/bar". Non-matching pathes (e.g. D:\data\foo\bar\) will be attempted
     * to be adjusted to match the criterion.
     *
     * If the file (including the directory) is already containing, the method
     * will return false and the file will not be added.
     *
     * @param   fileName    Name and path to the file which shall be added to the
     +                      archive.
     * @param   internalDirectory   Directory which shall contain the file within
     *                              the archive.
     * @param   replaceExistingFile Determines whether to replace a possibly existing
     *                              file in the archive by the given one (forces
     *                              adding the file)
     * @return  true if the file has been successfully added, false otherwise.
     */
    bool addFile(const std::string& fileName, const std::string& internalDirectory = "",
        const bool replaceExistingFile = false);

    bool addFile(tgt::File* const file, const std::string& internalDirectory = "",
        const bool replaceExistingFile = false);

    /**
     * Indicates whether the file containing the archive exists or not.
     */
    bool archiveExists() { return archiveExists_; }

    /**
     * Closes this archive and frees the internal resources and file handles.
     * This method is called from dtor automatically, but it should be called
     * to close archives which have been openend by <code>open()</code>.
     * After the archive has been closed, it can be re-opened via <code>open()
     * </code>.
     */
    void close();

    /**
     * Indicates whether a file of the given filename is contained within tihs
     * archive or not.
     *
     * @param   fileName name of the file of which its presence is to be verified.
     * @return  true if the file is contained, false if no such file is contained.
     */
    bool containsFile(const std::string& fileName) const;

    /**
     * Extracts the file of the given filename to the given target (memory or disk) and
     * return a handle to that file if it was extraced successfully. Otherwise the handle
     * is NULL.
     * If the file has been added to this archive, but <code>save()</code> has not
     * been called when the file is about to be extracted, the method will also fail
     * for it can only extract files which have already been saved with the archive.
     *
     * NOTE: The caller has to free the returned handle by deleting it using C++ operator
     * <code>delete</code> to avoid memory leaks, especially if the target is
     * TARGET_MEMORY!
     *
     * @param   fileName    Name of the file within the archive to extract.
     * @param   target  Destition for the file. The file can be extracted to disk
     *                  or memory.
     * @param   targetDirectory Name of the directory in which the file shall be extracted.
     *                          If string is empty, the current working directory will be used.
     * @param   keepDirectoryStructure  Determines whether the internal directory structure
     *              from the archive is used and eventually created on extracting the file
     *              or not.
     * @param   replaceExistingFile    Determines whether to replace an existing file while
     *              extracting it from the archive.
     * @return  An openend handle to the currently extracted file. The caller is reponsible
     *          for closing the handle using <code>close</code> and deleting it using C++
     *          delete operator.
     */
    tgt::File* extractFile(const std::string& fileName, const ZipArchive::ArchiveTarget target,
        const std::string& targetDirectory = "", bool keepDirectoryStructure = true,
        const bool replaceExistingFile = false);

    /**
     * Extracts all files within this archive to the directory given in dirName. The
     * archive's internal directory structure will be kept. Files already exsiting
     * in the destination directory will not be replaced and left unchanged.
     *
     * @param   dirName The name of the directory into which the files shall be
     *                  extracted.
     * @param   replaceExistingFiles    Determines whether to replace existing files while
     *              extracting from the archive.
     * @return  Returns the number of files which have been extracted.
     */
    size_t extractFilesToDirectory(const std::string& dirName,
        const bool replaceExistingFiles = false);

    /**
     * Returns the (internal) names (including possible directory names) of
     * all files which already exists within this archive or which have been
     * added but not been saved yet.
     */
    std::vector<std::string> getContainedFileNames() const;

    /**
     * Returns the number files within this archive, including those which have
     * been added but not yet saved.
     */
    size_t getNumFilesInArchive() const;

    /**
     * Indicates whether this archive is open or not.
     */
    bool isOpen() const;

    /**
     * Opens the archive. This usually done by the ctor, but the archive can be
     * closed by calling <code>close()</code> and re-opened using this method.
     *
     * @return  True if the archive was opened successfully, false otherwise.
     */
    bool open();

    /**
     * Removes the file given by fileName from the archive.
     *
     * NOTE: the archive will only be removed physically if <code>save()</code>
     * is called before closing the archive via <code>close()</code> or the dtor
     * is called! Files which were added, but not saved, will act as if they were
     * never present when removed.
     *
     * @return  True if the file was removed successfully, false otherwise. None-
     *          existing files cannot be removed and false will be returned.
     */
    bool removeFile(const std::string& fileName);

    /**
     * Saves the archive physically to disk by (re-)writing the archive given
     * in archiveName_ which is set by the ctor.
     *
     * @return  True if the archive was saved successfully, false otherwise.
     */
    bool save();

private:
    enum {
        SIZE_ZIPLOCALFILEHEADER = 30,
        SIZE_ZIPDATADESCRIPTOR = 12,
        SIZE_ZIPFILEHEADER = 46,
        SIZE_ZIPEOCDRECORD = 22
    };

    enum {
        SIGNATURE_ZIPLOCALFILEHEADER = 0x04034b50,
        SIGNATURE_ZIPFILEHEADER = 0x02014b50,
        SIGNATURE_ZIPEOCDHEADERRECORD = 0x06054b50
    };

#pragma pack(push, 2)
    struct ZipLocalFileHeader {
        uint32_t signature;             // 0x04034b50
        uint16_t versionNeeded;
        uint16_t generalPurposeFlag;    // usually 0 for our purposes
        uint16_t compressionMethod;     // 0 = none, 8 = deflate (common method)
        uint16_t lastModTime;
        uint16_t lastModDate;
        uint32_t crc32;
        uint32_t compressedSize;
        uint32_t uncompressedSize;
        uint16_t filenameLength;
        uint16_t extraFieldLength;
    };

    struct ZipFileHeader {
        uint32_t signature;             // 0x02014b50
        uint16_t versionMadeBy;
        uint16_t versionNeeded;
        uint16_t generalPurposeFlag;    // usually 0
        uint16_t compressionMethod;
        uint16_t lastModTime;
        uint16_t lastModDate;
        uint32_t crc32;
        uint32_t compressedSize;
        uint32_t uncompressedSize;
        uint16_t filenameLength;
        uint16_t extraFieldLength;
        uint16_t fileCommentLength;
        uint16_t diskNumberStart;
        uint16_t internalFileAttributes;
        uint32_t externalFileAttributes;
        uint32_t localHeaderOffset;
    };

    struct ZipEOCDHeaderRecord {    // EOCD = end of central directory
        uint32_t signature;         // 0x06054b50
        uint16_t numberOfDisk;
        uint16_t numberOfDiskWithStartOfCD;
        uint16_t numberOfEntriesInThisCD;
        uint16_t numberOfEntriesInCD;
        uint32_t sizeOfCD;
        uint32_t offsetStartCD;
        uint16_t commmentLength;
    };
#pragma pack(pop)

    struct VRN_CORE_API ArchivedFile {
        std::string fileName_;      // the name of the file
        std::string fileExtra_;     // extra field of file (usually empty)
        std::string fileComment_;   // comment on this file
        std::string extFileName_;   // external file name for new files
        tgt::File* extHandle_;   // handle for external files (alternative to file names, e.g. for mmapped files)
        bool isNewInArchive_;       // indicates whether this file was alreay in archive or not
        size_t localHeaderOffset_;  // offset of the Local File Header structure within the archive
        ZipFileHeader zipFileHader_;    // the File Header structure for the Central Directory
        ZipLocalFileHeader zipLocalFileHeader_; // the Local File Header structure for the file

        ArchivedFile()
            : fileName_("")
            , fileExtra_("")
            , fileComment_("")
            , extFileName_("")
            , extHandle_(0)
            , isNewInArchive_(true)
            , localHeaderOffset_(0)
        {}

        /**
         * Used to sort existing files according to their offset within the archive.
         */
        static bool smallerArchiveOffset(const ArchivedFile& lhs, const ArchivedFile& rhs) {
            return (lhs.localHeaderOffset_ < rhs.localHeaderOffset_);
        }
    };

private:
    tgt::File* extractUncompressedToDisk(const std::string& outFileName,
        const size_t uncompressedSize, const size_t archiveOffset);

    tgt::File* extractUncompressedToMemory(const std::string& outFileName,
        const size_t uncompressedSize, const size_t archiveOffset);

    size_t deflateToDisk(tgt::File& inFile, std::ofstream& archive, unsigned long& crc);

    tgt::File* inflateToDisk(const std::string& outFileName, const size_t compressedSize,
        const size_t uncompressedSize, size_t archiveOffset);

    tgt::File* inflateToMemory(const std::string& outFileName, size_t compressedSize,
        const size_t uncompressedSize, size_t archiveOffset);

    /**
     * Converts the directory names and separators as required for zip files. In
     * detail this means that backslashes '\' are converted to slashes '/', drive
     * letters and colons ':' are removed, leading slashes are removed, the names
     * are converted to lower case and a training slash is added.
     */
    std::string adjustDirectoryStyle(const std::string& directory) const;

    /**
     * Internal method which verifies whether the archive's file handle is open and valid.
     * This method should be called at first in every method attempting to access the
     * archive's files handle.
     *
     * @return false if the file handle is invalid, true if it is valid. If the handle
     * is invalid, the caller is supposed to cancel, because proceeding is impossible.
     */
    bool checkFileHandleValid() const;

    /**
     * Copies files which exist in the this archive to the one which is about to
     * replace it. This is necessary when saving an archive which already contains
     * a set of files but which might have been modified or have been deleted.
     */
    size_t copyExistingFiles(std::ofstream& ofs,
        std::vector<ZipArchive::ArchivedFile>& existingFiles);

    /**
     * Prepares the directory structure for the archived file. This used when extracting
     * a file which contains a directory structure. The structure is relative to the
     * destination directory in which the file is about to be extracted. Existing
     * directories will be detected.
     *
     * @param   zipFileName Complete name of the file within the archive, including
     *                      directory information and file name.
     */
    bool prepareDirectories(const std::string& zipFileName) const;

    bool readEOCDHeaderRecord(ZipArchive::ZipEOCDHeaderRecord& eocdHeaderRec);

    bool readFileHeader(ZipArchive::ZipFileHeader& fileHeader, size_t fileOffset);

    bool readLocalFileHeader(ZipArchive::ZipLocalFileHeader& localFileHeader, size_t fileOffset);

    std::string readString(size_t fileOffset, size_t numChars);

    bool readZipFile();

    bool writeCentralDirectory(std::ofstream& ofs, const std::vector<ZipArchive::ArchivedFile>& files);

    /**
     * Writes files, which are new to his archive. The new archive is supposed to be
     * opened in ofs and existing files should have already been written by a call to
     * <code>copyExisitingFiles()</code>.
     */
    size_t writeNewFiles(std::ofstream& ofs, std::vector<ZipArchive::ArchivedFile>& newFiles);

private:
    static const std::string loggerCat_;
    static const size_t MAX_BUFFER_SIZE;    /**< Controls memory consumption during (de-)compression */
    static const uint16_t ZIP_VERSION;      /**< Version of zip format this archive can understand (2.0). */

    tgt::File* archive_;                     /**< Handle to the archived if opened */
    const std::string archiveName_;     /**< The archive's names */
    bool archiveAltered_;               /**< Indicates whether the archive was altered. */
    bool archiveExists_;                /**< Indicates whether this archive's file existed. */

    typedef std::map<std::string, ArchivedFile> ArchiveMap;
    ArchiveMap files_;
};

} // namespace voreen

#endif
