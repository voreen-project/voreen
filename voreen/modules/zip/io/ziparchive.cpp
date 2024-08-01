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

#include "ziparchive.h"

#include <algorithm>
#include <cstring>

#include <zlib.h>

using tgt::FileSystem;
using tgt::File;
using tgt::RegularFile;
using tgt::MemoryFile;

namespace voreen {

const std::string ZipArchive::loggerCat_ = "tgt.ZipArchive";
const size_t ZipArchive::MAX_BUFFER_SIZE = 4096;
const uint16_t ZipArchive::ZIP_VERSION = 0x0014;

ZipArchive::ZipArchive(const std::string& archiveName, const bool autoOpen)
    : archive_(0)
    , archiveName_(archiveName)
    , archiveAltered_(false)
    , archiveExists_(FileSystem::fileExists(archiveName))
{
    if (autoOpen == true)
        open();
}

ZipArchive::~ZipArchive() {
    close();
}

bool ZipArchive::addFile(const std::string& fileName, const std::string& internalDirectory,
                         const bool replaceExistingFile)
{
    ArchivedFile af;
    af.extFileName_ = fileName;
    af.fileComment_ = "";
    af.fileExtra_ = "";
    af.isNewInArchive_ = true;
    af.localHeaderOffset_ = 0;
    af.fileName_ = adjustDirectoryStyle(internalDirectory);
    af.fileName_ += FileSystem::fileName(fileName);

    std::pair<ArchiveMap::iterator, bool> res = files_.insert(std::make_pair(af.fileName_, af));
    if (res.second == true)
        archiveAltered_ = true;
    else if (replaceExistingFile == true) {
        files_[af.extFileName_] = af;
        archiveAltered_ = true;
    } else
        LINFO("addFile(): file '" << af.fileName_ << "' is already present in this archive!");

    return res.second;
}

bool ZipArchive::addFile(File* const file, const std::string& internalDirectory,
                         const bool replaceExistingFile)
{
    if (file == 0) {
        LERROR("addFile(): file handle is NULL!");
        return false;
    }

    ArchivedFile af;
    af.extFileName_ = "";
    af.extHandle_ = file;
    af.fileComment_ = "";
    af.fileExtra_ = "";
    af.isNewInArchive_ = true;
    af.localHeaderOffset_ = 0;
    af.fileName_ = adjustDirectoryStyle(internalDirectory);
    af.fileName_ += FileSystem::fileName(file->getName());

    std::pair<ArchiveMap::iterator, bool> res = files_.insert(std::make_pair(af.fileName_, af));
    if (res.second == true)
        archiveAltered_ = true;
    else if (replaceExistingFile == true) {
        files_[af.fileName_] = af;
        archiveAltered_ = true;
    } else
        LINFO("addFile(): file '" << af.fileName_ << "' is already present in this archive!");

    return res.second;
}

void ZipArchive::close() {
    if ((archive_ != 0) && archive_->isOpen())
        archive_->close();
    files_.clear();
    delete archive_;
    archive_ = 0;
}

bool ZipArchive::containsFile(const std::string& fileName) const {
    return (files_.find(fileName) != files_.end());
}

File* ZipArchive::extractFile(const std::string& fileName, const ZipArchive::ArchiveTarget target,
                              const std::string& targetDirectory, const bool keepDirectoryStructure,
                              const bool replaceExistingFile)
{
    ArchiveMap::iterator it = files_.find(fileName);
    if (it == files_.end()) {
        LERROR("extractFile(): there is no such file named '" << fileName << "' to extract!");
        return 0;
    }

    ArchivedFile& af = it->second;
    if (af.isNewInArchive_ == true) {
        LERROR("extractFile(): the file '" << fileName << "' has just been added to the archive");
        LERROR(" and cannot be extracted, because the archive has not been saved yet. Call save() first.");
        return 0;
    }

    if (readLocalFileHeader(af.zipLocalFileHeader_, af.localHeaderOffset_) == false) {
        LERROR("extractFile(): could not read LocalFileHeader for file '" << fileName << "'!");
        return 0;
    }

    ZipLocalFileHeader& lfh = af.zipLocalFileHeader_;
    size_t bufferSize = (lfh.compressionMethod == 0) ? lfh.uncompressedSize : lfh.compressedSize;
    if (bufferSize == 0) {
        LERROR("File size of '" << af.fileName_ << "' is 0.");
        return 0;
    }

    if (lfh.generalPurposeFlag != 0x0000) {
        LERROR("The file " << af.fileName_ << " seems to make " <<
            "use of advanced features this reader cannot deal with");
        return 0;
    }

    size_t fileOffset = (af.localHeaderOffset_ + SIZE_ZIPLOCALFILEHEADER
        + lfh.filenameLength + lfh.extraFieldLength);
    std::string zipFileName = ((keepDirectoryStructure == true) ? af.fileName_
        : FileSystem::fileName(af.fileName_));

    if (targetDirectory.empty() == false) {
        if (targetDirectory[targetDirectory.size() - 1] != '/')
            zipFileName = targetDirectory + "/" + zipFileName;
        else
            zipFileName = targetDirectory + zipFileName;
    }

    LINFO("Extracting file " << af.fileName_ << " to " << zipFileName);
    if (FileSystem::fileExists(zipFileName)) {
        if (replaceExistingFile) {
            LWARNING("Overwriting existing file " << zipFileName);
        }
        else {
            LERROR("File " << zipFileName << " already exists");
            return 0;
        }
    }


    if ((keepDirectoryStructure == true) && (prepareDirectories(af.fileName_) == false)) {
        LERROR("Failed to prepare directroy structure");
        return 0;
    }

    switch (lfh.compressionMethod) {
        case 0:
            LDEBUG("compression method: archive (no compression)");
            switch (target) {
                case TARGET_DISK:
                    return extractUncompressedToDisk(zipFileName, lfh.uncompressedSize, fileOffset);
                case TARGET_MEMORY:
                    return extractUncompressedToMemory(zipFileName, lfh.uncompressedSize, fileOffset);
            }
            break;
        case 8:
            LDEBUG("compression method: deflate");
            switch (target) {
                case TARGET_DISK:
                    return inflateToDisk(zipFileName, lfh.compressedSize,
                        lfh.uncompressedSize, fileOffset);
                case TARGET_MEMORY:
                    return inflateToMemory(zipFileName, lfh.compressedSize,
                        lfh.uncompressedSize, fileOffset);
            }
            break;
        default:
            LERROR("unsupported compression method (code = " << lfh.compressionMethod << ")!");
            LERROR("This archiver currently supports only method 'deflate' and 'uncompressed'.");
            break;
    }
    return 0;
}

size_t ZipArchive::extractFilesToDirectory(const std::string& dirName,
                                           const bool replaceExistingFiles)
{
    std::string currentDir = FileSystem::currentDirectory();
    if (FileSystem::changeDirectory(dirName) == false) {
        LERROR("extractFilesToDirectory(): failed to locate target directory '" << dirName << "'!");
        return 0;
    }

    LINFO("extracting following files to directory '" << dirName << "'...");
    size_t counter = 0;
    for (ArchiveMap::const_iterator it = files_.begin(); it != files_.end(); ++it) {
        File* xFile = extractFile(it->second.fileName_, TARGET_DISK, "", true, replaceExistingFiles);
        if (xFile != 0) {
            xFile->close();
            delete xFile;
            ++counter;
        }
    }

    FileSystem::changeDirectory(currentDir);
    return counter;
}

std::vector<std::string> ZipArchive::getContainedFileNames() const {
    std::vector<std::string> fileNames;
    for (ArchiveMap::const_iterator it = files_.begin(); it != files_.end(); ++it)
        fileNames.push_back(it->first);

    return fileNames;
}

size_t ZipArchive::getNumFilesInArchive() const {
    return files_.size();
}

bool ZipArchive::isOpen() const {
    return (archive_ && archive_->isOpen());
}

bool ZipArchive::open() {
    if (isOpen() == true) {
        LWARNING("open(): the archive is already open!");
        return false;
    }

    archiveExists_ = FileSystem::fileExists(archiveName_);
    if (archiveExists_ == true) {
        archive_ = FileSys.open(archiveName_);
        if ((archive_ != 0) && (archive_->good() == true))
            return readZipFile();
        else {
            LERROR("open(): failed to open archive '" << archiveName_ << "' for reading!");
            return false;
        }
    }
    return true;
}

bool ZipArchive::removeFile(const std::string& fileName) {
    bool res = (files_.erase(fileName) > 0);
    if (res == true)
        archiveAltered_ = true;
    return res;
}

bool ZipArchive::save() {
    if (archiveAltered_ == false) {
        LINFO("save(): The archive has not been changed since reading it.");
        return true;
    }

    if ((archiveExists_ == true) && (checkFileHandleValid() == false))
        return false;

    std::vector<ArchivedFile> existingFiles;
    std::vector<ArchivedFile> newFiles;
    for (ArchiveMap::iterator it = files_.begin(); it != files_.end(); ++it) {
        if (it->second.isNewInArchive_ == true)
            newFiles.push_back(it->second);
        else
            existingFiles.push_back(it->second);
    }

    std::sort(existingFiles.begin(), existingFiles.end(), ArchivedFile::smallerArchiveOffset);

    // First rewrite all files which all already exist in the archive and which shall be kep
    // to the new output archive.
    //
    std::string tmpArchiveName = (archiveExists_ == true) ? (archiveName_ + ".tmp") : archiveName_;
    std::ofstream ofs(tmpArchiveName.c_str(), std::ios_base::out | std::ios_base::binary);
    if (ofs.fail() == true) {
        LERROR("save(): could not create an output file handle!");
        ofs.close();
        return false;
    }

    size_t numExistingFiles = copyExistingFiles(ofs, existingFiles);
    if (numExistingFiles != existingFiles.size()) {
        LERROR("save(): could not copy all previously existing files!");
        ofs.close();
        FileSystem::deleteFile(tmpArchiveName);
        return false;
    }

    size_t numNewFiles = writeNewFiles(ofs, newFiles);
    if (numNewFiles != newFiles.size()) {
        LERROR("save(): could not write all new files!");
        ofs.close();
        FileSystem::deleteFile(tmpArchiveName);
        return false;
    }

    // Merge the vector containing the formerly new files with the
    // vector of the existing files.
    //
    existingFiles.insert(existingFiles.end(), newFiles.begin(), newFiles.end());

    bool res = writeCentralDirectory(ofs, existingFiles);
    ofs.close();
    if (res == false) {
        LERROR("save(): failed to write EOCDHeaderRecord!");
        FileSystem::deleteFile(tmpArchiveName);
        return false;
    }

    // Close existing archive, rename tempory archive to archives file name
    // and re-open the handle
    //
    if (archiveExists_ == true) {
        close();
        res = FileSystem::deleteFile(archiveName_);
        if (res == true)
            FileSystem::renameFile(tmpArchiveName, archiveName_);
        else {
            LERROR("save(): failed to replace existing archive '" << archiveName_ << "' with new one ");
            LERROR("stored in '" << tmpArchiveName << "'!");
        }
    } else
        archiveExists_ = true;

    if (res == true) {
        archive_ = FileSys.open(archiveName_);
        if (checkFileHandleValid() == true) {
            for (size_t i = 0; i < existingFiles.size(); ++i)
                files_.insert(std::make_pair(existingFiles[i].fileName_, existingFiles[i]));
        }
    }
    return res;
}

// private methods
//

File* ZipArchive::extractUncompressedToDisk(const std::string& outFileName,
                                            const size_t uncompressedSize,
                                            const size_t archiveOffset)
{
    if (checkFileHandleValid() == false)
        return 0;

    std::ofstream ofs(outFileName.c_str(), std::ios_base::binary | std::ios_base::out);
    if (ofs.fail()) {
        LERROR("failed to open file '" << outFileName << "' for writing decompressed data!");
        ofs.close();
        return 0;
    }

    // Allocate buffers for in- and output
    //
    size_t inBufferSize = (uncompressedSize < MAX_BUFFER_SIZE) ? uncompressedSize : MAX_BUFFER_SIZE;
    char* inBuffer = 0;

    try {
        inBuffer = new char[inBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("inflateToMemory(): failed to allocate " << inBufferSize << " Bytes for input buffer!");
        return 0;
    }

    size_t readTotal = 0;
    archive_->seek(archiveOffset, File::BEGIN);
    do {
        // If the number of bytes left to read becomes smaller than the buffer's size
        // reduce the number of bytes to be read.
        //
        size_t left = uncompressedSize - readTotal;
        if (left < inBufferSize)
            inBufferSize = left;

        size_t read = archive_->read(inBuffer, inBufferSize);
        readTotal += read;
        ofs.write(inBuffer, read);
        if (ofs.fail() == true) {
            LERROR("Failed to write to output file '" << outFileName <<"'!");
            break;
        }
    } while ((archive_->eof() == false) && (readTotal < uncompressedSize));

    delete [] inBuffer;
    ofs.close();

    return new RegularFile(outFileName);
}

File* ZipArchive::extractUncompressedToMemory(const std::string& outFileName,
                                              const size_t uncompressedSize,
                                              const size_t archiveOffset)
{
    if (checkFileHandleValid() == false)
        return 0;

    // Allocate buffers for in- and output
    //
    size_t inBufferSize = uncompressedSize;
    char* inBuffer = 0;

    try {
        inBuffer = new char[inBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("inflateToMemory(): failed to allocate " << inBufferSize << " Bytes for input buffer!");
        return 0;
    }

    archive_->seek(archiveOffset, File::BEGIN);
    size_t read = archive_->read(inBuffer, inBufferSize);
    return new MemoryFile(inBuffer, read, outFileName, true);
}

size_t ZipArchive::deflateToDisk(File& inFile, std::ofstream& archive, unsigned long& crc) {
    if ((inFile.isOpen() == false) || (archive.is_open() == false)) {
        LERROR("deflateToDisk(): erroneous parameters! Handles might be closed.");
        return 0;
    }

    // Prepare z_stream structure for inflating
    //
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        LERROR("deflateToDisk(): delfateInit2 failed!");
        return 0;
    }

    size_t inBufferSize = (inFile.size() < MAX_BUFFER_SIZE) ? inFile.size() : MAX_BUFFER_SIZE;
    size_t outBufferSize = MAX_BUFFER_SIZE;
    char* inBuffer = 0;
    char* outBuffer = 0;

    try {
        inBuffer = new char[inBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("deflateToDisk(): failed to allocate " << inBufferSize << " Bytes for input buffer!");
        return 0;
    }

    try {
        outBuffer = new char[outBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("deflateToDisk(): failed to allocate " << outBufferSize << " Bytes for output buffer!");
        delete [] inBuffer;
        return 0;
    }

    // While we are reading the input file, we can also compute the CRC32...
    //
    crc = crc32(crc, Z_NULL, 0);
    size_t readTotal = 0;
    size_t writtenTotal = 0;
    bool error = false;
    inFile.seek(0, File::BEGIN);
    do {
        size_t read = inFile.read(inBuffer, inBufferSize);
        crc = crc32(crc, reinterpret_cast<Bytef*>(inBuffer), static_cast<uInt>(read));
        readTotal += read;
        strm.avail_in = static_cast<uInt>(read);
        strm.next_in = reinterpret_cast<Bytef*>(inBuffer);
        do {
            strm.avail_out = static_cast<uInt>(outBufferSize);
            strm.next_out = reinterpret_cast<Bytef*>(outBuffer);
            int res = deflate(&strm, ((inFile.eof() == true) ? Z_FINISH : Z_NO_FLUSH));
            if ((res != Z_OK) && (res != Z_STREAM_END)) {
                LERROR("deflateToDisk(): deflate returned code " << res << "!");
                error = true;
                break;
            }

            size_t write = outBufferSize - strm.avail_out;
            archive.write(outBuffer, write);
            if (archive.fail() == true) {
                LERROR("deflateToDisk(): failed to write to output archive!");
                error = true;
                break;
            } else
                writtenTotal += write;
        } while (strm.avail_out == 0);

        if (error == true)
            break;
    } while ((inFile.eof() == false) && (readTotal < inFile.size()));

    delete [] inBuffer;
    delete [] outBuffer;

    deflateEnd(&strm);
    return writtenTotal;
}

File* ZipArchive::inflateToDisk(const std::string& outFileName, const size_t compressedSize,
        const size_t uncompressedSize, size_t archiveOffset)
{
    if (checkFileHandleValid() == false)
        return 0;

    // Prepare z_stream structure for inflating
    //
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    strm.avail_out = 0;
    strm.next_out = Z_NULL;

    // Initialize inflate:
    //
    // NOTE: inflateInit2() is a badly documented macro. The second parameter ("windowBits")
    // suppresses header and CRC32 check by the inflate() function if it is negative!
    //
    if (inflateInit2(&strm, -15) != Z_OK) {
        LERROR("inflateToDisk(): call to zlib function inflateInit() failed!");
        return 0;
    }

    std::ofstream ofs(outFileName.c_str(), std::ios_base::binary | std::ios_base::out);
    if (ofs.fail()) {
        LERROR("failed to open file '" << outFileName << "' for writing decompressed data!");
        ofs.close();
        return 0;
    }

    // Allocate buffers for in- and output
    //
    size_t inBufferSize = (compressedSize < MAX_BUFFER_SIZE) ? compressedSize : MAX_BUFFER_SIZE;
    size_t outBufferSize = (uncompressedSize < MAX_BUFFER_SIZE) ? uncompressedSize : MAX_BUFFER_SIZE;
    char* inBuffer = 0;
    char* outBuffer = 0;

    try {
        inBuffer = new char[inBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("inflateToMemory(): failed to allocate " << inBufferSize << " Bytes for input buffer!");
        return 0;
    }

    try {
        outBuffer = new char[outBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("inflateToMemory(): failed to allocate " << outBufferSize << " Bytes for output buffer!");
        delete [] inBuffer;
        return 0;
    }

    size_t readTotal = 0;
    size_t writtenTotal = 0;
    archive_->seek(archiveOffset, File::BEGIN);
    do {
        // If the number of bytes left to read becomes smaller than the buffer's size
        // reduce the number of bytes to be read.
        //
        size_t left = compressedSize - readTotal;
        if (left < inBufferSize)
            inBufferSize = left;

        size_t read = archive_->read(inBuffer, inBufferSize);
        readTotal += read;
        strm.avail_in = static_cast<uInt>(read);
        strm.next_in = reinterpret_cast<Bytef*>(inBuffer);

        do {
            strm.avail_out = static_cast<uInt>(outBufferSize);
            strm.next_out = reinterpret_cast<Bytef*>(outBuffer);
            int res = inflate(&strm, Z_NO_FLUSH);
            if ((res == Z_OK) || (res == Z_STREAM_END)) {
                size_t written = outBufferSize - strm.avail_out;
                ofs.write(outBuffer, written);
                if (ofs.fail() == true) {
                    LERROR("Failed to write to output file '" << outFileName <<"'!");
                    break;
                }
                writtenTotal += written;
            } else
                break;
        } while ((strm.avail_out == 0) && (writtenTotal < uncompressedSize));
    } while ((archive_->eof() == false) && (readTotal < compressedSize));

    inflateEnd(&strm);
    delete [] inBuffer;
    delete [] outBuffer;
    ofs.close();

    return new RegularFile(outFileName);
}

File* ZipArchive::inflateToMemory(const std::string& outFileName, const size_t compressedSize,
        const size_t uncompressedSize, size_t archiveOffset)
{
    if (checkFileHandleValid() == false)
        return 0;

    // Prepare z_stream structure for inflating
    //
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    strm.avail_out = 0;
    strm.next_out = Z_NULL;

    // Initialize inflate:
    //
    // NOTE: inflateInit2() is an undocumented macro. The second parameter ("windowBits")
    // suppresses header and CRC32 check by the inflate() function if it is negative!
    //
    if (inflateInit2(&strm, -15) != Z_OK) {
        LERROR("inflateToDisk(): call to zlib function inflateInit() failed!");
        return 0;
    }

    // Allocate buffers for in- and output
    //
    size_t inBufferSize = (compressedSize < MAX_BUFFER_SIZE) ? compressedSize : MAX_BUFFER_SIZE;
    size_t outBufferSize = uncompressedSize;
    char* inBuffer = 0;
    char* outBuffer = 0;

    try {
        inBuffer = new char[inBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("inflateToMemory(): failed to allocate " << inBufferSize << " Bytes for input buffer!");
        return 0;
    }

    try {
        outBuffer = new char[outBufferSize];
    } catch (std::bad_alloc&) {
        LERROR("inflateToMemory(): failed to allocate " << outBufferSize << " Bytes for memory file!");
        delete [] inBuffer;
        return 0;
    }

    size_t readTotal = 0;
    size_t writtenTotal = 0;
    strm.avail_out = static_cast<uInt>(outBufferSize);
    strm.next_out = reinterpret_cast<Bytef*>(outBuffer);

    archive_->seek(archiveOffset, File::BEGIN);
    do {
        // If the number of bytes left to read becomes smaller than the buffer's size
        // reduce the number of bytes to be read.
        //
        size_t left = compressedSize - readTotal;
        if (left < inBufferSize)
            inBufferSize = left;

        size_t read = archive_->read(inBuffer, inBufferSize);
        readTotal += read;
        strm.avail_in = static_cast<uInt>(read);
        strm.next_in = reinterpret_cast<Bytef*>(inBuffer);

        // As the entire uncompressed file size is known, the buffer has
        // been allocated to be large enough to hold all inflated data and
        // should therefore not need to be resized.
        //
        int res = inflate(&strm, Z_NO_FLUSH);
        if ((res == Z_OK) || (res == Z_STREAM_END))
            writtenTotal += outBufferSize - strm.avail_out;
        if ((res == Z_STREAM_END) || (strm.avail_out > 0) || (writtenTotal >= uncompressedSize))
            break;

    } while ((archive_->eof() == false) && (readTotal < compressedSize));

    inflateEnd(&strm);
    delete [] inBuffer;

    return new MemoryFile(outBuffer, outBufferSize, outFileName, true);
}

// TODO: similar functions already exist in .cpp file for tgt::FileSystem within an
// anonymous namespace. Keeping these functions somewhere where they would be more
// accessible and recyclable would be preferable... (dirk)
//
std::string ZipArchive::adjustDirectoryStyle(const std::string& directory) const {
    if (directory.empty() == true)
        return "";

    std::string converted(directory);
    size_t offsetStart = 0;
    size_t numChars = directory.size();

    // Convert backslashes '\\' to slashes '/'
    //
    for (size_t i = 0; i < directory.size(); ++i) {
        if (converted[i] == '\\')
            converted[i] = '/';
        else if (converted[i] == ':')
            offsetStart = (i + 2);

        // Remove leading slashes
        //
        if ((converted[i] == '/') && (i == offsetStart))
            ++offsetStart;
    }

    // Append trailing slash if not present.
    //
    if ((numChars > 0) && (converted[numChars - 1] != '/')) {
        converted += "/";
        ++numChars;
    }

    // Remove drive letter, colon and slash if present in windows OS
    //
    if ((offsetStart > 0) && (offsetStart < numChars))
        converted = converted.substr(offsetStart, numChars);

    // Convert to lower case
    //
    std::transform(converted.begin(), converted.end(), converted.begin(), tolower);
    return converted;
}

bool ZipArchive::checkFileHandleValid() const {
    if (((archiveExists_ == true) && (archive_ == 0)) || (archive_->isOpen() == false)) {
        std::cout << "\teof()?: " << archive_->eof() << "\n";
        std::cout << "\topen()?: " << archive_->isOpen() << "\n";
        LERROR("Archive is not opened for reading! Call open() first.");
        return false;
    }
    return true;
}

size_t ZipArchive::copyExistingFiles(std::ofstream& ofs,
                                     std::vector<ZipArchive::ArchivedFile>& existingFiles)
{
    size_t counter = 0;
    for (size_t i = 0; i < existingFiles.size(); ++i) {
        if (existingFiles[i].isNewInArchive_ == true) {
            LWARNING("copyExistingFiles(): attempting to copy non-existing file '");
            LWARNING(existingFiles[i].fileName_ << "'!");
            continue;
        }

        std::streampos outFileOffset = ofs.tellp();

        bool error = false;
        ZipLocalFileHeader& lfh = existingFiles[i].zipLocalFileHeader_;
        if (lfh.signature != SIGNATURE_ZIPLOCALFILEHEADER)
            error = ! readLocalFileHeader(existingFiles[i].zipLocalFileHeader_,
                existingFiles[i].localHeaderOffset_);
        else
            archive_->seek(existingFiles[i].localHeaderOffset_ + SIZE_ZIPLOCALFILEHEADER, File::BEGIN);

        size_t read = SIZE_ZIPLOCALFILEHEADER;
        if (error == false) {
            ofs.write(reinterpret_cast<char*>(&lfh), SIZE_ZIPLOCALFILEHEADER);
            error = ofs.fail();
        }

        // Read and write file name...
        //
        if ((error == false) && (lfh.filenameLength > 0)) {
            char* buffer = new char[lfh.filenameLength];
            read += archive_->read(buffer, lfh.filenameLength);
            ofs.write(buffer, lfh.filenameLength);
            error = ofs.fail();
            delete [] buffer;
        }

        // Read and write extra file field...
        //
        if ((error == false) && (lfh.extraFieldLength > 0)) {
            char* buffer = new char[lfh.extraFieldLength];
            read += archive_->read(buffer, lfh.extraFieldLength);
            ofs.write(buffer, lfh.extraFieldLength);
            error = ofs.fail();
            delete [] buffer;
        }

        if (error == false) {
            size_t bufferSize = (lfh.compressedSize < MAX_BUFFER_SIZE) ? lfh.compressedSize : MAX_BUFFER_SIZE;
            char* buffer = new char[bufferSize];
            size_t readTotal = 0;
            do {
                size_t r = archive_->read(buffer, bufferSize);
                readTotal += r;
                ofs.write(buffer, r);
                error = ofs.fail();
                if (error == true)
                    break;
            } while ((readTotal < lfh.compressedSize) && (archive_->eof() == false));
            delete [] buffer;
            read += lfh.compressedSize;
        }

        if (error == true) {
            LERROR("copyExistingFiles(): failed to copy data!");
            break;
        }

        // Adjust the offset of the LocalFileHeader structure for the new output file...
        //
        ZipFileHeader& fileHeader = existingFiles[i].zipFileHader_;
        fileHeader.localHeaderOffset = (outFileOffset >= 0) ? static_cast<uint32_t>(outFileOffset) : 0xFFFFFFFF;
        std::cout << "\tlocalHeaderOffset = " << outFileOffset << "\n";
        ++counter;
    }   // for (i
    return counter;
}

bool ZipArchive::prepareDirectories(const std::string& zipFileName) const {
    for (size_t lastOccurance = 0, strOffset = 0; strOffset != std::string::npos; )
    {
        strOffset = zipFileName.find('/', lastOccurance);
        if ((strOffset != std::string::npos) && (lastOccurance < strOffset)) {
            std::string dirName = zipFileName.substr(0, strOffset);
            if (FileSystem::dirExists(dirName) == false) {
                if (FileSystem::createDirectory(dirName) == false) {
                    LERROR("prepareDirectories(): failed to create directory '" << dirName << "'!");
                    return false;
                }
            }
            lastOccurance = strOffset + 1;
        }
    }
    return true;
}

bool ZipArchive::readEOCDHeaderRecord(ZipArchive::ZipEOCDHeaderRecord& eocdHeaderRec)
{
    if (checkFileHandleValid() == false)
        return false;

    // If no archive comment is present, the End Of Central Directory Record
    // is the last entry within the file.
    //
    archive_->seek(-SIZE_ZIPEOCDRECORD, File::END);
    archive_->read(&eocdHeaderRec, SIZE_ZIPEOCDRECORD);
    if (eocdHeaderRec.signature == SIGNATURE_ZIPEOCDHEADERRECORD)
        return true;

    memset(&eocdHeaderRec, 0, SIZE_ZIPEOCDRECORD);

    // The comment may be of max. 64 KByte length, so search for a the EOCDHeaderRecord
    // in a block of 64 KB + sizeof(EOCDHeaderRecord) starting from the end of the file
    //
    size_t fileSize = archive_->size();
    size_t maxBufferSize = static_cast<size_t>(0xFFFF) + SIZE_ZIPEOCDRECORD;
    long bufferSize = static_cast<long>((maxBufferSize < fileSize) ? maxBufferSize : fileSize);
    char* buffer = new char[bufferSize];
    archive_->seek(-bufferSize, File::END);
    archive_->read(buffer, bufferSize);


    // Wrap a std::string around the buffer in order to facilitate searching for the
    // signature.
    //
    std::string content(buffer, bufferSize);
    char signature[4] = {'P', 'K', 0x05, 0x06};
    size_t pos = content.rfind(signature);
    if (pos == std::string::npos) {
        delete [] buffer;
        return false;
    }

    // If the signature could be found, copy the EOCDHeaderRecord
    //
    memcpy(&eocdHeaderRec, (buffer + pos), SIZE_ZIPEOCDRECORD);
    delete [] buffer;
    return true;
}

bool ZipArchive::readFileHeader(ZipArchive::ZipFileHeader& fileHeader, size_t fileOffset) {
    if (checkFileHandleValid() == false)
        return false;
    archive_->seek(fileOffset, File::BEGIN);
    archive_->read(&fileHeader, SIZE_ZIPFILEHEADER);
    if (fileHeader.signature == SIGNATURE_ZIPFILEHEADER)
        return true;

    memset(&fileHeader, 0, SIZE_ZIPFILEHEADER);
    return false;
}

bool ZipArchive::readLocalFileHeader(ZipArchive::ZipLocalFileHeader& localFileHeader,
                                         size_t fileOffset)
{
    if (checkFileHandleValid() == false)
        return false;
    archive_->seek(fileOffset, File::BEGIN);
    archive_->read(&localFileHeader, SIZE_ZIPLOCALFILEHEADER);
    if (localFileHeader.signature == SIGNATURE_ZIPLOCALFILEHEADER)
        return true;

    memset(&localFileHeader, 0, SIZE_ZIPLOCALFILEHEADER);
    return false;
}

std::string ZipArchive::readString(size_t fileOffset, size_t numChars) {
    if ((checkFileHandleValid() == false) || (numChars == 0))
        return "";

    char* buffer = new char[numChars + 1];
    archive_->seek(fileOffset, File::BEGIN);
    archive_->read(buffer, numChars);
    buffer[numChars] = '\0';

    std::string str(buffer, numChars);
    delete [] buffer;
    return str;
}

bool ZipArchive::readZipFile() {
    if (checkFileHandleValid() == false)
        return false;

    ZipEOCDHeaderRecord eocdHeaderRec;
    if (readEOCDHeaderRecord(eocdHeaderRec) == false)
        return false;

    if (eocdHeaderRec.numberOfDisk != eocdHeaderRec.numberOfDiskWithStartOfCD) {
        LERROR("Need disk #" << eocdHeaderRec.numberOfDiskWithStartOfCD << " to read central directory!");
        return false;
    }

    if (eocdHeaderRec.numberOfEntriesInCD <= 0) {
        LERROR("Error: central directory does not contain any entries!");
        return false;
    }

    size_t offset = eocdHeaderRec.offsetStartCD;
    for (uint16_t i = 0; i < eocdHeaderRec.numberOfEntriesInCD; ++i) {
        ArchivedFile af;
        ZipFileHeader& fileHeader = af.zipFileHader_;
        bool res = readFileHeader(af.zipFileHader_, offset);
        if (res == true) {
            if (fileHeader.generalPurposeFlag != 0x0000) {
                LWARNING("A file seems to make use of advanced features ");
                LWARNING("this reader is unable to understand. skipping...");
                continue;
            }

            bool validCompression = false;
            switch (fileHeader.compressionMethod) {
                case 0:
                case 8:
                    validCompression = true;
                    break;
                default:
                    break;
            }

            if (validCompression == false) {
                LWARNING("Unsupported compression method (code = " << fileHeader.compressionMethod << ")!");
                LWARNING("This archiver only supports standard deflate (code 8) and uncompressed (code 0).");
                continue;
            }

            size_t fileNameOffset = (offset + SIZE_ZIPFILEHEADER);
            size_t extraFieldOffset = fileNameOffset + fileHeader.filenameLength;
            size_t commentOffset = extraFieldOffset + fileHeader.extraFieldLength;

            if (fileHeader.filenameLength > 0)
                af.fileName_ = readString(fileNameOffset, fileHeader.filenameLength);

            if (fileHeader.extraFieldLength > 0)
                af.fileExtra_ = readString(extraFieldOffset, fileHeader.extraFieldLength);

            if (fileHeader.fileCommentLength > 0)
                af.fileComment_ = readString(commentOffset, fileHeader.fileCommentLength);

            af.isNewInArchive_ = false;
            af.localHeaderOffset_ = fileHeader.localHeaderOffset;
            std::pair<ArchiveMap::iterator, bool> rs =
                files_.insert(std::make_pair(af.fileName_, af));
            if (rs.second == false)
                LERROR("failed to extract directory entry for file '" << af.fileName_ << "'!");
        } else
            LERROR("Failed to read directory entry for entry #" << i << "!");

        // The next file header is found at the current offset plus the size of
        // this FileHeader struct (constant) and plus the 3 variable length fields
        //
        offset += (SIZE_ZIPFILEHEADER + fileHeader.filenameLength
            + fileHeader.extraFieldLength + fileHeader.fileCommentLength);
    }

    return true;
}

bool ZipArchive::writeCentralDirectory(std::ofstream& ofs,
                                       const std::vector<ZipArchive::ArchivedFile>& files)
{
    std::streampos offsetCD = ofs.tellp();
    size_t sizeCD = 0;
    bool error = false;
    for (size_t i = 0; i < files.size(); ++i) {
        const ZipFileHeader& fileHeader = files[i].zipFileHader_;
        ofs.write(reinterpret_cast<const char*>(&fileHeader), SIZE_ZIPFILEHEADER);
        sizeCD += SIZE_ZIPFILEHEADER;
        error = ofs.fail();

        size_t stringSize = files[i].fileName_.size();
        sizeCD += stringSize;
        if ((error == false) && (stringSize > 0) && (stringSize == fileHeader.filenameLength)) {
            ofs.write(files[i].fileName_.c_str(), stringSize);
            error = ofs.fail();
        }

        stringSize = files[i].fileExtra_.size();
        sizeCD += stringSize;
        if ((error == false) && (stringSize > 0) && (stringSize == fileHeader.extraFieldLength)) {
            ofs.write(files[i].fileExtra_.c_str(), stringSize);
            error = ofs.fail();
        }

        stringSize = files[i].fileComment_.size();
        sizeCD += stringSize;
        if ((error == false) && (stringSize > 0) && (stringSize == fileHeader.fileCommentLength)) {
            ofs.write(files[i].fileComment_.c_str(), stringSize);
            error = ofs.fail();
        }

        if (error == true) {
            LERROR("writeCentralDirectory(): failed to write Central Directory Entry in zip archive!");
            break;
        }
    }   // for

    ZipEOCDHeaderRecord eocd;
    eocd.signature = SIGNATURE_ZIPEOCDHEADERRECORD;
    eocd.numberOfDisk = 0;
    eocd.numberOfDiskWithStartOfCD = 0;
    eocd.numberOfEntriesInThisCD = static_cast<uint16_t>(files.size());
    eocd.numberOfEntriesInCD = static_cast<uint16_t>(files.size());
    eocd.sizeOfCD = static_cast<uint32_t>(sizeCD);
    eocd.offsetStartCD = (offsetCD >= 0) ? static_cast<uint32_t>(offsetCD) : 0;
    eocd.commmentLength = 0;

    if (error == false) {
        ofs.write(reinterpret_cast<char*>(&eocd), SIZE_ZIPEOCDRECORD);
        error = ofs.fail();
    }

    return (! error);
}

size_t ZipArchive::writeNewFiles(std::ofstream& ofs,
                                 std::vector<ZipArchive::ArchivedFile>& newFiles)
{
    size_t counter = 0;
    bool error = false;
    for (size_t i = 0; i < newFiles.size(); ++i) {
        if (newFiles[i].isNewInArchive_ == false) {
            LWARNING("writeNewFiles(): attempting to add an existing or non-physical file '");
            LWARNING(newFiles[i].fileName_ << "' as new file!\n");
            continue;
        }

        if ((newFiles[i].extFileName_.empty() == true) && (newFiles[i].extHandle_ == 0)) {
            LERROR("wirteNewFiles(): the external file name for the new file is empty!");
            continue;
        }

        if (newFiles[i].fileName_.size() > 0xFFFF) {
            LERROR("writeNewFiles(): the file name may not exceed (2^16 - 1) bytes!");
            continue;
        }

        if ((FileSystem::fileExists(newFiles[i].extFileName_) == false) && (newFiles[i].extHandle_ == 0)) {
            LERROR("file '" << newFiles[i].extFileName_ << "' does not exist!");
            continue;
        } else
            LINFO("Writing file " << newFiles[i].fileName_);

        // Seek forward from the current position to skip the LocalFileHeader struct and the file
        // name. The space has to be kept free until we know all details which are available
        // when deflateToDisk() has finished.
        //
        std::streampos archiveOffset = ofs.tellp();
        ofs.seekp(SIZE_ZIPLOCALFILEHEADER + newFiles[i].fileName_.size(), std::ios_base::cur);

        unsigned long crc32 = 0;
        size_t compressedSize = 0;
        size_t uncompressedSize = 0;
        if ((newFiles[i].extFileName_.empty() == false) && (newFiles[i].extHandle_ == 0)) {
            RegularFile inFile(newFiles[i].extFileName_);
            compressedSize = deflateToDisk(inFile, ofs, crc32);
            uncompressedSize = inFile.size();
        } else if (newFiles[i].extHandle_ != 0) {
            compressedSize = deflateToDisk(*(newFiles[i].extHandle_), ofs, crc32);
            uncompressedSize = newFiles[i].extHandle_->size();
        }

        error = (compressedSize == 0);

        // The compressed data now have been written and the header + file name
        // can be inserted.
        //
        if (error == false) {
            ZipLocalFileHeader& lfh = newFiles[i].zipLocalFileHeader_;
            lfh.signature = 0x04034b50;
            lfh.versionNeeded = ZIP_VERSION;
            lfh.generalPurposeFlag = 0;
            lfh.compressionMethod = 0x0008;
            lfh.lastModTime = 0;    // Argh! Need last modification time of file in DOS (!) style...
            lfh.lastModDate = 0;    // Argh! Need last modification date of file in DOS (!) style...
            lfh.crc32 = crc32;
            lfh.compressedSize = static_cast<uint32_t>(compressedSize);
            lfh.uncompressedSize = static_cast<uint32_t>(uncompressedSize);
            lfh.filenameLength = static_cast<uint16_t>(newFiles[i].fileName_.size());
            lfh.extraFieldLength = 0;

            // Seek the position where the LocalFileHeader has to be placed
            //
            ofs.seekp(archiveOffset, std::ios_base::beg);
            ofs.write(reinterpret_cast<char*>(&lfh), SIZE_ZIPLOCALFILEHEADER);
            error = ofs.fail();
        }

        if (error == false) {
            ofs.write(newFiles[i].fileName_.c_str(), newFiles[i].fileName_.size());
            error = ofs.fail();
        }

        if (error == true) {
            LERROR("writeNewFiles(): an error has occured while writing to archive!");
            break;
        } else {
            // The stream position is now again at the beginning of the block containing
            // the compressed data, so seek forward to find the position where the
            // next file can be written.
            //
            ofs.seekp(compressedSize, std::ios_base::cur);

            // Adjust the state of the file within this archive: it is not new any longer
            // and the file name has been adjusted. Furthermore the file now has a
            // offset for the LocalFileHeader within this new archive...
            //
            newFiles[i].isNewInArchive_ = false;
            newFiles[i].localHeaderOffset_ = static_cast<uint32_t>(archiveOffset);
            newFiles[i].fileExtra_ = "";
            newFiles[i].fileComment_ = "";

            // In addition, the entry for the Central Directory can be created.
            //
            ZipFileHeader& fileHeader = newFiles[i].zipFileHader_;
            fileHeader.signature = SIGNATURE_ZIPFILEHEADER;
            fileHeader.versionMadeBy = ZIP_VERSION;
            fileHeader.versionNeeded = ZIP_VERSION;
            fileHeader.generalPurposeFlag = 0;
            fileHeader.compressionMethod = 0x0008;
            fileHeader.lastModTime = 0;
            fileHeader.lastModDate = 0;
            fileHeader.crc32 = crc32;
            fileHeader.compressedSize = static_cast<uint32_t>(compressedSize);
            fileHeader.uncompressedSize = static_cast<uint32_t>(uncompressedSize);
            fileHeader.filenameLength = static_cast<uint16_t>(newFiles[i].fileName_.size());
            fileHeader.extraFieldLength = static_cast<uint16_t>(newFiles[i].fileExtra_.size());
            fileHeader.fileCommentLength = static_cast<uint16_t>(newFiles[i].fileComment_.size());
            fileHeader.diskNumberStart = 0;
            fileHeader.internalFileAttributes = 0;
            fileHeader.externalFileAttributes = 0;
            fileHeader.localHeaderOffset = static_cast<uint32_t>(newFiles[i].localHeaderOffset_);
        }
        ++counter;
    }   // for (i
    return counter;
}

}   // namespace tgt
