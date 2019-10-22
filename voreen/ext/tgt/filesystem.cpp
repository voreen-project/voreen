/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/filesystem.h"

#include "tgt/types.h"

#include <algorithm>
#include <string>
#include <cstring>
#include <cctype> // std::tolower
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <sys/stat.h>
#include <time.h>

#ifndef WIN32
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <utime.h>
#include <dirent.h>
#else
#include <windows.h>
#include <tchar.h>
#include <direct.h>
#include <sys/utime.h>
#endif

using std::string;

//TODO: should really switch to exceptions instead of just logging with LERROR. joerg
// see http://en.wikipedia.org/wiki/Aliasing_%28computing%29#Conflicts_with_optimization

// anonymous namespace
namespace {

enum {
    TEMP_BUFFERSIZE = 128,
    LINE_BUFFERSIZE = 512,

    /*
     * Standard Archive Format - Standard TAR - USTAR
     */
    RECORDSIZE = 512,
    NAMSIZ = 100,
    TUNMLEN = 32,
    TGNMLEN = 32
};

//TODO: these definitions with arrays of chars break strict-aliasing rules for gcc

union record {
    char charptr[RECORDSIZE];
    struct header {
        char    name[NAMSIZ];
        char    mode[8];
        char    uid[8];
        char    gid[8];
        char    size[12];
        char    mtime[12];
        char    chksum[8];
        char    linkflag;
        char    linkname[NAMSIZ];
        char    magic[8];
        char    uname[TUNMLEN];
        char    gname[TGNMLEN];
        char    devmajor[8];
        char    devminor[8];
    } header;
};

size_t findDelim(char* buf, size_t size, char delim) {
    size_t pos = 0;
    while (pos < size) {
        if (buf[pos] == delim)
            return pos;
        ++pos;
    }
    return size+1;
}

} // anonymous namespace

namespace tgt {

File::File()
    : size_(0)
{}

File::File(const std::string& name)
    : name_(name),
      size_(0)
{}

File::~File()
{}

size_t File::readLine(char *buf, size_t maxCount, char delim) {
    size_t bytesread = 0;
    size_t start = tell();
    char* p = buf;

    if (tell()+maxCount > size())
        maxCount = size()-tell();
    while (bytesread < maxCount) {
        size_t bufsize = TEMP_BUFFERSIZE;
        if ((bytesread + TEMP_BUFFERSIZE) > maxCount)
            bufsize = (maxCount - bytesread);
        read(p, bufsize);
        size_t d = findDelim(p, bufsize, delim);
        if (d < bufsize) {
            bytesread += d;
            p[d] = 0;
            seek(start + bytesread + 1);
            return bytesread;
        }
        p += bufsize;
        bytesread += bufsize;
    }
    return bytesread;
}

std::string File::getLine(char delim) {
    char* buf = new char[LINE_BUFFERSIZE];
    size_t br = readLine(buf, LINE_BUFFERSIZE, delim);
    std::string ret(buf, br);
    while (br == LINE_BUFFERSIZE) {
        br = readLine(buf, LINE_BUFFERSIZE, delim);
        ret.append(buf, br);
    }
    delete[] buf;
    return ret;
}

std::string File::getAsString(void) {
    char* buf = new char[size()];
    read(buf, size());
    std::string ret(buf, size());
    delete[] buf;
    return ret;
}

size_t File::skipLine(char delim) {
    size_t bytesread = 0;
    size_t start = tell();
    char* p = new char[TEMP_BUFFERSIZE];
    while (true) {
        read(p, TEMP_BUFFERSIZE);
        size_t d = findDelim(p, TEMP_BUFFERSIZE, delim);
        if (d < TEMP_BUFFERSIZE) {
            bytesread += d-1;
            delete[] p;
            seek(start + bytesread + 2);
            return bytesread;
        }
        bytesread += TEMP_BUFFERSIZE;
    }
}

//-----------------------------------------------------------------------------

const std::string RegularFile::loggerCat_("tgt.RegularFile");

RegularFile::RegularFile(const std::string& filename)
    : File(filename)
{
    file_.open(filename.c_str(), std::ios::binary);

    if (!file_) {
        LDEBUG("Cannot open file: " << filename);
        return;
    }

    // get file size
    file_.seekg(0, std::ios::end);
    size_ = tell();
    file_.seekg(std::ios::beg);
}

RegularFile::~RegularFile() {
    close();
}

void RegularFile::close() {
    if (isOpen())
        file_.close();
}

size_t RegularFile::read(void* buf, size_t count) {
    if (!isOpen() || eof())
        return 0;

    file_.read(static_cast<char*>(buf), count);
    return static_cast<size_t>(file_.gcount());
}

void RegularFile::skip(size_t count) {
    if (isOpen())
        file_.seekg(count, std::ios::cur);
}

void RegularFile::seek(std::streamoff pos) {
    if (isOpen())
        file_.seekg(pos);
}

void RegularFile::seek(std::streamoff offset, File::SeekDir seekDir) {
    if (isOpen()) {
        switch (seekDir)  {
        case File::BEGIN:
            file_.seekg(offset, std::ios_base::beg);
            break;
        case File::CURRENT:
            file_.seekg(offset, std::ios_base::cur);
            break;
        case File::END:
            file_.seekg(offset, std::ios_base::end);
            break;
        }
    }
}

size_t RegularFile::tell() {
    if (isOpen()) {
        std::streampos p = file_.tellg();
        if (p >= 0)
            return static_cast<size_t>(p);
    }
    return 0;
}

bool RegularFile::eof() {
    if (isOpen()) {
        if (tell() == size())
            return true;
        else
            return file_.eof();
    } else {
        return true;
    }
}

bool RegularFile::isOpen() {
    return file_.is_open();
}

bool RegularFile::good() {
    if (isOpen())
        return file_.good();
    else
        return false;
}

//-----------------------------------------------------------------------------

MemoryFile::MemoryFile(const char* data, size_t size, const std::string& filename, bool deleteData)
  : File(filename),
  data_(data),
  pos_(0),
  deleteData_(deleteData)
{
    size_ = size;
}

MemoryFile::~MemoryFile(){
    if (data_ && deleteData_)
        delete[] data_;
}

void MemoryFile::close() {
    pos_ = 0;
}

size_t MemoryFile::read(void* buf, size_t count) {
    if ((pos_ + count) < size_) {
        memcpy(buf, (const void*) (data_ + pos_), count);
        pos_ += count;
        return count;
    }
    else {
        if (eof())
            return 0;
        size_t r = size_ - pos_;
        memcpy(buf, (const void*) (data_ + pos_), r);
        pos_ = size_ - 1;
        return r;
    }
}

void MemoryFile::skip(size_t count) {
    pos_ += count;
    if (pos_ >= size_)
        pos_ = size_ - 1;
}

void MemoryFile::seek(std::streamoff pos) {
    pos_ = static_cast<size_t>(pos);
    if (pos_ >= size_)
        pos_ = size_ - 1;
}

void MemoryFile::seek(std::streamoff offset, File::SeekDir seekDir) {
    switch (seekDir) {
        case File::BEGIN:
            seek(offset);
            break;
        case File::CURRENT:
            seek(pos_+offset);
            break;
        case File::END:
            pos_ = size_ - 1;
            break;
    }
}

size_t MemoryFile::tell() {
    return pos_;
}

bool MemoryFile::eof() {
    return (pos_ >= size_ - 1);
}

bool MemoryFile::isOpen() {
    return true;
}

bool MemoryFile::good() {
    // only possibility not to be good is on EOF
    return !eof();
}

//-----------------------------------------------------------------------------

const std::string TarFile::loggerCat_("tgt.TarFile");

TarFile::TarFile(const std::string& filename, const std::string& tarfilename, size_t offset, size_t size)
    : File(filename)
{
    file_ = FileSys.open(tarfilename.c_str());
    name_ = filename;

    // Check if file is open:
    if (!file_) {
        std::cout << "Failed to open tar file: " << tarfilename << std::endl;
        return;
    }

    if (!file_->good()) {
        std::cout << "Bad tar file: " << tarfilename << std::endl;
        file_->close();
        delete file_;
        file_ = 0;
        return;
    }

#ifdef TGT_DEBUG
    file_->seek(0, File::END);
    size_ = file_->tell();
    file_->seek(File::BEGIN);

    if ((offset + size) > size_) {
        std::cout << "Wrong tar file size/offset! file: " << filename << " in " << tarfilename << std::endl;
        std::cout << "size_ " << size_ << " offset " << offset << " size " << size << std::endl;
        file_->close();
        delete file_;
        file_ = 0;
        return;
    }
#endif // TGT_DEBUG

    size_ = size;
    offset_ = offset;
    file_->seek(offset_, File::BEGIN);
}

TarFile::~TarFile() {
    close();
}

void TarFile::close() {
    if (file_) {
        file_->close();
        delete file_;
        file_ = 0;
    }
}

size_t TarFile::read(void* buf, size_t count) {
    if (tell() + count > size_)
        count = size_ - tell();
    if (!file_)
        return 0;
    file_->read(static_cast<char*>(buf), count);
    return count;
}

void TarFile::skip(size_t count) {
    if (tell() + count > size_)
        count = size_ - tell();

    if (file_)
        file_->seek(count, File::CURRENT);
}

void TarFile::seek(std::streamoff pos) {
    if (file_)
        file_->seek(offset_ + pos);
}

void TarFile::seek(std::streamoff offset, File::SeekDir seekDir) {
    if (file_)
        file_->seek(offset, seekDir);
}

size_t TarFile::tell() {
    if (file_)
        return ((size_t)(file_->tell())) - offset_;
    else
        return 0;
}

bool TarFile::eof() {
    if (file_)
       return (tell() == size_);
    else
        return true;
}

bool TarFile::isOpen() {
    return (file_ != 0);
}

bool TarFile::good() {
    if (eof())
        return false;

    if (file_)
        return file_->good();
    else
        return false;
}


//-----------------------------------------------------------------------------

MemoryFileFactory::MemoryFileFactory(const std::string& filename, char* data, size_t size) {
    data_ = data;
    size_ = size;
    filename_ = filename;
}

MemoryFileFactory::MemoryFileFactory(const std::string& filename, const std::string& data) {
    size_ = data.length();
    data_ = new char[size_];
    memcpy(data_, data.c_str(), size_);
    filename_ = filename;
}

File* MemoryFileFactory::open(const std::string& filename) {
    if (filename == filename_)
        return new MemoryFile(data_, size_, filename_);
    else
        return 0;
}

std::vector<std::string> MemoryFileFactory::getFilenames() {
     std::vector<std::string> files;
     files.push_back(filename_);
     return files;
}

//-----------------------------------------------------------------------------

const std::string TarFileFactory::loggerCat_("tgt.TarFileFactory");

TarFileFactory::TarFileFactory(const std::string& filename, const std::string& rootpath) {
    filename_ = filename;
    File* file = FileSys.open(filename);

    // Check if file is open
    if (!file) {
        LERROR("Failed to open TAR archive " << filename);
        return;
    }

    if (!file->good()) {
        LERROR("Error opening TAR archive " << filename);
        return;
    }

    bool end = false;

    record* currecord = new record;

    while (!end) {
        file->read(currecord->charptr, RECORDSIZE);

        bool allzero = true;
        //check if record is all zero (=>eof)
        for (int i = 0; i < 512; i++) {
            if (currecord->charptr[i] != 0)
                allzero = false;
        }

        if (!allzero) {
            //convert from ascii-octal to normal int
            size_t filesize = 0;
            size_t mul = 1;
            for (int i = 10; i >= 0; i--) {
                filesize += ((currecord->header.size[i] - 48) * mul);
                mul *= 8;
            }

            if (filesize > 0) {
                size_t blocks = (filesize / RECORDSIZE);
                if ((filesize % RECORDSIZE) != 0)
                    ++blocks;
                LDEBUG("Found file " << currecord->header.name);
                LDEBUG("  Add it as " << rootpath+currecord->header.name
                       << "; filesize: " << filesize << "; offset " << file->tell());
                files_[rootpath+currecord->header.name].size_ = filesize;
                files_[rootpath+currecord->header.name].offset_ = file->tell();

                file->seek(RECORDSIZE*blocks, File::CURRENT);
            }

            if (file->eof())
                end = true;
        }
        else
            end = true;
    }
    file->close();
    delete file;
}

File* TarFileFactory::open(const std::string& filename) {
    LDEBUG("Loading file " << filename << " from " << filename_);
    LDEBUG("  Offset: " << files_[filename].offset_ << "; size: " << files_[filename].size_);
    TarFile* tf = new TarFile(filename, filename_, files_[filename].offset_,
                              files_[filename].size_);
    return tf;
}

std::vector<std::string> TarFileFactory::getFilenames() {
    std::vector<std::string> files;

    std::map<std::string, ArchivedFile>::iterator theIterator;
    for (theIterator = files_.begin(); theIterator != files_.end(); theIterator++)
        files.push_back((*theIterator).first);

    return files;
}


//-----------------------------------------------------------------------------

const std::string FileSystem::loggerCat_("tgt.FileSystem.FileSystem");
#ifdef WIN32
const char FileSystem::goodSlash_ = '\\';
const char FileSystem::badSlash_ = '/';
#else
const char FileSystem::goodSlash_ = '/';
const char FileSystem::badSlash_ = '\\';
#endif

FileSystem::FileSystem() {
}

FileSystem::~FileSystem() {
    while (!factories_.empty()) {
        delete (factories_.back());
        factories_.pop_back();
    }
}

File* FileSystem::open(const std::string& filename) {
    if (virtualFS_.find(filename) != virtualFS_.end()) {
        LDEBUG("Opening file " << filename << " from virtualFS");
        return virtualFS_[filename]->open(filename);
    } else {
        LDEBUG("Opening file " << filename << " from real FS");

        RegularFile* f = new RegularFile(filename);
        if (f->isOpen()) {
            return f;
        } else {
            delete f;
            return 0;
        }
    }
}

bool FileSystem::exists(const std::string& filename) {
    if (virtualFS_.find(filename) != virtualFS_.end()) {
        LDEBUG("Not checking if file " << filename << " exists in virtualFS, unimplemented.");
        //FIXME: this does not check if file exists in virtualFS, unimplemented! joerg
        return true;
    } else {
        std::ifstream file(filename.c_str(), std::ios::binary);
        return file.good();
    }
}

void FileSystem::addFactory(FileFactory* ff) {
    std::vector<std::string> files = ff->getFilenames();
    std::vector<std::string>::iterator theIterator;
    for (theIterator = files.begin(); theIterator != files.end(); theIterator++)
        virtualFS_[*theIterator] = ff;

    factories_.push_back(ff);
}

void FileSystem::addMemoryFile(const std::string& filename, char* data, size_t size) {
    FileFactory* ff = new MemoryFileFactory(filename, data, size);
    addFactory(ff);
}

void FileSystem::addMemoryFile(const std::string& filename, const std::string& data) {
    FileFactory* ff = new MemoryFileFactory(filename, data);
    addFactory(ff);
}

void FileSystem::addPackage(const std::string& filename, const std::string& rootpath) {
    LINFO("adding package " << filename);
    LDEBUG("root path: " << rootpath);

    //TODO: do something like TextureReader in TexMgr...
    std::string::size_type loc = filename.find(".tar", 0);
    if (loc != std::string::npos) {
        LDEBUG("Recognized tar file");
        addFactory(new TarFileFactory(filename, rootpath));
    }
}

namespace {

const string PATH_SEPARATORS = "/\\";

std::string replaceAllCharacters(const std::string& name, const char oldChar,
                                 const char newChar)
{
    // the simplest algorithm but it should be fast enough
    //
    std::string conv(name);
    for (size_t i = 0; i < conv.size(); ++i) {
        if (conv[i] == oldChar)
            conv[i] = newChar;
    }
    return conv;
}

size_t removeTrailingCharacters(std::string& str, const char trailer) {
    if (str.empty())
        return 0;

    size_t pos = str.find_last_not_of(trailer);
    if (pos != std::string::npos) {
        size_t count = str.size() - (pos + 1);
        str.resize(pos + 1);
        return count;
    }

    return 0;
}

} // namespace


//
// static methods for the regular filesystem
//

bool FileSystem::isAbsolutePath(const std::string& path) {
    return ((path.size() > 0 && path.at(0) == '/') || (path.size() > 1 && path.at(1) ==':'));
}

static bool tryMakeAbsolute(const string& path, string& output) {
    char* buffer;
#ifdef WIN32
    buffer = static_cast<char*>(malloc(4096));
    buffer[0] = 0;
    if (GetFullPathName(path.c_str(), 4096, buffer, 0) == 0) {
        free(buffer);
        buffer = 0;
    }
#elif _POSIX_VERSION >= 200809L || defined (linux)
    // use safe realpath if available
    buffer = realpath(path.c_str(), 0);
#else
    char* resolvedPath = static_cast<char*>(malloc(4096));
    buffer = realpath(path.c_str(), resolvedPath);
    // on success buffer is equal to resolvedPath and gets freed later
    // on failure we have to free resolvedPath
    if (!buffer)
        free(resolvedPath);
#endif

    if (buffer) {
        string result(buffer);
        free(buffer);
        output = tgt::FileSystem::cleanupPath(result);
        return true;
    } else {
        output = path;
        return false;
    }
}

string FileSystem::absolutePath(const string& path) {
    string absprefix;
    string prefix = path;
    string suffix = "";
    while(!tryMakeAbsolute(prefix, absprefix)) {
        if(prefix.empty() || prefix == "/") {
          return path;
        }
        suffix = fileName(prefix) + "/" + suffix;
        prefix = dirName(prefix);
    }
    return cleanupPath(absprefix + "/" + suffix);
}

std::string FileSystem::relativePath(const std::string& path, const std::string& dir) {
    // when there is no dir or path we just return the path
    if (dir.empty() || path.empty())
        return path;

    // make paths absolute and add trailing separator
    string abspath = cleanupPath(absolutePath(path)) + "/";
    string absdir = cleanupPath(absolutePath(dir)) + "/";

    // if absolute path could not be obtained, just return the input path
    if(!isAbsolutePath(abspath) || !isAbsolutePath(absdir))
        return path;

    // if both directories are the same we return an empty relative path
    if (abspath.compare(absdir) == 0)
        return "";

    /*
    // catch differing DOS-style drive names
    if (abspath.size() < 1 || abspath.size() < 1 || abspath[0] != absdir[0]) {
        std::transform(abspath.begin(), abspath.begin()+1, abspath.begin(), static_cast<int (*)(int)>(std::tolower));
        std::transform(absdir.begin(), absdir.begin()+1, absdir.begin(), static_cast<int (*)(int)>(std::tolower));
    }*/
#ifdef WIN32
    // do not generate a relative patch across different Windows drives
    if (!abspath.empty() && !absdir.empty() && (abspath[0] != absdir[0])) {
        return path;
    }
#endif

    // find common part in path and dir string
    string::size_type pospath = abspath.find_first_of(PATH_SEPARATORS);
    string::size_type posdir = absdir.find_first_of(PATH_SEPARATORS);
    size_t i = 0;
    while (abspath.compare(0, pospath, absdir, 0, posdir) == 0) {
        i = pospath;
        pospath = abspath.find_first_of(PATH_SEPARATORS, pospath + 1);
        posdir = absdir.find_first_of(PATH_SEPARATORS, posdir + 1);
    }

    // now we have remaining then non-common parts of both paths
    string restpath = abspath.substr(i + 1);
    string restdir = absdir.substr(i + 1);

    // the remaining path is our initial relative path
    string relative = restpath;

    // add ".." for each path separator in the remaining part of dir
    string::size_type pos = restdir.find_first_of(PATH_SEPARATORS);
    while (pos != string::npos) {
        relative = "../" + relative;
        pos = restdir.find_first_of(PATH_SEPARATORS, pos + 1);
    }

    // cleanup and return result
    return cleanupPath(relative);
}

string FileSystem::fileName(const string& filepath) {
    string::size_type separator = filepath.find_last_of("/\\");
    if (separator != string::npos)
        return filepath.substr(separator + 1);
    else
        return filepath;
}

string FileSystem::baseName(const string& filepath) {
    string filename = fileName(filepath);
    string::size_type dot = filename.rfind(".");

    if (dot != string::npos)
        return filename.substr(0, dot);
    else
        return filename;
}

string FileSystem::fullBaseName(const string& filepath) {
    string::size_type dot = filepath.rfind(".");

    if (dot != string::npos)
        return filepath.substr(0, dot);
    else
        return filepath;
}

string FileSystem::dirName(const std::string& filepath) {
    if (dirExists(filepath)) //< directory passed
        return filepath;

    string::size_type separator = filepath.find_last_of("/\\");
    if (separator != string::npos)
        return filepath.substr(0, separator);
    else
        return "";
}

std::string FileSystem::parentDir(const std::string& dir) {
    std::string curDir = cleanupPath(dir);

    string::size_type separator = curDir.find_last_of("/\\");
    if (separator != string::npos)
        return curDir.substr(0, separator);
    else
        return "";
}

string FileSystem::fileExtension(const string& path, bool lowercase) {
    string filename = fileName(path);

    string::size_type dot = filename.rfind(".");
    string extension;
    if (dot != string::npos)
        extension = filename.substr(dot + 1);

    if (lowercase)
        std::transform(extension.begin(), extension.end(), extension.begin(), tolower);

    return extension;
}

bool FileSystem::comparePaths(const std::string& path1, const std::string& path2) {

    std::string pathAbs1 = cleanupPath(absolutePath(path1));
    std::string pathAbs2 = cleanupPath(absolutePath(path2));

    return (pathAbs1 == pathAbs2);
}

string FileSystem::cleanupPath(std::string path, bool native /*= true*/) {

    if (path.empty())
        return "";

    // remove double path separators
    string::size_type p_start = 0;
#ifdef WIN32
    // do not convert windows network path prefix "\\", if on windows
    if (path.length() >= 2 && path.substr(0, 2) == "\\\\")
        p_start = 2;
#endif
    string::size_type p = p_start;
    while (p != string::npos) {
        string::size_type p_next;
        p_next = path.find("//", p);
        if (p_next == string::npos)
            p_next = path.find("\\\\", p);
        if (p_next == string::npos)
            p_next = path.find("\\/", p);
        if (p_next == string::npos)
            p_next = path.find("/\\", p);

        if (p_next != string::npos)
            path = path.substr(0, p_next) + path.substr(p_next + 1);

        p = p_next;
    }

    // remove './' - note that path must not be empty.
    if(path.length() > 1) {
        string::iterator prev = path.begin(), iter = path.begin();
        while(iter+1 != path.end()) {
            string::iterator next = iter+1;
            if(*iter == '.' && *next == goodSlash_ && *prev != '.') {
                // prev does not change
                iter = path.erase(path.erase(iter));
            }
            else {
                prev = iter;
                iter++;
            }
        }

        // remove leading './' separately
        if(path.length() > 2 && path[0] == '.' && path[1] == goodSlash_)
            path = path.substr(2);
    }

    // remove trailing separator
    if (path.find_last_of(PATH_SEPARATORS) == path.size() - 1)
        path = path.substr(0, path.size() - 1);

    // convert to win or unix separators
#ifdef WIN32
    bool convertToWin = native;
#else
    bool convertToWin = false;
#endif
    if (convertToWin) {
        // convert to native windows separators
        path = replaceAllCharacters(path, '/', '\\');

        // convert drive letter to uppercase
        if ((path.size()>1) && isalpha(path[0]) && (path[1] == ':'))
            std::transform(path.begin(), path.begin()+1, path.begin(), toupper);
    }
    else {
        // convert to unix separators
        tgtAssert(p_start <= path.length(), "invalid p_start param"); // see above
        path = path.substr(0, p_start) + replaceAllCharacters(path.substr(p_start), '\\', '/');
    }

    return path;
}

std::vector<std::string> FileSystem::splitPath(std::string path) {
    path = cleanupPath(path);

    // split by goodSlash_
    std::vector<std::string> components;
    std::stringstream stream(path);
    std::string item;
    while(std::getline(stream, item, goodSlash_)) {
        components.push_back(item);
    }
    return components;
}

std::string FileSystem::currentDirectory() {
#ifdef WIN32
    char* buffer = new char[MAX_PATH + 1];
    memset(buffer, 0, MAX_PATH + 1);

    DWORD size = GetCurrentDirectory(MAX_PATH, buffer);
    if (size >= MAX_PATH) {
        delete [] buffer;
        buffer = new char[size + 1];
        memset(buffer, 0, size + 1);
        size = GetCurrentDirectory(size, buffer);
    }
    std::string dir(buffer);
    dir.resize(size);
    delete[] buffer;

    return dir;
#else
    char path[FILENAME_MAX];
    if (getcwd(path, FILENAME_MAX) == 0)
        return "";
    else
        return std::string(path);
#endif
}

bool FileSystem::changeDirectory(const std::string& directory) {
    if (directory.empty())
        return false;

    std::string converted = replaceAllCharacters(directory, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);
#ifdef WIN32
    return (_chdir(converted.c_str()) == 0);
#else
    return (chdir(converted.c_str()) == 0);
#endif
}

bool FileSystem::createDirectory(const std::string& directory) {
    if (directory.empty()) {
        return false;
    }

    std::string converted = replaceAllCharacters(directory, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

#ifdef WIN32
    return CreateDirectory(converted.c_str(), 0) != 0;
#else
    return mkdir(converted.c_str(), 0777) == 0;
#endif
}

bool FileSystem::createDirectoryRecursive(const std::string& directory) {
    if(dirExists(directory))
        return true;
    else {
        std::string dir = parentDir(directory);
        if(dir == directory) {
            // Reached root of (potentially relative path), assuming / or current directory exists.
            return true;
        }

        if(createDirectoryRecursive(dir))
            return createDirectory(directory);

        return false;
    }
}

bool FileSystem::deleteDirectory(const std::string& directory) {
    if (directory.empty())
        return false;

    std::string converted = replaceAllCharacters(directory, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);
#ifdef WIN32
    return (_rmdir(converted.c_str()) == 0);
#else
    return (rmdir(converted.c_str()) == 0);
#endif
}

bool FileSystem::deleteDirectoryRecursive(const std::string& directory) {
    if (directory.empty())
        return false;

    bool success = clearDirectory(directory);
    success &= deleteDirectory(directory);

    return success;
}

bool FileSystem::clearDirectory(const std::string& directory) {
    if (directory.empty())
        return false;

    std::string converted = replaceAllCharacters(directory, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

    bool success = true;

    //recursively delete all subdirectories
    std::vector<std::string> subDirs = listSubDirectories(converted);
    for(size_t i=0; i<subDirs.size(); i++) {
        success &= deleteDirectoryRecursive(converted + goodSlash_ + subDirs[i]);
    }
    //delete all files
    std::vector<std::string> files = listFiles(converted);
    for(size_t i=0; i<files.size(); i++) {
        success &= deleteFile(converted + goodSlash_ + files[i]);
    }

    return success;
}

bool FileSystem::deleteFile(const std::string& filename) {
    if (filename.empty())
        return false;

    std::string converted = replaceAllCharacters(filename, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

#ifdef WIN32
    return (DeleteFile(converted.c_str()) != 0);
#else
    return (remove(converted.c_str()) == 0);
#endif
}

bool FileSystem::renameFile(const std::string& filename, const std::string& newName,
                            bool ignorePath)
{
    if (FileSystem::fileExists(filename) == false)
        return false;

    std::string converted = replaceAllCharacters(filename, badSlash_, goodSlash_);
    std::string convertedNew = replaceAllCharacters(newName, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);
    removeTrailingCharacters(convertedNew, goodSlash_);

    bool res = false;
    if (ignorePath == true) {
        std::string name = (FileSystem::dirName(filename) + std::string(&goodSlash_, 1)
            + FileSystem::fileName(newName));
        res = (rename(converted.c_str(), name.c_str()) == 0);
    } else {
        res = (rename(converted.c_str(), convertedNew.c_str()) == 0);
    }
    return res;
}

void FileSystem::copyFile(const std::string& srcFile, const std::string& destFile) {
    // check if input file/output dir exists
    if (!fileExists(srcFile))
        throw tgt::FileNotFoundException("srcFile not found", srcFile);
    std::string destDir = dirName(destFile);
    if (destDir != "" && !dirExists(destDir))
        throw tgt::FileNotFoundException("destination directory not found", destDir);

    // open input and output files
    std::ifstream src(srcFile.c_str(), std::ios::in | std::ios::binary);
    if (!src.good())
        throw tgt::FileAccessException("failed to open srcFile fir reading", srcFile);
    std::ofstream dest(destFile.c_str(), std::ios::trunc | std::ios::binary);
    if (!dest.good())
        throw tgt::FileAccessException("failed to open dest file for writing", destFile);

    // copy data
    dest << src.rdbuf();

    // check if copying succeeded
    if (!src.good())
        throw tgt::FileException("reading from source file failed", srcFile);
    if (!dest.good())
        throw tgt::FileException("writing to dest file failed", destFile);
}

void FileSystem::updateFileTime(const std::string& filename) {
    if (!fileExists(filename))
        throw tgt::Exception("File does not exist: " + filename);

    struct utimbuf newTime;
    newTime.actime = time(NULL);       /* set atime to current time */
    newTime.modtime = time(NULL);      /* set mtime to current time */
    utime(filename.c_str(), &newTime);
}

bool FileSystem::fileExists(const std::string& filename) {
    if (filename.empty())
        return false;

    std::string converted = replaceAllCharacters(filename, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

#ifdef WIN32
    WIN32_FILE_ATTRIBUTE_DATA fad;
    return (GetFileAttributesEx(converted.c_str(), GetFileExInfoStandard, &fad) != 0);
#else
    struct stat st;
    return (stat(converted.c_str(), &st) == 0);
#endif
}

uint64_t FileSystem::dirSize(const std::string& directory, const bool recursive) {
    std::vector<std::string> files;
    if(recursive)
        files = listFilesRecursive(directory, false);
    else
        files = listFiles(directory, false);

    uint64_t size = 0;
    for(size_t i=0; i<files.size(); i++) {
        size += fileSize(directory + goodSlash_ + files[i]);
    }
    return size;
}

uint64_t FileSystem::fileSize(const std::string& filename) {
    if (filename.empty())
        return 0;

    std::string converted = replaceAllCharacters(filename, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

#ifdef WIN32
    WIN32_FILE_ATTRIBUTE_DATA fad;
    if (GetFileAttributesEx(converted.c_str(), GetFileExInfoStandard, &fad)) {
        LARGE_INTEGER size;
        size.HighPart = fad.nFileSizeHigh;
        size.LowPart = fad.nFileSizeLow;
        return size.QuadPart;
    }
#else
    struct stat st;
    if(stat(converted.c_str(), &st) == 0) {
        return st.st_size;
    }
#endif
    else {
        // error
        return 0;
    }
}

time_t FileSystem::fileTime(const std::string& filename) {
    if (filename.empty())
        return 0;

    std::string converted = replaceAllCharacters(filename, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

#ifdef WIN32
    WIN32_FILE_ATTRIBUTE_DATA fad;
    if (GetFileAttributesEx(converted.c_str(), GetFileExInfoStandard, &fad)) {
        LARGE_INTEGER size;
        size.HighPart = fad.ftLastWriteTime.dwHighDateTime;
        size.LowPart = fad.ftLastWriteTime.dwLowDateTime;
        return size.QuadPart;
    }
#else
    struct stat st;
    if(stat(converted.c_str(), &st) == 0) {
        return st.st_mtime;
    }
#endif
    else {
        // error
        return 0;
    }
}

time_t FileSystem::fileAccessTime(const std::string& filename) {
    if (filename.empty())
        return 0;

    std::string converted = replaceAllCharacters(filename, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

#ifdef WIN32
    WIN32_FILE_ATTRIBUTE_DATA fad;
    if (GetFileAttributesEx(converted.c_str(), GetFileExInfoStandard, &fad)) {
        LARGE_INTEGER size;
        size.HighPart = fad.ftLastAccessTime.dwHighDateTime;
        size.LowPart = fad.ftLastAccessTime.dwLowDateTime;
        return size.QuadPart;
    }
#else
    struct stat st;
    if(stat(converted.c_str(), &st) == 0) {
        return st.st_atime;
    }
#endif
    else {
        // error
        return 0;
    }
}

bool FileSystem::dirExists(const string& dirpath) {
#ifdef WIN32
    DWORD result = GetFileAttributes(dirpath.c_str());
    return (result != INVALID_FILE_ATTRIBUTES) && (result & FILE_ATTRIBUTE_DIRECTORY);
#else
    DIR* dir = opendir(dirpath.c_str());
    if (dir) {
        closedir(dir);
        return true;
    }
    else
        return false;
#endif
}

std::vector<std::string> FileSystem::readDirectory(const std::string& directory, const bool sort,
                                                   const bool recursiveSearch) {
    if(recursiveSearch)
        return listFilesRecursive(directory, sort);
    else
        return listFiles(directory, sort);
}

#ifdef WIN32
std::vector<std::string> FileSystem::listFiles(const std::string& directory, const bool sort) {
    std::vector<std::string> result;
    if (directory.empty())
        return result;

    std::string converted = replaceAllCharacters(directory, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

    WIN32_FIND_DATA findFileData = {0};
    HANDLE hFind = 0;

    std::stack<std::string> stackDirs;
    std::string dir(converted + "\\*");

    hFind = FindFirstFile(dir.c_str(), &findFileData);
    if (hFind != INVALID_HANDLE_VALUE) {
        do {
            std::string file(findFileData.cFileName);
            if (! (findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
                    result.push_back(file);
            }
        } while (FindNextFile(hFind, &findFileData) != 0);
    }
    FindClose(hFind);

    if (sort)
        std::sort(result.begin(), result.end());

    return result;
}
#else
std::vector<std::string> FileSystem::listFiles(const std::string& directory, const bool sort) {
    std::vector<std::string> result;

    // POSIX directory listing
    DIR *dir = NULL;
    struct dirent *ent;
    if ((dir= opendir(directory.c_str())) != NULL) {
        std::string name;
        while ((ent = readdir(dir))) {
            name = ent->d_name;
            if ((name != ".") && (name != "..")) {
                // This does not work on all filesystems!!! Therefore replaced by stat
                //if(ent->d_type != DT_DIR)
                //    result.push_back(ent->d_name);
                std::string completeName = directory + "/" + name;
                struct stat st;
                if(stat(completeName.c_str(), &st) == 0 && !(S_ISDIR(st.st_mode)))
                    result.push_back(name);
            }
        }
    }
    if (dir) {
        closedir(dir);
    }

    if (sort)
        std::sort(result.begin(), result.end());

    return result;
}
#endif

std::vector<std::string> FileSystem::listFilesRecursive(const std::string& directory, const bool sort) {
    std::vector<std::string> result;
    if (directory.empty())
        return result;

    result = listFiles(directory, false);

    std::vector<std::string> subDirs = listSubDirectories(directory, false);
    for(size_t i=0; i<subDirs.size(); i++) {
        std::vector<std::string> files = listFilesRecursive(directory + "/" + subDirs[i], false);
        for(size_t j=0; j<files.size(); j++) {
            result.push_back(subDirs[i] + "/" +files[j]);
        }
    }

    if (sort)
        std::sort(result.begin(), result.end());

    return result;
}

#ifdef WIN32
std::vector<std::string> FileSystem::listSubDirectories(const std::string& directory, const bool sort) {
    std::vector<std::string> result;
    if (directory.empty())
        return result;

    std::string converted = replaceAllCharacters(directory, badSlash_, goodSlash_);
    removeTrailingCharacters(converted, goodSlash_);

    WIN32_FIND_DATA findFileData = {0};
    HANDLE hFind = 0;

    std::stack<std::string> stackDirs;
    std::string dir(converted + "\\*");

    hFind = FindFirstFile(dir.c_str(), &findFileData);
    if (hFind != INVALID_HANDLE_VALUE) {
        do {
            std::string file(findFileData.cFileName);
            if ((file != ".") && (file != "..")) {
                if (findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                        result.push_back(file);
                }
            }
        } while (FindNextFile(hFind, &findFileData) != 0);
    }
    FindClose(hFind);

    if (sort)
        std::sort(result.begin(), result.end());

    return result;
}
#else
std::vector<std::string> FileSystem::listSubDirectories(const std::string& directory, const bool sort) {
    std::vector<std::string> result;

    // POSIX directory listing
    DIR *dir = NULL;
    struct dirent *ent;
    if ((dir= opendir(directory.c_str())) != NULL) {
        std::string name;
        while ((ent = readdir(dir))) {
            name = ent->d_name;
            if ((name != ".") && (name != "..")) {
                // This does not work on all filesystems!!! Therefore replaced by stat
                //if(ent->d_type == DT_DIR)
                //    result.push_back(ent->d_name);
                std::string completeName = directory + "/" + name;
                struct stat st;
                if(stat(completeName.c_str(), &st) == 0 && S_ISDIR(st.st_mode))
                    result.push_back(name);
            }
        }
    }
    if (dir) {
        closedir(dir);
    }


    if (sort)
        std::sort(result.begin(), result.end());

    return result;
}
#endif

std::vector<std::string> FileSystem::listSubDirectoriesRecursive(const std::string& directory, const bool sort) {
    std::vector<std::string> result;
    if (directory.empty())
        return result;

    result = listSubDirectories(directory, false);

    std::vector<std::string> subDirs = result;
    for(size_t i=0; i<subDirs.size(); i++) {
        std::vector<std::string> dirs = listSubDirectoriesRecursive(directory + "/" + subDirs[i], false);
        for(size_t j=0; j<dirs.size(); j++) {
            result.push_back(subDirs[i] + "/" + dirs[j]);
        }
    }

    if (sort)
        std::sort(result.begin(), result.end());

    return result;
}


} // namespace
