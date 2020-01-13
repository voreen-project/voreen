/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#include "tgt/exception.h"

namespace tgt {

using std::string;

//------------------------------------------------------------------------------

const char* Exception::what() const throw() {
    return what_.c_str();
};

//------------------------------------------------------------------------------

std::string FileException::getErrorMessage(const string& exception, const string& what,
                                           const string& filename) const
    throw()
{
    string s = exception;

    if (!what.empty())
        s += ": " + what;

    if (!filename.empty())
        s += " (" + filename + ")";

    return s;
}

//------------------------------------------------------------------------------

FileNotFoundException::FileNotFoundException(const std::string& what, const std::string& filename)
    : FileException(getErrorMessage("FileNotFoundException", what, filename), filename)
{}

//------------------------------------------------------------------------------

FileAccessException::FileAccessException(const std::string& what, const std::string& filename)
    : FileException(getErrorMessage("FileAccessException", what, filename), filename)
{}

//------------------------------------------------------------------------------

CorruptedFileException::CorruptedFileException(const std::string& what, const std::string& filename)
    : FileException(getErrorMessage("CorruptedFileException", what, filename), filename)
{}

//------------------------------------------------------------------------------

UnsupportedFormatException::UnsupportedFormatException(const std::string& extension, const std::string& filename)
    : FileException(getErrorMessage("UnsupportedFormatException", "unsupported extension: " + extension, filename),
                    filename)
{}

//------------------------------------------------------------------------------

IOException::IOException(const std::string& what, const std::string& filename)
    : FileException(getErrorMessage("IOException", what, filename), filename)
{}

} // namespace tgt
