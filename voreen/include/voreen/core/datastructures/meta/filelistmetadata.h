/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_FILELISTMETADATA_H
#define VRN_FILELISTMETADATA_H

#include "voreen/core/io/serialization/serialization.h"
#include "templatemetadata.h"

namespace voreen {

///Metadata encapsulating a list of files
class VRN_CORE_API FileListMetaData : public TemplateMetaData<std::vector<std::string> > {
public:
    FileListMetaData() : TemplateMetaData<std::vector<std::string> >() {}
    FileListMetaData(std::vector<std::string> files) : TemplateMetaData<std::vector<std::string> >(files) {}

    virtual MetaDataBase* clone() const      { return new FileListMetaData(getValue()); }
    virtual std::string getClassName() const { return "FileListMetaData"; }
    virtual MetaDataBase* create() const { return new FileListMetaData(); }

    virtual std::string toString() const {
        std::vector<std::string>::const_iterator i;
        std::stringstream s;
        for (i = getValue().begin(); i != getValue().end(); ++i) {
            s << std::endl << *i;
        }
        return s.str();
    }
};

} // namespace

#endif
