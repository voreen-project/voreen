/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_XMLSERIALIZATIONCONSTANTS_H
#define VRN_XMLSERIALIZATIONCONSTANTS_H

#include "voreen/core/voreencoreapi.h"
#include <string>
#include <iostream>

namespace voreen {

/**
 * The @c XmlSerializationConstants class defines serialization constants, which are shared
 * by @c XmlSerializer and @c XmlDeserializer.
 *
 * @see XmlSerializer
 * @see XmlDeserializer
 */
class VRN_CORE_API XmlSerializationConstants {
public:
    /**
     * XML version attribute in XML declaration.
     */
    static const std::string XMLVERSION;

    /**
     * XML encoding attribute in XML declaration.
     */
    static const std::string XMLENCODING;

    /**
     * XML standalone attribute in XML declaration.
     */
    static const std::string XMLSTANDALONE;

    /**
     * XML root node name of serialized XML document.
     */
    static const std::string ROOTNODE;

    /**
     * Version attribute name.
     *
     * @note To obtain finer version control, you can set a version attribute for every XML node.
     */
    static const std::string VERSIONATTRIBUTE;

    /**
     * Version of @c XmlSerializer and @c XmlDeserializer.
     */
    static const std::string VERSION;

    /**
     * Id attribute name for reference resolving purposes.
     */
    static const std::string IDATTRIBUTE;

    /**
     * Reference attribute name for reference resolving purposes.
     */
    static const std::string REFERENCEATTRIBUTE;

    /**
     * Type attribute name for polymorphism purposes.
     */
    static const std::string TYPEATTRIBUTE;

    /**
     * Name of the value attribute which contains the content of simple types.
     */
    static const std::string VALUEATTRIBUTE;

    /**
     * Id prefix for reference resolving purposes.
     */
    static const std::string IDPREFIX;

    /**
     * Name of item nodes for collection serialization.
     */
    static const std::string ITEMNODE;

    /**
     * Name of key nodes for @c key/value pair serialization of maps.
     */
    static const std::string KEYNODE;

    /**
     * Name of value nodes for @c key/value pair serialization of maps.
     */
    static const std::string VALUENODE;

    /**
     * Name of temporary XML nodes.
     */
    static const std::string TEMPNODE;

    /**
     * Attribute name of vector's x element.
     */
    static const std::string VECTORXATTRIBUTE;

    /**
     * Attribute name of vector's y element.
     */
    static const std::string VECTORYATTRIBUTE;

    /**
     * Attribute name of vector's z element.
     */
    static const std::string VECTORZATTRIBUTE;

    /**
     * Attribute name of vector's w element.
     */
    static const std::string VECTORWATTRIBUTE;

    /**
     * Attribute name of color's r element.
     */
    static const std::string COLORRATTRIBUTE;

    /**
     * Attribute name of color's g element.
     */
    static const std::string COLORGATTRIBUTE;

    /**
     * Attribute name of color's b element.
     */
    static const std::string COLORBATTRIBUTE;

    /**
     * Attribute name of color's a element.
     */
    static const std::string COLORAATTRIBUTE;
};

} // namespace

#endif // VRN_XMLSERIALIZATIONCONSTANTS_H
