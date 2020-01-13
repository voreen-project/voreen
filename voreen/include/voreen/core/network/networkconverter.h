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

#ifndef VRN_NETWORKCONVERTER_H
#define VRN_NETWORKCONVERTER_H

#include <string>

class TiXmlElement;
class TiXmlNode;

namespace voreen {

/**
 * This abstract class serves as the base of several incrementally working NetworkConverters. Each NetworkConverter should change
 * the Network only by one version. This way it isn't necessary to write n new Converters if a network change occurs, but
 * only one - from the latest but one to the now latest.</br>
 * In this class there are some basic tools necessary to manipulate the TiXmlElement itself. A subclass should only
 * overwrite the NetworkConverter::convert(TiXmlElement*) method and don't add any other public methods.</br>
 * some outside method should care about specific version numbers and which NetworkConverter's method is to be called.</br>
 * Note: As far as this system is implemented now, the NetworkConverter doesn't check if the version is really the one
 * it is meant for. So it might get jumbled up if the supervising class makes an error.</br>
 */
class NetworkConverter {
public:
    virtual ~NetworkConverter() {}

    /**
     * <b>The</b> method of this class. It takes the whole document and manipulates it, so it is viable for the new version it is
     * intended for. Increase the counter 'NETWORK_VERSION' in processornetwork.cpp to the new version, so the TiXmlElement will be consistent
     * \param elem The TiXmlElement which should be converted to a new version
     */
    virtual void convert(TiXmlElement* elem) = 0;

protected:
    void changePropertyType(TiXmlElement* elem, const std::string& from, const std::string& to);
    void changePropertyName(TiXmlElement* elem, const std::string& ofProcessor, const std::string& from, const std::string& to);
    void changePropertyValue(TiXmlElement* elem, const std::string& ofProcessor, const std::string& from, const std::string& to);

    void changeProcessorType(TiXmlElement* elem, const std::string& from, const std::string& to);

    void changeLinkEvaluatorType(TiXmlElement* elem, const std::string& from, const std::string& to);
    static const std::string loggerCat_; ///< The logging information used in LINFO etc.
};

//-------------------------------------------------------------------------------------------------

/**
 * The @c NetworkConverter4to5 converts the serialized XML data for a
 * version 4 @c ProcessorNetwork to XML data for a version 5 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter4to5 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * The @c NetworkConverter5to6 converts the serialized XML data for a
 * version 5 @c ProcessorNetwork to XML data for a version 6 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter5to6 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * The @c NetworkConverter6to7 converts the serialized XML data for a
 * version 6 @c ProcessorNetwork to XML data for a version 7 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter6to7 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * The @c NetworkConverter7to8 converts the serialized XML data for a
 * version 7 @c ProcessorNetwork to XML data for a version 8 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter7to8 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * The @c NetworkConverter8to9 converts the serialized XML data for a
 * version 8 @c ProcessorNetwork to XML data for a version 9 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter8to9 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * The @c NetworkConverter9to10 converts the serialized XML data for a
 * version 9 @c ProcessorNetwork to XML data for a version 10 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter9to10 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * The @c NetworkConverter10to11 converts the serialized XML data for a
 * version 10 @c ProcessorNetwork to XML data for a version 11 @c ProcessorNetwork.
 *
 * @see NetworkConverter
 */
class NetworkConverter10to11 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * Moves the VolumeContainer entries to the references VolumeSources
 */
class NetworkConverter11to12 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
    void convertVolumeContainer(TiXmlElement* workspaceNode);
};

/**
 * Converts unary- and binaryimageprocessor to use shader property in parent classe.
 */
class NetworkConverter12to13 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * PropertyOwner uses attribute 'mapKey' and not 'name' as before.
 */
class NetworkConverter13to14 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * VolumeCollection -> VolumeList renaming.
 */
class NetworkConverter14to15 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * OptimizedProxyGeometry modes renaming.
 */
class NetworkConverter15to16 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * New Converter for Voreen 4.5
 * - Rename tf evaluators and properties
 * - rename meta data
 */
class NetworkConverter16to17 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * StreamlineCreator:
 * - Vec2Property -> IntervalProperty
 * - minNumStreamlinesPerBundle: absolute -> relative
 */
class NetworkConverter17to18 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * New Converter for Voreen 5.0
 * - Remove 'none' classification modus
 * - ...
 */
class NetworkConverter18to19 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

/**
 * Converter for Voreen 5.1
 * - Make OctreeCreator inherit AsyncComputeProcessor
 */
class NetworkConverter19to20 : public NetworkConverter {
public:
    virtual void convert(TiXmlElement* networkElem);
};

} // namespace

#endif // VRN_NETWORKCONVERTER_H
