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

#ifndef VRN_PORTCONNECTION_H
#define VRN_PORTCONNECTION_H

#include "voreen/core/io/serialization/serializable.h"

namespace voreen {

class Port;

/**
 * The @c PortConnection is responsible for serialization and deserialization
 * of connections between different processor ports.
 *
 * @see Serializable
 */
class VRN_CORE_API PortConnection : public Serializable {
public:
    /**
     * Creates a @c PortConnection from given outport to given inport.
     *
     * @param outport the outport
     * @param inport the inport
     */
    PortConnection(Port* outport, Port* inport);

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Sets the outport of the connection.
     *
     * @param value the outport
     */
    void setOutport(Port* value);

    /**
     * Returns the outport of the connection.
     *
     * @return the outport
     */
    Port* getOutport() const;

    /**
     * Sets the inport of the connection.
     *
     * @param value the inport
     */
    void setInport(Port* value);

    /**
     * Returns the inport of the connection.
     *
     * @return the inport
     */
    Port* getInport() const;

private:
    /**
     * The @c PortEntry class processes the port data
     * to get a well readable XML file.
     *
     * @see PortConnection
     * @see Serializable
     */
    class PortEntry : public Serializable {
    public:
        /**
         * Creates a @c PortEntry for given @c Port.
         *
         * @param port the port
         */
        PortEntry(Port* port);

        /**
        * @see Serializable::serialize
        */
        virtual void serialize(Serializer& s) const;

        /**
        * @see Serializable::deserialize
        */
        virtual void deserialize(Deserializer& s);

        /**
         * Returns the port which data will be processed.
         *
         * @return the port
         */
        Port* getPort() const;

    private:
        /**
         * Port which data will be processed.
         */
        Port* port_;
    };

    friend class Deserializer;

    /**
     * Default constructor for serialization purposes.
     */
    PortConnection();

    /**
     * Outport entry of connection.
     */
    PortEntry outport_;

    /**
     * Inport entry of connection.
     */
    PortEntry inport_;
};

} // namespace

#endif // VRN_PORTCONNECTION_H
