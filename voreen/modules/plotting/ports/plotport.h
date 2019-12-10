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

#ifndef VRN_PLOTPORT_H
#define VRN_PLOTPORT_H

#include "voreen/core/ports/genericport.h"
#include "../datastructures/plotbase.h"

namespace voreen {

/**
 * \brief This class describes the Port who transports PlotBase Objects
 */
class VRN_CORE_API PlotPort : public GenericPort<PlotBase> {
public:

    PlotPort(PortDirection direction, const std::string& name, const std::string& guiName = "",
                        bool allowMultipleConnections = false,
                        Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);
    ~PlotPort();

    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new PlotPort(direction,id,guiName);}
    virtual std::string getClassName() const {return "PlotPort";}

    /**
     * Assigns the passed PlotBase handle to the port and
     * deletes the currently assigned one, if deletePrevious
     * is specified.
     */
    void setData(const PlotBase* handle, bool deletePrevious = false);

    /**
     * Deletes the assigned PlotBase handle and
     * assigns the null pointer to itself.
     */
    void deleteData();

    /**
     * Returns true, if the port contains a PlotBaseHandle object
     * and the PlotBaseHandle has a valid PlotBase.
     */
     virtual bool isReady() const;

};

} // namespace

#endif // VRN_PLOTPORT_H
