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

#ifndef VRN_PARALLELCOORDINATESSECTIONSPROPERTY_H
#define VRN_PARALLELCOORDINATESSECTIONSPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/link/linkevaluatorid.h"

namespace voreen {

struct ParallelCoordinatesSectionsPropertyData {
    std::string ensembleHash;
    tgt::Bounds bounds;
    std::string member;
    float time; // Time has a value between 0 and 1 that maps linearly to the common time domain of the ensemble.
    std::vector<std::pair<std::string, int>> fields;
    std::vector<std::list<tgt::vec2>> sections;

    ParallelCoordinatesSectionsPropertyData() = default;
    ParallelCoordinatesSectionsPropertyData( std::string ensembleHash,
                                             tgt::Bounds bounds,
                                             std::string member,
                                             float time,
                                             std::vector<std::pair<std::string, int>> fields,
                                             std::vector<std::list<tgt::vec2>> sections )
        : ensembleHash(std::move(ensembleHash))
        , bounds(bounds)
        , member(std::move(member) )
        , time( time )
        , fields( std::move( fields ) )
        , sections( std::move( sections ) )
    {}

    bool operator!=( const ParallelCoordinatesSectionsPropertyData& other ) const {
        return (ensembleHash != other.ensembleHash) || (bounds != other.bounds) || ( member != other.member ) || ( time != other.time ) || ( fields != other.fields ) || ( sections != other.sections );
    }
};

class ParallelCoordinatesSectionsProperty : public TemplateProperty<ParallelCoordinatesSectionsPropertyData> {
public:
    ParallelCoordinatesSectionsProperty( const std::string& id, const std::string& guiText, int invalidationLevel = Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT );
    ParallelCoordinatesSectionsProperty() = default;
    virtual ~ParallelCoordinatesSectionsProperty() = default;

    virtual Property* create() const override;
    virtual std::string getClassName() const override;
    virtual std::string getTypeDescription() const override;
};

class LinkEvaluatorParallelCoordinatesSectionsId : public LinkEvaluatorIdGeneric<ParallelCoordinatesSectionsPropertyData> {
public:
    virtual std::string getClassName() const override { return "LinkEvaluatorParallelCoordinatesSectionsId"; }
    virtual LinkEvaluatorBase* create() const override { return new LinkEvaluatorParallelCoordinatesSectionsId(); }
};

}

#endif // VRN_PARALLELCOORDINATESSECTIONSPROPERTY_H