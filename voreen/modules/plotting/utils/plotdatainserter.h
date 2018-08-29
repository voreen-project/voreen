/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTDATAINSERTER_H
#define VRN_PLOTDATAINSERTER_H

#include "../datastructures/plotdata.h"
#include "../datastructures/plotcell.h"

#include <string>
#include <vector>
#include <set>

namespace voreen {

/**
 * Helper class offering easier access to PlotData::insert methods.
 **/
class VRN_CORE_API PlotDataInserter {
public:
    /**
     * Inserts a row containing \a arg1 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     **/
    static bool insert(PlotData& pData, const PlotCellValue& arg1);

    /**
     * Inserts a row containing \a arg1 to \a arg2 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     **/
    static bool insert(PlotData& pData, const PlotCellValue& arg1, const PlotCellValue& arg2);

    /**
     * Inserts a row containing \a arg1 to \a arg3 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg3    data in third column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3);

    /**
     * Inserts a row containing \a arg1 to \a arg4 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg3    data in third column
     * \param   arg4    data in fourth column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4);

    /**
     * Inserts a row containing \a arg1 to \a arg5 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg3    data in third column
     * \param   arg4    data in fourth column
     * \param   arg5    data in fifth column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4,
                       const PlotCellValue& arg5);

    /**
     * Inserts a row containing \a arg1 to \a arg6 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg3    data in third column
     * \param   arg4    data in fourth column
     * \param   arg5    data in fifth column
     * \param   arg6    data in sixth column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4,
                       const PlotCellValue& arg5,
                       const PlotCellValue& arg6);

    /**
     * Inserts a row containing \a arg1 to \a arg7 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg4    data in fourth column
     * \param   arg5    data in fifth column
     * \param   arg6    data in sixth column
     * \param   arg7    data in seventh column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4,
                       const PlotCellValue& arg5,
                       const PlotCellValue& arg6,
                       const PlotCellValue& arg7);

    /**
     * Inserts a row containing \a arg1 to \a arg8 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg4    data in fourth column
     * \param   arg5    data in fifth column
     * \param   arg6    data in sixth column
     * \param   arg7    data in seventh column
     * \param   arg8    data in eighth column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4,
                       const PlotCellValue& arg5,
                       const PlotCellValue& arg6,
                       const PlotCellValue& arg7,
                       const PlotCellValue& arg8);

    /**
     * Inserts a row containing \a arg1 to \a arg9 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg4    data in fourth column
     * \param   arg5    data in fifth column
     * \param   arg6    data in sixth column
     * \param   arg7    data in seventh column
     * \param   arg8    data in eighth column
     * \param   arg9    data in ninth column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4,
                       const PlotCellValue& arg5,
                       const PlotCellValue& arg6,
                       const PlotCellValue& arg7,
                       const PlotCellValue& arg8,
                       const PlotCellValue& arg9);

    /**
     * Inserts a row containing \a arg1 to \a arg10 into \a pData.
     *
     * \param   pData   target PlotData for insertion
     * \param   arg1    data in first column
     * \param   arg2    data in second column
     * \param   arg3    data in third column
     * \param   arg4    data in fourth column
     * \param   arg5    data in fifth column
     * \param   arg6    data in sixth column
     * \param   arg7    data in seventh column
     * \param   arg8    data in eighth column
     * \param   arg9    data in ninth column
     * \param   arg10   data in tenth column
     **/
    static bool insert(PlotData& pData,
                       const PlotCellValue& arg1,
                       const PlotCellValue& arg2,
                       const PlotCellValue& arg3,
                       const PlotCellValue& arg4,
                       const PlotCellValue& arg5,
                       const PlotCellValue& arg6,
                       const PlotCellValue& arg7,
                       const PlotCellValue& arg8,
                       const PlotCellValue& arg9,
                       const PlotCellValue& arg10);

private:
    // no instaces of this static class allowed
    PlotDataInserter() {};
};

} // namespace voreen

#endif // VRN_PLOTDATA_H
