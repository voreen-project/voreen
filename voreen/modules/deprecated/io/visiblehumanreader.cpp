/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "visiblehumanreader.h"

#include <iostream>
#include <fstream>

namespace voreen {

    VisibleHumanReader::VisibleHumanReader() {
        colorDataset_ = true;
    }

    void VisibleHumanReader::createDataset(std::vector<std::string> filenames, std::string outputFile, int byteSkip, bool color) {
        std::fstream outputStream(outputFile.c_str(), std::ios::out |std::ios::binary);

        if (color == true) {

            char*** data = new char**[1216];
            for (int i=0;i<1216; i++) {
                data[i] = new char*[2048];
                for (int j=0; j<2048; j++) {
                    data[i][j] = new char[3];
                }
            }

            for (size_t i=0; i<filenames.size(); i++) {

                std::fstream inputStream(filenames.at(i).c_str(), std::ios::in |std::ios::binary);
                std::cout << "reading: " << filenames.at(i) << std::endl;

                char* charArrayBuffer = new char[1];

                int cutLeft,cutRight,cutTop,cutBottom;

                cutLeft = 0;
                cutRight=0;
                cutTop=0;
                cutBottom=0;

                int sizex = 512;
                int sizey = 512;

                for (int k=0;k<3; k++) {
                    for (int i=0; i<sizey;i++) {
                        for (int j=0; j<sizex; j++) {
                            inputStream.read(charArrayBuffer,1);
                            data[i][j][k] = charArrayBuffer[0];
                        }
                    }
                }

                for (int i=cutTop; i<sizey-cutBottom; i++) {
                    for (int j=cutLeft; j<sizex-cutRight; j++) {
                        for (int k=0; k<3; k++) {
                            charArrayBuffer[0] = data[i][j][k];
                            outputStream.write(charArrayBuffer,1);
                        }
                    }
                }

                inputStream.close();
                delete[] charArrayBuffer;
            }
            outputStream.close();
            delete[] data;
        }

        else {
            char* charArrayBuffer = new char[2];
            for (size_t i=0; i<filenames.size(); i++) {
                std::fstream inputStream(filenames.at(i).c_str(), std::ios::in |std::ios::binary);
                std::cout << "reading: " << filenames.at(i) << std::endl;
                if (byteSkip != 0)
                    inputStream.seekg(byteSkip+1,std::ios::beg);
                long max=512*512;
                for (long j=0; j<max; j++) {
                    inputStream.read(charArrayBuffer, 2);
                    outputStream.write(charArrayBuffer,2);
                }

                inputStream.close();
            }
            delete[] charArrayBuffer;
            outputStream.close();
        }

    }

} //namespace voreen
