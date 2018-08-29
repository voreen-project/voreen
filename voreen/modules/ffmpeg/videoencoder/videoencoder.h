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

#ifndef VRN_VIDEOENCODER_H
#define VRN_VIDEOENCODER_H

#include "tgt/exception.h"
#include "tgt/tgt_gl.h"

#include "voreen/core/voreencoreapi.h"

#include <string>
#include <vector>

namespace voreen {

/**
 *
 * <emph>VideoEncoder</emph> is ready to receive frames by {@link #nextFrame(GLvoid*)}
 *  after call of {@link #startVideoEncoding(char*, unsigned int, unsigned int, unsigned int)}<br />
 * calling {@link #stopVideoEncoding()} will complete the encoding process (in any case)
 *
 * <h3>references:</h3>
 * <ul>
 * <li><a href="http://www.ffmpeg.com.cn/index.php/Output_example.c#write_video_frame">http://www.ffmpeg.com.cn/index.php/Output_example.c#write_video_frame</a></li>
 * </ul>
 */
class VRN_CORE_API VideoEncoder {
public:
   /**
    * @see videoencoder.cpp#containerCodecPairNames
    * @see videoencoder.cpp#containerAppendix
    * @see videoencoder.cpp#containerCodecId
    * FIXME: WMVWMV crashes at av_set_parameter
    * FIXME: OGGTHEORA results in 0 second video
    */
    enum ContainerCodecPair {
        GUESS, // let ffmpeg guess by filename or default to mpeg4
        MPEG4AVI,
        WMVWMV,
        FLVFLV,
        HUFFYUVAVI,
        OGGTHEORA,
        LAST // dummy
    };

    VideoEncoder();
    ~VideoEncoder();

    /**
     * @param filePath target video file's full path (absolute or relative)
     * @param fps target video's fps
     * @param width source=target-video's width
     * @param height video's height
     * @param pixelFormat input pixel format, supported: GL_RGBA
     * @param pixelType input pixel type, supported: GL_FLOAT, GL_UNSIGNED_SHORT, GL_UNSIGNED_BYTE
     */
    void startVideoEncoding(std::string filePath, const int fps,
                            const int width, const int height,
                            GLint pixelFormat, GLenum pixelType);
    /**
     * finish encoding
     */
    void stopVideoEncoding();

    /**
     * @param pixels passed frame must match <emph>width</emph>, <emph>height</emph> and <emph>pixelFormat</emph>, <emph>-Type</emph><br />
     * latest two assumed as GL_RGBA and GL_FLOAT
     */
    void nextFrame(GLvoid* pixels);

    void setup(int preset, int bitrate);

    ContainerCodecPair getPreset(){
        return preset_;
    }

    const char* getContainerAppendix();

    int getBitrate(){
        return bitrate_;
    }

    static std::vector<std::string> getSupportedFormatsByFileEnding();
    static std::vector<std::string> getSupportedFormatDescriptions();
    static const char** getContainerCodecPairNames();

protected:
    /**
     * tgt logging category
     */
    static std::string loggerCat_;

    int bitrate_;
    ContainerCodecPair preset_;
};

} // namespace

#endif // VRN_VIDEOENCODER_H

