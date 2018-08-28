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

#include "videoencoder.h"

#include "tgt/logmanager.h"
#include "tgt/vector.h"

#include <sstream>

extern "C" {
#include <libavformat/avformat.h>
#include <libavformat/avio.h>
#include <libavutil/mathematics.h>
#include <libswscale/swscale.h>
}

namespace voreen {

namespace {

enum PixelOrderManipulation {
    NONE, FLIP_VERT, REVERSE
};

/**
 * @see ContainerCodecPair
 * @see VideoEncoder#getContainerCodecPairNames()
 */
const char *containerCodecPairNames[] = { "auto", "mpeg4 in avi", "wmv in wmv",
        "flv in flv", "huffyuv in avi", "ogg in ogg" };

/**
 * float to uint8_t conversion
 * TODO we have to put down an assumption about the pixel-order expected or let it be specified for the input data
 * @param pixels as floats
 * @param size number of pixels
 * @param lineSize (ignored but for FLIP_VERT)
 * @param manip (defaults: NONE)
 * @return pixels of type uint8_t
 */
uint8_t* floatToUint8_tRGB(float* pixels, const int size, const int lineSize,
        enum PixelOrderManipulation manip = NONE) {
    tgt::vec4* pixels_f = reinterpret_cast<tgt::vec4*> (pixels);
    tgt::col3* pixels_b = new tgt::col3[size];
    int j = 0;
    for (int i = 0; i < size; i++) {
        /*
         * 0 ... lineSize-1
         * lineSize ... 2*lineSize-1
         * ...
         * size-lineSize ... size-1
         */
        switch (manip) {
        case NONE:
            j = i; // up to compiler optimisation
            break;
        case REVERSE:
            j = size - i - 1;
            break;
        case FLIP_VERT:
            int curLine = (size / lineSize) - (i / lineSize) - 1;
            j = curLine * lineSize + (i % lineSize); // assume int div truncates
            break;
        }

        // dropping alpha channel and clamp to [0.0; 1.0]
        pixels_b[i] = tgt::clamp(pixels_f[j].xyz(), 0.f, 1.f) * 255.f;
    }
    return reinterpret_cast<uint8_t*> (pixels_b);
}

/**
 * uint16_t to uint8_t conversion
 * TODO we have to put down an assumption about the pixel-order expected or let it be specified for the input data
 * @param pixels as uint16_t
 * @param size number of pixels
 * @param lineSize (ignored but for FLIP_VERT)
 * @param manip (defaults: NONE)
 * @return pixels of type uint8_t
 */
uint8_t* uint16ToUint8_tRGB(uint16_t* pixels, const int size, const int lineSize,
        enum PixelOrderManipulation manip = NONE) {
    tgt::Vector4<uint16_t>* pixels_16 = reinterpret_cast<tgt::Vector4<uint16_t>*> (pixels);
    tgt::col3* pixels_b = new tgt::col3[size];
    int j = 0;
    for (int i = 0; i < size; i++) {
        /*
         * 0 ... lineSize-1
         * lineSize ... 2*lineSize-1
         * ...
         * size-lineSize ... size-1
         */
        switch (manip) {
        case NONE:
            j = i; // up to compiler optimisation
            break;
        case REVERSE:
            j = size - i - 1;
            break;
        case FLIP_VERT:
            int curLine = (size / lineSize) - (i / lineSize) - 1;
            j = curLine * lineSize + (i % lineSize); // assume int div truncates
            break;
        }

        // dropping alpha channel and converting to 8 bit
        pixels_b[i] = tgt::col3(pixels_16[j].x >> 8, pixels_16[j].y >> 8, pixels_16[j].z >> 8);
    }
    return reinterpret_cast<uint8_t*> (pixels_b);
}

/**
 * uint8_t RGBA to uint8_t RGB conversion
 * TODO we have to put down an assumption about the pixel-order expected or let it be specified for the input data
 * @param pixels as uint8_t
 * @param size number of pixels
 * @param lineSize (ignored but for FLIP_VERT)
 * @param manip (defaults: NONE)
 * @return pixels of type uint8_t
 */
uint8_t* uint8ToUint8_tRGB(uint8_t* pixels, const int size, const int lineSize,
        enum PixelOrderManipulation manip = NONE) {
    tgt::col4* pixels_col = reinterpret_cast<tgt::col4*> (pixels);
    tgt::col3* pixels_b = new tgt::col3[size];
    int j = 0;
    for (int i = 0; i < size; i++) {
        /*
         * 0 ... lineSize-1
         * lineSize ... 2*lineSize-1
         * ...
         * size-lineSize ... size-1
         */
        switch (manip) {
        case NONE:
            j = i; // up to compiler optimisation
            break;
        case REVERSE:
            j = size - i - 1;
            break;
        case FLIP_VERT:
            int curLine = (size / lineSize) - (i / lineSize) - 1;
            j = curLine * lineSize + (i % lineSize); // assume int div truncates
            break;
        }

        // dropping alpha channel
        pixels_b[i] = tgt::col3(pixels_col[j].x, pixels_col[j].y, pixels_col[j].z);
    }
    return reinterpret_cast<uint8_t*> (pixels_b);
}

/**
 * size of encoder's output-buffer and half of IO-buffer's size
 */
unsigned int videoOutBufferSize;

/*
 * some default codec parameters commonly used
 */
const int defaultGopSize = 10;
const int defaultMaxBFrame = 1;
const int minimalFps = 1;
const int minimalBitrate = 400000;

/*
 * containerAppendices
 * @see tgt::ContainerCodecPair
 * @see VEAVContext#initFormatContext(tgt::ContainerCodecPair)
 */
const char* mpegAppendix = "mpg";
const char* aviAppendix = "avi";
const char* flvAppendix = "flv";
const char* wmvAppendix = "asf";
const char* oggAppendix = "ogg";

const char* containerAppendix(VideoEncoder::ContainerCodecPair preset) {
    switch (preset) {
    case VideoEncoder::MPEG4AVI:
    case VideoEncoder::HUFFYUVAVI:
        return aviAppendix;
    case VideoEncoder::WMVWMV:
        return wmvAppendix;
    case VideoEncoder::FLVFLV:
        return flvAppendix;
    case VideoEncoder::OGGTHEORA:
        return oggAppendix;
    default:
        return mpegAppendix;
    }
}

#if (LIBAVFORMAT_VERSION_MAJOR < 54)
CodecID containerCodecId(VideoEncoder::ContainerCodecPair preset) {
    switch (preset) {
    default: //TODO: why??? crash in windows
        return CODEC_ID_MPEG2VIDEO;
    case VideoEncoder::FLVFLV:
        return CODEC_ID_FLV1;
    case VideoEncoder::MPEG4AVI:
        return CODEC_ID_MPEG4;
    case VideoEncoder::WMVWMV:
        return CODEC_ID_WMV2;
    case VideoEncoder::HUFFYUVAVI:
        return CODEC_ID_HUFFYUV;
    case VideoEncoder::OGGTHEORA:
        return CODEC_ID_THEORA;
    }
}
#else
AVCodecID containerCodecId(VideoEncoder::ContainerCodecPair preset) {
    switch (preset) {
    default: //TODO: why??? crash in windows
        return AV_CODEC_ID_MPEG2VIDEO;
    case VideoEncoder::FLVFLV:
        return AV_CODEC_ID_FLV1;
    case VideoEncoder::MPEG4AVI:
        return AV_CODEC_ID_MPEG4;
    case VideoEncoder::WMVWMV:
        return AV_CODEC_ID_WMV2;
    case VideoEncoder::HUFFYUVAVI:
        return AV_CODEC_ID_HUFFYUV;
    case VideoEncoder::OGGTHEORA:
        return AV_CODEC_ID_THEORA;
    }
}
#endif

AVPixelFormat codecPixFmt(VideoEncoder::ContainerCodecPair preset) {
    switch (preset) {
    default:
        return AV_PIX_FMT_YUV420P;
    case VideoEncoder::HUFFYUVAVI:
        return AV_PIX_FMT_YUV422P;
    }
}

struct VEAVContext {

    /**
     * log category for tgt-logging
     */
    std::string loggerCat_;

    // Encoder format-configuration and stream-state
    AVFormatContext* formatContext_;
    AVStream* videoStream;

    SwsContext* swsContext; // converts inputFrame to videoFrame, which is ready for encoding
    AVFrame* videoFrame; // encoded with videoFormat pix_fmt

    int frameCount;

    /**
     * encoder's output buffer
     * @see #videoOutBufferSize
     */
    uint8_t* videoBuffer;

    GLenum pixelType_; // input specification
    GLint pixelFormat_;

    /**
     * init picture- & video-buffer
     * @see #initCodec(int,int,int,int) should be run before
     */
    bool initBuffers() {
        AVCodecContext* codecContext = videoStream->codec;
        AVPixelFormat pixelFormat = codecContext->pix_fmt;
        const int width = codecContext->width;
        const int height = codecContext->height;

        // alloc frame/picture buffer
#if LIBAVCODEC_VERSION_INT >= AV_VERSION_INT(55, 45, 101)
        videoFrame = av_frame_alloc();
#else
        videoFrame = avcodec_alloc_frame();
#endif

        // HACK
#if (LIBAVFORMAT_VERSION_MAJOR < 54)
        if(formatContext_->video_codec_id == CODEC_ID_HUFFYUV || formatContext_->video_codec_id == CODEC_ID_THEORA)
#else
        if(formatContext_->video_codec_id == AV_CODEC_ID_HUFFYUV || formatContext_->video_codec_id == AV_CODEC_ID_THEORA)
#endif
            videoFrame->pts = 0;

        videoOutBufferSize = avpicture_get_size(pixelFormat, width, height);
        uint8_t* pictureBuffer = static_cast<uint8_t*> (av_malloc(
                videoOutBufferSize));
        if (!pictureBuffer) {
            LERRORC(loggerCat_, "unable to alloc picture buffer");
            return false;
        }
        avpicture_fill(reinterpret_cast<AVPicture*> (videoFrame),
                pictureBuffer, pixelFormat, width, height);

        // alloc video buffer
        videoBuffer = static_cast<uint8_t*> (av_malloc(videoOutBufferSize));

        return true;
    }

    /**
     * init video codec
     */
    bool initCodec(VideoEncoder::ContainerCodecPair preset, const int bitRate,
                   const int width, const int height, const int fps)
    {
        AVCodecContext* codecContext = videoStream->codec;
        codecContext->codec_id = formatContext_->video_codec_id;
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        codecContext->codec_type = CODEC_TYPE_VIDEO;
#else
        codecContext->codec_type = AVMEDIA_TYPE_VIDEO;
#endif

        codecContext->bit_rate = bitRate;
        codecContext->width = width;
        codecContext->height = height;

        codecContext->time_base.num = 1;
        codecContext->time_base.den = fps;

        codecContext->gop_size = defaultGopSize;
        codecContext->pix_fmt = codecPixFmt(preset);
#if (LIBAVFORMAT_VERSION_MAJOR < 54)
        if (codecContext->codec_id == CODEC_ID_MPEG2VIDEO) {
#else
        if (codecContext->codec_id == AV_CODEC_ID_MPEG2VIDEO) {
#endif
            codecContext->max_b_frames = defaultMaxBFrame;
        }
#if (LIBAVFORMAT_VERSION_MAJOR < 54)
        if (codecContext->codec_id == CODEC_ID_MPEG1VIDEO) {
#else
        if (codecContext->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
#endif
            codecContext->mb_decision = FF_MB_DECISION_RD;
        }
#if (LIBAVFORMAT_VERSION_MAJOR < 54)
        if (codecContext->codec_id == CODEC_ID_MPEG4) {
#else
        if (codecContext->codec_id == AV_CODEC_ID_MPEG4) {
#endif
            const char xvid[] = "DX50";
            codecContext->codec_tag = (xvid[3] << 24) + (xvid[2] << 16)
                    + (xvid[1] << 8) + xvid[0];
        }

        codecContext->strict_std_compliance = FF_COMPLIANCE_VERY_STRICT; // we are the good ones

        AVCodec* codec = avcodec_find_encoder(codecContext->codec_id);

        if (!codec) {
            LERRORC(loggerCat_, "Could not find Codec " << codecContext->codec_id);
            return false;
        }

#if (LIBAVCODEC_VERSION_MAJOR < 53)
        if (avcodec_open(codecContext, codec) < 0) {
#else
        if (avcodec_open2(codecContext, codec, 0) < 0) {
#endif
            LERRORC(loggerCat_, "Could not load Codec " << codec->id);
            return false;
        }

        return true;
        /* codecContext == videoStream->codec */
    }

    bool checkParams(const int fps, const int width, const int height,
            GLenum pixelFormat, GLenum pixelType, const int bitRate) {

        // TODO check limits of swscale conversion
        if (GL_RGBA != pixelFormat || (GL_FLOAT != pixelType && GL_UNSIGNED_SHORT != pixelType && GL_UNSIGNED_BYTE != pixelType)) {
            LERRORC(loggerCat_,"only GL_RGBA + GL_FLOAT/GL_UNSIGNED_SHORT/GL_UNSIGNED_BYTE as input supported");
            return false;
        }

        if (width <= 0 || height <= 0) {
            LERRORC(loggerCat_, "illegal dimensions: (" << width << "x" << height << ")");
            return false;
        }

        if ((width % 2) + (height % 2) != 0) {
            LERRORC(loggerCat_, "width and height must be multiple of 2 - " << width << ":" << height);
            return false;
        }

        if (fps < minimalFps) {
            LERRORC(loggerCat_, "low fps: " << fps << "; minimum: " << minimalFps);
            return false;
        }

        if (bitRate < minimalBitrate) {
            LERRORC(loggerCat_, "low bitrate: " << bitRate << "; minimum: " << minimalBitrate);
            return false;
        }

        return true;
    }

    std::string initFormatContext(VideoEncoder::ContainerCodecPair preset, std::string filePath) {

        // init format-context
        formatContext_ = avformat_alloc_context();
#if (LIBAVFORMAT_VERSION_MAJOR < 54)
        formatContext_->audio_codec_id = CODEC_ID_NONE;
#else
        formatContext_->audio_codec_id = AV_CODEC_ID_NONE;
#endif
        // let ffmpeg guess what stream-format matches to container-format
        AVOutputFormat* outputFormat;

        if (preset == VideoEncoder::GUESS) {
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
            outputFormat = guess_format(0, filePath.c_str(), 0);
#else
            outputFormat = av_guess_format(0, filePath.c_str(), 0);
#endif

            if (!outputFormat) {
                LWARNINGC(loggerCat_, "Could not guess output format, using MPEG");
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
                outputFormat = guess_format(mpegAppendix, 0, 0);
#else
                outputFormat = av_guess_format(mpegAppendix, 0, 0);
#endif
                filePath.append(".");
                filePath.append(mpegAppendix);
            }

            if (!outputFormat) {
                LERRORC(loggerCat_, "Could not find suitable output format");
                return "";
            }
            formatContext_->oformat = outputFormat;
            formatContext_->video_codec_id = outputFormat->video_codec;
        } else {
            const char* appendix = containerAppendix(preset);
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
            outputFormat = guess_format(appendix, 0, 0);
#else
            outputFormat = av_guess_format(appendix, 0, 0);
#endif
            formatContext_->oformat = outputFormat;
            if (filePath.substr(filePath.size() - strlen(appendix) - 1)
                    == std::string(appendix)) {
                filePath.append(".");
                filePath.append(appendix);
            }
            formatContext_->video_codec_id = containerCodecId(preset);
        }
        return filePath;
    }

    VEAVContext(std::string filePath, const int fps, const int width,
                const int height, GLenum pixelFormat, GLenum pixelType,
                const int bitRate, VideoEncoder::ContainerCodecPair preset) {
        // loggerCat_
        // std::string-s concat with "+", but iostream-string-s concat with "<<"
        loggerCat_ = "";
        loggerCat_ += "tgt.VideoEncoder.AVContext::";
        loggerCat_ += filePath;

#ifndef TGT_DEBUG
        //av_log_set_level(AV_LOG_QUIET);
#endif

        if (!checkParams(fps, width, height, pixelFormat, pixelType, bitRate)) {
            throw tgt::Exception("invalid parameters, see log for details");
        }

        //
        swsContext = 0;
        frameCount = 0;

        pixelFormat_ = pixelFormat;
        pixelType_ = pixelType;

        // avcodec init
        av_register_all();

        filePath = initFormatContext(preset, filePath);
        if (filePath.size() == 0) {
            LERRORC(loggerCat_,"could not determine formatContext");
            return;
        }

        // init video-stream
#if (LIBAVFORMAT_VERSION_MAJOR < 54)
        if (formatContext_->video_codec_id != CODEC_ID_NONE) {
#else
        if (formatContext_->video_codec_id != AV_CODEC_ID_NONE) {
#endif
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
            videoStream = av_new_stream(formatContext_, 1);
#else
            videoStream = avformat_new_stream(formatContext_, NULL);
            videoStream->id = 1;
#endif

            if (!videoStream) {
                throw tgt::Exception("Could not allocate stream");
            }
            if (!(initCodec(preset, bitRate, width, height, fps) && initBuffers())) {
                throw tgt::Exception("Could not init codec, see log for details");
            }
        } else {
            throw tgt::Exception("Failed to allocate format-context");
            return;
        }

#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        // must set parameters
        if (av_set_parameters(formatContext_, 0) < 0) {
            LERRORC(loggerCat_, "invalid format parameters");
            return;
        }
#endif

#ifdef TGT_DEBUG
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        dump_format(formatContext_, 0, filePath.c_str(), 1);//outputs format specs
#else
        av_dump_format(formatContext_, 0, filePath.c_str(), 1);//outputs format specs
#endif
#endif

        // open video's file (will overwrite)
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        if (url_fopen(&formatContext_->pb, filePath.c_str(), URL_WRONLY) < 0) {
#else
        if (avio_open(&formatContext_->pb, filePath.c_str(), AVIO_FLAG_WRITE) < 0) {
#endif
            throw tgt::Exception("file " + filePath + " is not accessible for write-access");
        }

#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        /*
         * must be called before any I/O
         * double videoOutBufferSize for slow fs
         */
        url_setbufsize(formatContext_->pb, 2 * videoOutBufferSize);
#endif

#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        url_resetbuf(formatContext_->pb, URL_WRONLY);
#endif

        // write format's header (if any)
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        av_write_header(formatContext_);
#else
        avformat_write_header(formatContext_, 0);
#endif
    }

    /**
     * assume frames are passed sequential
     */
    void nextFrame(uint8_t* pixels) {
        AVCodecContext* codecContext = videoStream->codec;

        // since pixels are filled inside, this is rather a lightweight allocation
#if LIBAVCODEC_VERSION_INT >= AV_VERSION_INT(55, 45, 101)
        AVFrame* inputFrame = av_frame_alloc();
#else
        AVFrame* inputFrame = avcodec_alloc_frame();
#endif

        /*
         * fill input into inputFrame
         * inputFrame->data ^= pixels is freed outside this scope
         */
        avpicture_fill((AVPicture*) inputFrame, pixels, AV_PIX_FMT_RGB24,
                codecContext->width, codecContext->height);

        // convert RGB24 -> YUV420P
        if (swsContext == 0) {
            swsContext = sws_getContext(codecContext->width,
                    codecContext->height, AV_PIX_FMT_RGB24, codecContext->width,
                    codecContext->height, codecContext->pix_fmt, SWS_BICUBIC,
                    0, 0, 0);
        }
        sws_scale(swsContext, inputFrame->data, inputFrame->linesize, 0,
                codecContext->height, videoFrame->data, videoFrame->linesize);

        av_free(inputFrame);

#if (LIBAVFORMAT_VERSION_MAJOR < 56)
        int outputSize = avcodec_encode_video(codecContext, videoBuffer,
                videoOutBufferSize, videoFrame);
        if (outputSize > 0) {
            AVPacket packet;
            av_init_packet(&packet);

            packet.pts = av_rescale_q(codecContext->coded_frame->pts,
                    codecContext->time_base, videoStream->time_base);

            if (codecContext->coded_frame->key_frame) {
#if (LIBAVFORMAT_VERSION_MAJOR < 53)
                packet.flags |= PKT_FLAG_KEY;
#else
                packet.flags |= AV_PKT_FLAG_KEY;
#endif

            }
            packet.stream_index = videoStream->index;
            packet.data = videoBuffer;
            packet.size = outputSize;

            if (av_write_frame(formatContext_, &packet) != 0) {
                LWARNINGC(loggerCat_, "error while writing frame to file");
            }

            // this would free videoBuffer, we want to reuse
            //av_destruct_packet(&packet);
        }
#else
    // TODO: in newer versions, avcodec_encode_video has been replaced by avcodec_encode_video2, which unfortunately does not work the same way. Since the videoBuffer is not used here, it should be completely removed from the code before, since it is a waste of memory
    AVPacket pkt;
    av_init_packet(&pkt);
    int got_packet;
    pkt.data = videoBuffer;
    pkt.size = videoOutBufferSize;

        int outputSize = avcodec_encode_video2(codecContext, &pkt, videoFrame, &got_packet);

        if (got_packet) {
            pkt.pts = av_rescale_q(pkt.pts, codecContext->time_base, videoStream->time_base);

            pkt.stream_index = videoStream->index;
            //packet.data = videoBuffer;
            //pkt.size = outputSize;

            if (av_write_frame(formatContext_, &pkt) != 0) {
                LWARNINGC(loggerCat_, "error while writing frame to file");
            }
        }

        //av_free_packet(&pkt);
#endif

        frameCount++;
    }

    /**
     * write file trailer and finish encoding
     */
    virtual ~VEAVContext() {

        av_write_trailer(formatContext_);

        avcodec_close(videoStream->codec);

        av_free(videoFrame->data[0]);
        av_free(videoFrame);

        av_free(videoBuffer);

#if (LIBAVFORMAT_VERSION_MAJOR < 53)
        url_fclose(formatContext_->pb);
#else
        avio_close(formatContext_->pb);
#endif

        av_free(formatContext_);
    }
};

VEAVContext* encoderContext;// singleton

} /* anonymous namespace */

std::string VideoEncoder::loggerCat_("tgt.VideoEncoder");

VideoEncoder::VideoEncoder() {
    encoderContext = 0;
    bitrate_ = 2000000;
    preset_ = GUESS;
}

VideoEncoder::~VideoEncoder() {
}

void VideoEncoder::startVideoEncoding(std::string filePath, const int fps,
                                      const int width, const int height,
                                      GLint pixelFormat, GLenum pixelType) {
    if (encoderContext != 0) {
        LWARNING("already encoding a video");
        return;
    }

    // we must not reuse this context because it still stores state of finished encoding
    // TODO consider reusing everything but videoStream
    encoderContext = new VEAVContext(filePath, fps, width, height, pixelFormat,
            pixelType, bitrate_, preset_);
}

void VideoEncoder::stopVideoEncoding() {
    if (encoderContext == 0) {
        LERROR("Encoding has already aborted");
        return;
    }

    AVCodecContext* codecContext = encoderContext->videoStream->codec;

    if (encoderContext->frameCount % codecContext->time_base.den != 0) {
        LWARNING("passed frames doesn't match fps: " << encoderContext->frameCount << "/" << codecContext->time_base.den);
    }

    delete encoderContext;
    encoderContext = 0;
}

/**
 * TODO assumption of float not applied here cause we could specify pixelComponents' type when starting animation
 */
void VideoEncoder::nextFrame(GLvoid* pixels) {
    if (pixels == 0)
        return;

    AVCodecContext* codecContext = encoderContext->videoStream->codec;

    if (encoderContext->pixelFormat_ != GL_RGBA) {
        LERROR("Only GL_RGBA as input format supported.");
        return;
    }

    uint8_t* pixels8 = 0;
    if (encoderContext->pixelType_ == GL_FLOAT) {
        pixels8 = floatToUint8_tRGB(reinterpret_cast<float*>(pixels), codecContext->width
                    * codecContext->height, codecContext->width, FLIP_VERT);
    }
    else if (encoderContext->pixelType_ == GL_UNSIGNED_SHORT) {
        pixels8 = uint16ToUint8_tRGB(reinterpret_cast<uint16_t*>(pixels), codecContext->width
                    * codecContext->height, codecContext->width, FLIP_VERT);
    }
    else if (encoderContext->pixelType_ == GL_UNSIGNED_BYTE) {
        pixels8 = uint8ToUint8_tRGB(reinterpret_cast<uint8_t*>(pixels), codecContext->width
                    * codecContext->height, codecContext->width, FLIP_VERT);
    }
    else {
        LERROR("Unknown pixel type. Only GL_FLOAT/GL_UNSIGNED_SHORT/GL_UNSIGNED_BYTE supported");
        return;
    }

    encoderContext->nextFrame(pixels8);
    delete [] pixels8;
}

void VideoEncoder::setup(int preset, int bitrate) {
    preset_ = (ContainerCodecPair) preset;
    bitrate_ = bitrate;
}

std::vector<std::string> VideoEncoder::getSupportedFormatsByFileEnding() {
    std::vector<std::string> formats;
    formats.push_back("avi");
    formats.push_back("mpg");
    formats.push_back("wmv");
    formats.push_back("flv");
    formats.push_back("ogg");

    return formats;
}

std::vector<std::string> VideoEncoder::getSupportedFormatDescriptions() {
    std::vector<std::string> formats;
    formats.push_back("AVI video files");
    formats.push_back("MPEG video files");
    formats.push_back("WMV video files");
    formats.push_back("Flash FLV video files");
    formats.push_back("OGG video files");

    return formats;
}

const char** VideoEncoder::getContainerCodecPairNames() {
    return containerCodecPairNames;
}

/**
 * static
 */
const char* VideoEncoder::getContainerAppendix() {
    return containerAppendix(preset_);
}

} // namespace tgt
