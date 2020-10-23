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

#include "analyzevolumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/quaternion.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

using tgt::vec3;
using tgt::mat4;
using tgt::ivec3;
using tgt::hor;
using tgt::lessThanEqual;

namespace voreen {

const std::string AnalyzeVolumeReader::loggerCat_ = "voreen.io.VolumeReader.analyze";

void swap_long(unsigned char *pntr) {
    unsigned char b0, b1, b2, b3;
    b0 = *pntr;
    b1 = *(pntr+1);
    b2 = *(pntr+2);
    b3 = *(pntr+3);
    *pntr = b3;
    *(pntr+1) = b2;
    *(pntr+2) = b1;
    *(pntr+3) = b0;
}

void swap_short(unsigned char *pntr) {
    unsigned char b0, b1;
    b0 = *pntr;
    b1 = *(pntr+1);
    *pntr = b1;
    *(pntr+1) = b0;
}

struct header_key {
    int32_t sizeof_hdr;
    char data_type[10];
    char db_name[18];
    int32_t extents;
    int16_t session_error;
    char regular;
    char hkey_un0;

    void swapEndianess() {
        swap_long(reinterpret_cast<unsigned char*>(&sizeof_hdr));
        swap_long(reinterpret_cast<unsigned char*>(&extents));
        swap_short(reinterpret_cast<unsigned char*>(&session_error));
    }
};
/* 40 bytes */

struct image_dimension
{
    int16_t dim[8];
    int16_t unused8;
    int16_t unused9;
    int16_t unused10;
    int16_t unused11;
    int16_t unused12;
    int16_t unused13;
    int16_t unused14;
    int16_t datatype;
    int16_t bitpix;
    int16_t dim_un0;
    float pixdim[8];
    float vox_offset;
    float funused1;
    float funused2;
    float funused3;
    float cal_max;
    float cal_min;
    float compressed;
    float verified;
    int32_t glmax,glmin;

    void swapEndianess() {
        swap_short(reinterpret_cast<unsigned char*>(&dim[0]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[1]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[2]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[3]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[4]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[5]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[6]));
        swap_short(reinterpret_cast<unsigned char*>(&dim[7]));
        swap_short(reinterpret_cast<unsigned char*>(&unused8));
        swap_short(reinterpret_cast<unsigned char*>(&unused9));
        swap_short(reinterpret_cast<unsigned char*>(&unused10));
        swap_short(reinterpret_cast<unsigned char*>(&unused11));
        swap_short(reinterpret_cast<unsigned char*>(&unused12));
        swap_short(reinterpret_cast<unsigned char*>(&unused13));
        swap_short(reinterpret_cast<unsigned char*>(&unused14));
        swap_short(reinterpret_cast<unsigned char*>(&datatype));
        swap_short(reinterpret_cast<unsigned char*>(&bitpix));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[0]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[1]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[2]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[3]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[4]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[5]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[5]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[6]));
        swap_long(reinterpret_cast<unsigned char*>(&pixdim[7]));
        swap_long(reinterpret_cast<unsigned char*>(&vox_offset));
        swap_long(reinterpret_cast<unsigned char*>(&funused1));
        swap_long(reinterpret_cast<unsigned char*>(&funused2));
        swap_long(reinterpret_cast<unsigned char*>(&cal_max));
        swap_long(reinterpret_cast<unsigned char*>(&cal_min));
        swap_long(reinterpret_cast<unsigned char*>(&compressed));
        swap_long(reinterpret_cast<unsigned char*>(&verified));
        swap_short(reinterpret_cast<unsigned char*>(&dim_un0));
        swap_long(reinterpret_cast<unsigned char*>(&glmax));
        swap_long(reinterpret_cast<unsigned char*>(&glmin));
    }
};
/* 108 bytes */

struct data_history
{
    char descrip[80];
    char aux_file[24];
    char orient;
    char originator[10];
    char generated[10];
    char scannum[10];
    char patient_id[10];
    char exp_date[10];
    char exp_time[10];
    char hist_un0[3];
    int32_t views;
    int32_t vols_added;
    int32_t start_field;
    int32_t field_skip;
    int32_t omax, omin;
    int32_t smax, smin;

    void swapEndianess() {
        swap_long(reinterpret_cast<unsigned char*>(&views));
        swap_long(reinterpret_cast<unsigned char*>(&vols_added));
        swap_long(reinterpret_cast<unsigned char*>(&start_field));
        swap_long(reinterpret_cast<unsigned char*>(&field_skip));
        swap_long(reinterpret_cast<unsigned char*>(&omax));
        swap_long(reinterpret_cast<unsigned char*>(&omin));
        swap_long(reinterpret_cast<unsigned char*>(&smax));
        swap_long(reinterpret_cast<unsigned char*>(&smin));
    }
};
/* 200 bytes */

struct analyze_header
{
    header_key hk;
    image_dimension dime;
    data_history hist;
};
/* 348 bytes */

/* Acceptable values for datatype */
#define VRN_ANALYZE_DT_NONE 0
#define VRN_ANALYZE_DT_UNKNOWN 0
#define VRN_ANALYZE_DT_BINARY 1
#define VRN_ANALYZE_DT_UNSIGNED_CHAR 2
#define VRN_ANALYZE_DT_SIGNED_SHORT 4
#define VRN_ANALYZE_DT_SIGNED_INT 8
#define VRN_ANALYZE_DT_FLOAT 16
#define VRN_ANALYZE_DT_COMPLEX 32
#define VRN_ANALYZE_DT_DOUBLE 64
#define VRN_ANALYZE_DT_RGB 128
#define VRN_ANALYZE_DT_ALL 255

                            /*------------------- new codes for NIFTI ---*/
#define VRN_ANALYZE_DT_INT8                  256     /* signed char (8 bits)         */
#define VRN_ANALYZE_DT_UINT16                512     /* unsigned short (16 bits)     */
#define VRN_ANALYZE_DT_UINT32                768     /* unsigned int (32 bits)       */
#define VRN_ANALYZE_DT_INT64                1024     /* long long (64 bits)          */
#define VRN_ANALYZE_DT_UINT64               1280     /* unsigned long long (64 bits) */
#define VRN_ANALYZE_DT_FLOAT128             1536     /* long double (128 bits)       */
#define VRN_ANALYZE_DT_COMPLEX128           1792     /* double pair (128 bits)       */
#define VRN_ANALYZE_DT_COMPLEX256           2048     /* long double pair (256 bits)  */
#define VRN_ANALYZE_DT_RGBA32               2304     /* 4 byte RGBA (32 bits/voxel)  */

/* Acceptable values for intent codes, see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1diagrams_v2.pdf */
#define IC_INTENT_NONE          0
#define IC_INTENT_ESTIMATE      1001    /* parameter for estimate in intent_name */
#define IC_INTENT_LABEL         1002    /* parameter at each voxel is index to label defined in aux_file */
#define IC_INTENT_NEURONAME     1003    /* parameter at each voxel is index to label in NeuroNames label set */
#define IC_INTENT_GENMATRIX     1004    /* parameter at each voxel is matrix */
#define IC_INTENT_SYMMATRIX     1005    /* parameter at each voxel is symmetrical matrix */
#define IC_INTENT_DISPVECT      1006    /* parameter at each voxel is displacement vector */
#define IC_INTENT_VECTOR        1007    /* parameter at each voxel is vector */
#define IC_INTENT_POINTSET      1008    /* value at each voxel is spatial coordinate (vertices/nodes of surface mesh) */
#define IC_INTENT_TRIANGLE      1009    /* value at each voxel is spatial coordinate (vertices/nodes of surface mesh) */
#define IC_INTENT_QUATERNION    1010
#define IC_INTENT_DIMLESS       1011    /* dimensionless value */

#define NIFTI_UNITS_UNKNOWN 0

/** Space codes are multiples of 1. **/
#define NIFTI_UNITS_METER   1 /*! NIFTI code for meters. */
#define NIFTI_UNITS_MM      2 /*! NIFTI code for millimeters. */
#define NIFTI_UNITS_MICRON  3 /*! NIFTI code for micrometers. */

/** Time codes are multiples of 8. **/
#define NIFTI_UNITS_SEC     8 /*! NIFTI code for seconds. */
#define NIFTI_UNITS_MSEC   16 /*! NIFTI code for milliseconds. */
#define NIFTI_UNITS_USEC   24 /*! NIFTI code for microseconds. */

/*** These units are for spectral data: ***/
#define NIFTI_UNITS_HZ     32 /*! NIFTI code for Hertz. */
#define NIFTI_UNITS_PPM    40 /*! NIFTI code for ppm. */
#define NIFTI_UNITS_RADS   48 /*! NIFTI code for radians per second. */

#define XYZT_TO_SPACE(xyzt)       ( (xyzt) & 0x07 )
#define XYZT_TO_TIME(xyzt)        ( (xyzt) & 0x38 )

//struct complex
//{
    //float real;
    //float imag;
//};

struct nifti_1_header { /* NIFTI-1 usage         */  /* ANALYZE 7.5 field(s) */
                        /*************************/  /************************/
                                           /*--- was header_key substruct ---*/
 int32_t   sizeof_hdr;    /*!< MUST be 348           */  /* int sizeof_hdr;      */
 char  data_type[10]; /*!< ++UNUSED++            */  /* char data_type[10];  */
 char  db_name[18];   /*!< ++UNUSED++            */  /* char db_name[18];    */
 int32_t   extents;       /*!< ++UNUSED++            */  /* int extents;         */
 int16_t session_error; /*!< ++UNUSED++            */  /* short session_error; */
 char  regular;       /*!< ++UNUSED++            */  /* char regular;        */
 char  dim_info;      /*!< MRI slice ordering.   */  /* char hkey_un0;       */

                                      /*--- was image_dimension substruct ---*/
 int16_t dim[8];        /*!< Data array dimensions.*/  /* short dim[8];        */
 float intent_p1 ;    /*!< 1st intent parameter. */  /* int16_t unused8;       */
                                                     /* int16_t unused9;       */
 float intent_p2 ;    /*!< 2nd intent parameter. */  /* int16_t unused10;      */
                                                     /* int16_t unused11;      */
 float intent_p3 ;    /*!< 3rd intent parameter. */  /* int16_t unused12;      */
                                                     /* int16_t unused13;      */
 int16_t intent_code ;  /*!< NIFTI_INTENT_* code.  */  /* short unused14;      */
 int16_t datatype;      /*!< Defines data type!    */  /* short datatype;      */
 int16_t bitpix;        /*!< Number bits/voxel.    */  /* short bitpix;        */
 int16_t slice_start;   /*!< First slice index.    */  /* short dim_un0;       */
 float pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];     */
 float vox_offset;    /*!< Offset into .nii file */  /* float vox_offset;    */
 float scl_slope ;    /*!< Data scaling: slope.  */  /* float funused1;      */
 float scl_inter ;    /*!< Data scaling: offset. */  /* float funused2;      */
 int16_t slice_end;     /*!< Last slice index.     */  /* float funused3;      */
 char  slice_code ;   /*!< Slice timing order.   */
 char  xyzt_units ;   /*!< Units of pixdim[1..4] */
 float cal_max;       /*!< Max display intensity */  /* float cal_max;       */
 float cal_min;       /*!< Min display intensity */  /* float cal_min;       */
 float slice_duration;/*!< Time for 1 slice.     */  /* float compressed;    */
 float toffset;       /*!< Time axis shift.      */  /* float verified;      */
 int32_t   glmax;         /*!< ++UNUSED++            */  /* int glmax;           */
 int32_t   glmin;         /*!< ++UNUSED++            */  /* int glmin;           */

                                         /*--- was data_history substruct ---*/
 char  descrip[80];   /*!< any text you like.    */  /* char descrip[80];    */
 char  aux_file[24];  /*!< auxiliary filename.   */  /* char aux_file[24];   */

 int16_t qform_code ;   /*!< NIFTI_XFORM_* code.   */  /*-- all ANALYZE 7.5 ---*/
 int16_t sform_code ;   /*!< NIFTI_XFORM_* code.   */  /*   fields below here  */
                                                     /*   are replaced       */
 float quatern_b ;    /*!< Quaternion b param.   */
 float quatern_c ;    /*!< Quaternion c param.   */
 float quatern_d ;    /*!< Quaternion d param.   */
 float qoffset_x ;    /*!< Quaternion x shift.   */
 float qoffset_y ;    /*!< Quaternion y shift.   */
 float qoffset_z ;    /*!< Quaternion z shift.   */

 float srow_x[4] ;    /*!< 1st row affine transform.   */
 float srow_y[4] ;    /*!< 2nd row affine transform.   */
 float srow_z[4] ;    /*!< 3rd row affine transform.   */

 char intent_name[16];/*!< 'name' or meaning of data.  */

 char magic[4] ;      /*!< MUST be "ni1\0" or "n+1\0". */

 void swapEndianess() {
     swap_long(reinterpret_cast<unsigned char*>(&sizeof_hdr));
     swap_long(reinterpret_cast<unsigned char*>(&extents));
     swap_short(reinterpret_cast<unsigned char*>(&session_error));
     swap_short(reinterpret_cast<unsigned char*>(&dim[0]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[1]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[2]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[3]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[4]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[5]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[6]));
     swap_short(reinterpret_cast<unsigned char*>(&dim[7]));
     swap_long(reinterpret_cast<unsigned char*>(&intent_p1));
     swap_long(reinterpret_cast<unsigned char*>(&intent_p2));
     swap_long(reinterpret_cast<unsigned char*>(&intent_p3));
     swap_short(reinterpret_cast<unsigned char*>(&intent_code));
     swap_short(reinterpret_cast<unsigned char*>(&datatype));
     swap_short(reinterpret_cast<unsigned char*>(&bitpix));
     swap_short(reinterpret_cast<unsigned char*>(&slice_start));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[0]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[1]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[2]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[3]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[4]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[5]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[6]));
     swap_long(reinterpret_cast<unsigned char*>(&pixdim[7]));
     swap_long(reinterpret_cast<unsigned char*>(&vox_offset));
     swap_long(reinterpret_cast<unsigned char*>(&scl_slope));
     swap_long(reinterpret_cast<unsigned char*>(&scl_inter));
     swap_short(reinterpret_cast<unsigned char*>(&slice_end));
     swap_long(reinterpret_cast<unsigned char*>(&cal_max));
     swap_long(reinterpret_cast<unsigned char*>(&cal_min));
     swap_long(reinterpret_cast<unsigned char*>(&slice_duration));
     swap_long(reinterpret_cast<unsigned char*>(&toffset));
     swap_long(reinterpret_cast<unsigned char*>(&glmax));
     swap_long(reinterpret_cast<unsigned char*>(&glmin));
     swap_short(reinterpret_cast<unsigned char*>(&qform_code));
     swap_short(reinterpret_cast<unsigned char*>(&sform_code));
     swap_long(reinterpret_cast<unsigned char*>(&quatern_b));
     swap_long(reinterpret_cast<unsigned char*>(&quatern_c));
     swap_long(reinterpret_cast<unsigned char*>(&quatern_d));
     swap_long(reinterpret_cast<unsigned char*>(&qoffset_x));
     swap_long(reinterpret_cast<unsigned char*>(&qoffset_y));
     swap_long(reinterpret_cast<unsigned char*>(&qoffset_z));

     swap_long(reinterpret_cast<unsigned char*>(&srow_x[0]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_x[1]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_x[2]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_x[3]));

     swap_long(reinterpret_cast<unsigned char*>(&srow_y[0]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_y[1]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_y[2]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_y[3]));

     swap_long(reinterpret_cast<unsigned char*>(&srow_z[0]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_z[1]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_z[2]));
     swap_long(reinterpret_cast<unsigned char*>(&srow_z[3]));
 }
};                   /**** 348 bytes total ****/

AnalyzeVolumeReader::AnalyzeVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("hdr");
    extensions_.push_back("nii");

    //make sure the header sizes are ok when compiled on all architectures:
    tgtAssert(sizeof(header_key) == 40, "Wrong header_key size!");
    tgtAssert(sizeof(image_dimension) == 108, "Wrong image_dimension size!");
    tgtAssert(sizeof(data_history) == 200, "Wrong data_history size!");
    tgtAssert(sizeof(analyze_header) == 348, "Wrong analyze_header size!");
    tgtAssert(sizeof(nifti_1_header ) == 348, "Wrong nifti_1_header size!");
}

std::string AnalyzeVolumeReader::getRelatedImgFileName(const std::string& fileName) {
    std::string imgFileName = fileName.substr(0, fileName.length()-3);
    imgFileName += "img";
    return imgFileName;
}

VolumeBase* AnalyzeVolumeReader::read(const VolumeURL& origin) {
    VolumeBase* result = 0;

    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (!tmp.empty())
        volumeId = stoi(tmp);

    int channel = 0;
    tmp = origin.getSearchParameter("channel");
    if(!tmp.empty())
        channel = stoi(tmp);

    std::unique_ptr<VolumeList> collection(read(origin.getPath(), volumeId, channel));

    if (collection && collection->size() == 1) {
        result = collection->first();
    }
    else if (collection && collection->size() > 1) {
        while(!collection->empty()) {
           VolumeBase* vh = collection->first();
           collection->remove(vh);
           delete vh;
        }
        throw tgt::FileException("Only one volume expected", origin.getPath());
    }

    return result;
}

VolumeList* AnalyzeVolumeReader::read(const std::string &url) {
    return read(url, -1, -1);
}

VolumeList* AnalyzeVolumeReader::read(const std::string &url, int volId, int channel) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    //check if we got analyze or nifti:

    //extension .nii => standalone nifti
    if(fileName.find(".nii") != std::string::npos)
        return readNifti(fileName, true, volId, channel);
    else {
        //check magic string:
        std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary);
        if(!file) {
            throw tgt::FileNotFoundException("Failed to open file: ", fileName);
        }

        file.seekg(0, std::ios::end);
        std::streamoff fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        if(fileSize < 348)
            throw tgt::CorruptedFileException("Analyze/Nifti file is smaller than 348 bytes!", fileName);

        nifti_1_header header;
        if (!file.read((char*)&header, sizeof(header))) {
            throw tgt::CorruptedFileException("Failed to read header!", fileName);
        }
        file.close();

        if( (header.magic[0] == 'n') && (header.magic[2] == '1') && (header.magic[3] == 0) ) {
            if(header.magic[1] == '+') {
                return readNifti(fileName, true, volId, channel);
            }
            else if(header.magic[1] == 'i') {
                return readNifti(fileName, false, volId, channel);
            }
            else
                return readAnalyze(fileName, volId, channel);
        }
        else
            return readAnalyze(fileName, volId, channel);
    }
    return 0;
}

VolumeList* AnalyzeVolumeReader::readNifti(const std::string &fileName, bool standalone, int volId, int channel) {
    LINFO("Loading nifti file " << fileName);

    std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary);
    if(!file) {
        throw tgt::FileNotFoundException("Failed to open file: ", fileName);
    }

    nifti_1_header header;
    if (!file.read((char*)&header, sizeof(header))) {
        throw tgt::CorruptedFileException("Failed to read header!", fileName);
    }

    file.close();

    bool bigEndian = false;
    //check if swap is necessary:
    if((header.dim[0] < 0) || (header.dim[0] > 15)) {
        bigEndian = true;
        header.swapEndianess();
    }

    if(header.sizeof_hdr != 348) {
        throw tgt::CorruptedFileException("Invalid header.sizeof_hdr", fileName);
    }

    if(!( (header.magic[0] == 'n') && (header.magic[2] == '1') && (header.magic[3] == 0) ))
        throw tgt::CorruptedFileException("Not a Nifti header!", fileName);

    if(header.magic[1] == '+') {
        if(!standalone)
            LWARNING("Tried to read standalone Nifti as hdr+img!");
        standalone = true;
    }
    else if(header.magic[1] == 'i') {
        if(!standalone)
            LWARNING("Tried to read hdr+img Nifti as standalone!");
        standalone = false;
    }
    else
        throw tgt::CorruptedFileException("Not a Nifti header!", fileName);

    ivec3 dimensions;
    dimensions.x = header.dim[1];
    dimensions.y = header.dim[2];
    dimensions.z = header.dim[3];
    LINFO("Resolution: " << dimensions);

    int numVolumes = header.dim[4];
    LINFO("Number of volumes: " << numVolumes);

    int numChannels = header.dim[0] < 5 ? 1 : header.dim[5];
    LINFO("Number of channels: " << numChannels);

    if (hor(lessThanEqual(dimensions, tgt::ivec3::zero))) {
        LERROR("Invalid resolution or resolution not specified: " << dimensions);
        throw tgt::CorruptedFileException("error while reading data", fileName);
    }

    vec3 spacing;
    spacing.x = header.pixdim[1];
    spacing.y = header.pixdim[2];
    spacing.z = header.pixdim[3];
    LINFO("Spacing: " << spacing);

    int timeunit = XYZT_TO_TIME(header.xyzt_units);
    int spaceunit = XYZT_TO_SPACE(header.xyzt_units);
    LINFO("timeunit: " << timeunit << " spaceunit: " << spaceunit);

    float dt = header.pixdim[4];
    float toffset = header.toffset;
    switch(timeunit) {
        case NIFTI_UNITS_SEC:
            dt *= 1000.0f;
            toffset *= 1000.0f;
            break;
        case NIFTI_UNITS_MSEC:
            //nothing to do
            break;
        case NIFTI_UNITS_USEC:
            dt /= 1000.0f;
            toffset /= 1000.0f;
            break;
    }

    switch(spaceunit) {
        case NIFTI_UNITS_MM:
            //nothing to do
            break;
        case NIFTI_UNITS_METER:
            spacing *= 1000.0f;
            LWARNING("Units: meter");
            break;
        case NIFTI_UNITS_MICRON:
            spacing /= 1000.0f;
            LWARNING("Units: micron");
            break;
        case NIFTI_UNITS_UNKNOWN:
        default:
            LWARNING("Unknown space unit!");
            break;
    }

    LINFO("Datatype: " << header.datatype);
    std::string baseType = "";
    RealWorldMapping denormalize;
    bool applyRWM = header.scl_slope != 0.0f;

    switch(header.intent_code) {
        case IC_INTENT_SYMMATRIX:  /* parameter at each voxel is symmetrical matrix */
            //TODO: should be relatively easy (=> tensors)
        case IC_INTENT_DISPVECT:   /* parameter at each voxel is displacement vector */
        case IC_INTENT_VECTOR:     /* parameter at each voxel is vector */
            //TODO: should be relatively easy
        case IC_INTENT_GENMATRIX:  /* parameter at each voxel is matrix */
            //TODO: should be relatively easy
        case IC_INTENT_POINTSET:   /* value at each voxel is spatial coordinate (vertices/nodes of surface mesh) */
        case IC_INTENT_TRIANGLE:   /* value at each voxel is spatial coordinate (vertices/nodes of surface mesh) */
        case IC_INTENT_QUATERNION:
            throw tgt::UnsupportedFormatException("Unsupported intent code!");
            break;
        case IC_INTENT_ESTIMATE:   /* parameter for estimate in intent_name */
        case IC_INTENT_LABEL:      /* parameter at each voxel is index to label defined in aux_file */
        case IC_INTENT_NEURONAME:  /* parameter at each voxel is index to label in NeuroNames label set */
        case IC_INTENT_DIMLESS:    /* dimensionless value */
        case IC_INTENT_NONE:
            break;
        default:
            LWARNING("Unhandled intent code");
            break;
    }
    //if (header.intent_code == IC_INTENT_SYMMATRIX) {
        //h.objectModel_ = "TENSOR_FUSION_LOW";
    //}
    if(baseType.empty()) {
        switch(header.datatype) {
            case VRN_ANALYZE_DT_UNSIGNED_CHAR:
                baseType = "uint8";
                denormalize = RealWorldMapping::createDenormalizingMapping<uint8_t>();
                break;
            case VRN_ANALYZE_DT_SIGNED_SHORT:
                baseType = "int16";
                denormalize = RealWorldMapping::createDenormalizingMapping<int16_t>();
                break;
            case VRN_ANALYZE_DT_SIGNED_INT:
                baseType = "int32";
                denormalize = RealWorldMapping::createDenormalizingMapping<int32_t>();
                break;
            case VRN_ANALYZE_DT_FLOAT:
                baseType = "float";
                break;
            case VRN_ANALYZE_DT_DOUBLE:
                baseType = "double";
                break;
            case VRN_ANALYZE_DT_RGB:
                baseType = "Vector3(uint8)";
                applyRWM = false;
                break;
            case VRN_ANALYZE_DT_RGBA32:         /* 4 byte RGBA (32 bits/voxel)  */
                baseType = "Vector4(uint8)";
                applyRWM = false;
                break;
            case VRN_ANALYZE_DT_INT8:           /* signed char (8 bits)         */
                baseType = "int8";
                break;
            case VRN_ANALYZE_DT_UINT16:         /* unsigned short (16 bits)     */
                baseType = "uint16";
                denormalize = RealWorldMapping::createDenormalizingMapping<uint16_t>();
                break;
            case VRN_ANALYZE_DT_UINT32:         /* unsigned int (32 bits)       */
                baseType = "uint32";
                denormalize = RealWorldMapping::createDenormalizingMapping<uint32_t>();
                break;
            case VRN_ANALYZE_DT_INT64:          /* long long (64 bits)          */
                baseType = "int64";
                break;
            case VRN_ANALYZE_DT_UINT64:         /* unsigned long long (64 bits) */
                baseType = "uint64";
                denormalize = RealWorldMapping::createDenormalizingMapping<uint64_t>();
                break;
            case VRN_ANALYZE_DT_FLOAT128:       /* long double (128 bits)       */
            case VRN_ANALYZE_DT_COMPLEX128:     /* double pair (128 bits)       */
            case VRN_ANALYZE_DT_COMPLEX256:     /* long double pair (256 bits)  */
            case VRN_ANALYZE_DT_ALL:
            case VRN_ANALYZE_DT_COMPLEX:
            case 0: //VRN_ANALYZE_DT_NONE/VRN_ANALYZE_DT_UNKNOWN
            case VRN_ANALYZE_DT_BINARY:
            default:
                throw tgt::UnsupportedFormatException("Unsupported datatype!");
        }
    }

    RealWorldMapping rwm(header.scl_slope, header.scl_inter, "");

    int headerskip = static_cast<uint16_t>(header.vox_offset);

    std::string rawFilename = fileName;
    if(!standalone)
        rawFilename = getRelatedImgFileName(fileName);

    mat4 pToW = mat4::identity;

    //Calculate transformation:
    if(header.sform_code > 0) {
        mat4 vToW(header.srow_x[0], header.srow_x[1], header.srow_x[2], header.srow_x[3],
                  header.srow_y[0], header.srow_y[1], header.srow_y[2], header.srow_y[3],
                  header.srow_z[0], header.srow_z[1], header.srow_z[2], header.srow_z[3],
                  0.0f, 0.0f, 0.0f, 1.0f);

        // When exporting, e.g. from matlab, the transformation might be present but invalid.
        if(vToW.getRotationalPartMat3() == tgt::mat3::zero) {
            vToW = mat4::identity;
            LWARNING("Invalid transformation, setting to identity");
        }

        mat4 pToV = mat4::createScale(1.0f/spacing); //no offset
        pToW = vToW * pToV;
    }
    else if(header.qform_code > 0) {
        float b = header.quatern_b;
        float c = header.quatern_c;
        float d = header.quatern_d;
        float a = static_cast<float>(sqrt(1.0-(b*b+c*c+d*d)));

        mat4 rot2(a*a+b*b-c*c-d*d,   2*b*c-2*a*d,       2*b*d+2*a*c,     0.0f,
                  2*b*c+2*a*d,       a*a+c*c-b*b-d*d,   2*c*d-2*a*b,     0.0f,
                  2*b*d-2*a*c,       2*c*d+2*a*b,       a*a+d*d-c*c-b*b, 0.0f,
                  0.0f,              0.0f,              0.0f,            1.0f);

        float qfac = header.pixdim[0];
        if(fabs(qfac) < 0.1f)
            qfac = 1.0f;
        mat4 sc = mat4::createScale(vec3(1.0f, 1.0f, qfac));

        mat4 os = mat4::createTranslation(vec3(header.qoffset_x, header.qoffset_y, header.qoffset_z));
        pToW = os * rot2 * sc;
    }

    VolumeList* vc = new VolumeList();
    size_t volSize = hmul(tgt::svec3(dimensions)) * (header.bitpix / 8);

    // Value -1 means loading all volumes / channels.
    int startId = 0;
    int stopId = numVolumes;
    if(volId != -1) {
        startId = volId;
        stopId = startId + 1;
    }
    int startChannel = 0;
    int stopChannel = numChannels;
    if(channel != -1) {
        startChannel = channel;
        stopChannel = startChannel +1;
    }

    for(int i=startId; i < stopId; i++) {
        for(int channelId=startChannel; channelId<stopChannel; channelId++) {
            size_t offset = headerskip + (channelId * numVolumes + i) * volSize;
            VolumeRepresentation* volume = new VolumeDiskRaw(rawFilename, baseType, dimensions, offset, bigEndian);
            Volume* vh = new Volume(volume, spacing, tgt::vec3::zero);

            VolumeURL origin(fileName);
            origin.addSearchParameter("volumeId", itos(i));
            origin.addSearchParameter("channel", itos(channelId));
            vh->setOrigin(origin);

            vh->setPhysicalToWorldMatrix(pToW);
            vh->setMetaDataValue<StringMetaData>("Description", std::string(header.descrip));
            vh->setMetaDataValue<IntMetaData>("FrameTime", static_cast<int>(toffset + i * dt));
            vh->setTimestep(toffset + i * dt);
            if (applyRWM)
                vh->setRealWorldMapping(RealWorldMapping::combine(denormalize, rwm));

            vc->add(vh);
        }
    }

    return vc;
}

VolumeList* AnalyzeVolumeReader::readAnalyze(const std::string &fileName, int volId, int channel) {
    LINFO("Loading analyze file " << fileName);
    LINFO("Related img file: " << getRelatedImgFileName(fileName));

    std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary);
    if(!file) {
        throw tgt::FileNotFoundException("Failed to open file: ", fileName);
    }

    file.seekg(0, std::ios::end);
    std::streamoff fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    if(fileSize != 348)
        LWARNING("Filesize != 348");

    header_key header;
    if (!file.read((char*)&header, sizeof(header))) {
        throw tgt::CorruptedFileException("Failed to read header!", fileName);
    }

    image_dimension dimension;
    if (!file.read((char*)&dimension, sizeof(dimension))) {
        throw tgt::CorruptedFileException("Failed to read dimensions!", fileName);
    }

    data_history history;
    if (!file.read((char*)&history, sizeof(history))) {
        throw tgt::CorruptedFileException("Failed to read history!", fileName);
    }

    bool bigEndian = false;
    //check if swap is necessary:
    if((dimension.dim[0] < 0) || (dimension.dim[0] > 15)) {
        bigEndian = true;
        header.swapEndianess();
        dimension.swapEndianess();
        history.swapEndianess();
    }

    ivec3 dimensions;
    dimensions.x = dimension.dim[1];
    dimensions.y = dimension.dim[2];
    dimensions.z = dimension.dim[3];
    LINFO("Resolution: " << dimensions);

    int numVolumes = dimension.dim[4];
    LINFO("Number of volumes: " << numVolumes);

    int numChannels = dimension.dim[0] < 5 ? 1 : dimension.dim[5];
    LINFO("Number of channels: " << numChannels);

    if (hor(lessThanEqual(dimensions, ivec3(0)))) {
        LERROR("Invalid resolution or resolution not specified: " << dimensions);
        throw tgt::CorruptedFileException("error while reading data", fileName);
    }

    vec3 spacing;
    spacing.x = dimension.pixdim[1];
    spacing.y = dimension.pixdim[2];
    spacing.z = dimension.pixdim[3];
    LINFO("Spacing: " << spacing);

    LINFO("Datatype: " << dimension.datatype);

    std::string baseType;
    switch(dimension.datatype) {
        case VRN_ANALYZE_DT_UNSIGNED_CHAR:
            baseType = "uint8";
            break;
        case VRN_ANALYZE_DT_SIGNED_SHORT:
            baseType = "int16";
            break;
        case VRN_ANALYZE_DT_SIGNED_INT:
            baseType = "int32";
            break;
        case VRN_ANALYZE_DT_FLOAT:
            baseType = "float";
            break;
        case VRN_ANALYZE_DT_DOUBLE:
            baseType = "double";
            break;
        case VRN_ANALYZE_DT_RGB:
            baseType = "Vector3(uint8)";
            break;
        case VRN_ANALYZE_DT_ALL:
        case VRN_ANALYZE_DT_COMPLEX:
        case 0: //VRN_ANALYZE_DT_NONE/VRN_ANALYZE_DT_UNKNOWN
        case VRN_ANALYZE_DT_BINARY:
        default:
            throw tgt::UnsupportedFormatException("Unsupported datatype!");
    }

    // Value -1 means loading all volumes / channels.
    int startId = 0;
    int stopId = numVolumes;
    if(volId != -1) {
        startId = volId;
        stopId = startId + 1;
    }
    int startChannel = 0;
    int stopChannel = numChannels;
    if(channel != -1) {
        startChannel = channel;
        stopChannel = startChannel +1;
    }

    // Nifti transformations give us the center of the first voxel, we translate to correct:
    mat4 pToW = mat4::createTranslation(-spacing * 0.5f);

    VolumeList* vc = new VolumeList();
    size_t volSize = hmul(tgt::svec3(dimensions)) * (dimension.bitpix / 8);
    for(int i=startId; i<stopId; i++) {
        for(int channelId=startChannel; channelId<stopChannel; channelId++) {
            size_t offset = (channelId * numVolumes + i) * volSize;
            VolumeRepresentation* volume = new VolumeDiskRaw(getRelatedImgFileName(fileName), baseType, dimensions,
                                                             offset, bigEndian);
            Volume* vh = new Volume(volume, spacing, tgt::vec3::zero);
            vh->setPhysicalToWorldMatrix(pToW);

            VolumeURL origin(fileName);
            origin.addSearchParameter("volumeId", itos(i));
            origin.addSearchParameter("channel", itos(channelId));
            vh->setOrigin(origin);

            vc->add(vh);
        }
    }

    return vc;
}

VolumeReader* AnalyzeVolumeReader::create(ProgressBar* progress) const {
    return new AnalyzeVolumeReader(progress);
}

std::vector<VolumeURL> AnalyzeVolumeReader::listVolumes(const std::string& url) const {
    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary);
    if(!file) {
        throw tgt::FileNotFoundException("Failed to open file: ", fileName);
    }

    nifti_1_header header;
    if (!file.read((char*)&header, sizeof(header))) {
        throw tgt::CorruptedFileException("Failed to read header!", fileName);
    }

    file.close();

    //bool bigEndian = false;
    //check if swap is necessary:
    if((header.dim[0] < 0) || (header.dim[0] > 15)) {
        //bigEndian = true;
        header.swapEndianess();
    }

    ivec3 dimensions;
    dimensions.x = header.dim[1];
    dimensions.y = header.dim[2];
    dimensions.z = header.dim[3];

    int numVolumes = header.dim[4];
    int numChannels = header.dim[0] < 5 ? 1 : header.dim[5];

    vec3 spacing;
    spacing.x = header.pixdim[1];
    spacing.y = header.pixdim[2];
    spacing.z = header.pixdim[3];

    std::vector<VolumeURL> result;

    for(int i=0; i<numVolumes; i++) {
        for(int channel=0; channel<numChannels; channel++) {
            VolumeURL origin(/*"nii", */fileName);
            origin.addSearchParameter("volumeId", itos(i));
            origin.addSearchParameter("channel", itos(channel));
            origin.getMetaDataContainer().addMetaData("channel", new IntMetaData(channel));
            origin.getMetaDataContainer().addMetaData("volumeId", new IntMetaData(i));
            origin.getMetaDataContainer().addMetaData("Spacing", new Vec3MetaData(spacing));
            origin.getMetaDataContainer().addMetaData("Dimensions", new IVec3MetaData(dimensions));
            result.push_back(origin);
        }
    }

    return result;
}

} // namespace voreen
