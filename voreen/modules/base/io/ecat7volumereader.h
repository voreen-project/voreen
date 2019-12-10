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

#ifndef VRN_ECAT7VOLUMEREADER_H
#define VRN_ECAT7VOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

#ifdef WIN32
#pragma pack(1)
#endif
struct ECAT7VolumeHeader {
  char magic_number[14];
  char original_file_name[32];
  int16_t sw_version;
  int16_t system_type;
  int16_t file_type;
  char serial_number[10];
  uint32_t scan_start_time;
  char isotope_name[8];
  float isotope_halflife;
  char radiopharmaceutical[32];
  float gantry_tilt;
  float gantry_rotation;
  float bed_elevation;
  float intrinsic_tilt;
  int16_t wobble_speed;
  int16_t transm_source_type;
  float distance_scanned;
  float transaxial_fov;
  int16_t angular_compression;
  int16_t coin_samp_mode;
  int16_t axial_samp_mode;
  float ecat_calibration_factor;
  int16_t calibration_units;
  int16_t calibration_units_label;
  int16_t compression_code;
  char study_type[12];
  char patient_id[16];
  char patient_name[32];
  char patient_sex[1];
  char patient_dexterity[1];
  float patient_age;
  float patient_height;
  float patient_weight;
  uint32_t patient_birth_date;
  char physician_name[32];
  char operator_name[32];
  char study_description[32];
  int16_t acquisition_type;
  int16_t patient_orientation;
  char facility_name[20];
  uint16_t num_planes;
  uint16_t num_frames;
  uint16_t num_gates;
  uint16_t num_bed_pos;
  float init_bed_position;
  float bed_position[15];
  float plane_separation;
  int16_t lwr_sctr_thres;
  int16_t lwr_true_thres;
  int16_t upr_true_thres;
  char user_process_code[10];
  int16_t acquisition_mode;
  float bin_size;
  float branching_fraction;
  uint32_t dose_start_time;
  float dosage;
  float well_counter_corr_factor;
  char data_units[32];
  int16_t septa_state;
  int16_t fill_cti[6];

  void swapEndianness();
#ifdef __unix__
} __attribute__((packed));
#else
};
#endif

struct ECAT7DirectoryEntry {
    int16_t unkown_;
    int16_t id_;
    uint32_t startBlock_;
    uint32_t endBlock_;
    int32_t status_;

    void swapEndianness();
#ifdef __unix__
} __attribute__((packed));
#else
};
#endif

struct ECAT7ImageHeader {
  int16_t data_type;
  int16_t num_dimensions;
  uint16_t x_dimension;
  uint16_t y_dimension;
  uint16_t z_dimension;
  float x_offset;
  float y_offset;
  float z_offset;
  float recon_zoom;
  float scale_factor;
  int16_t image_min;
  int16_t image_max;
  float x_pixel_size;
  float y_pixel_size;
  float z_pixel_size;
  uint32_t frame_duration;
  uint32_t frame_start_time;
  int16_t filter_code;
  float x_resolution;
  float y_resolution;
  float z_resolution;
  float num_r_elements;
  float num_angles;
  float z_rotation_angle;
  float decay_corr_fctr;
  int32_t processing_code;
  uint32_t gate_duration;
  int32_t r_wave_offset;
  int32_t num_accepted_beats;
  float filter_cutoff_frequency;
  float filter_resolution;
  float filter_ramp_slope;
  int16_t filter_order;
  float filter_scatter_fraction;
  float filter_scatter_slope;
  char annotation[40];
  float mt_1_1;
  float mt_1_2;
  float mt_1_3;
  float mt_2_1;
  float mt_2_2;
  float mt_2_3;
  float mt_3_1;
  float mt_3_2;
  float mt_3_3;
  float rfilter_cutoff;
  float rfilter_resolution;
  int16_t rfilter_code;
  int16_t rfilter_order;
  float zfilter_cutoff;
  float zfilter_resolution;
  int16_t zfilter_code;
  int16_t zfilter_order;
  float mt_1_4;
  float mt_2_4;
  float mt_3_4;
  int16_t scatter_type;
  int16_t recon_type;
  int16_t recon_views;
  int16_t fill_cti[87];
  int16_t fill_user[49];


  void swapEndianness();
#ifdef __unix__
} __attribute__((packed));
#else
};
#endif
#ifdef WIN32
#pragma pack()
#endif

/**
 * Reader for ECAT 7 files from CTI / Siemens scanners
 */
class VRN_CORE_API ECAT7VolumeReader : public VolumeReader {
    struct SubVolume {
        ECAT7DirectoryEntry de_;
        ECAT7ImageHeader ih_;

        tgt::svec3 getDimensions();
        tgt::vec3 getSpacing();
        tgt::vec3 getOffset();
        tgt::mat4 getTransformation();
        int getId();

        void transformMetaData(MetaDataContainer& mdc);

        void printSubHeader() const;
    };

    struct ECAT7Structure {
        ECAT7VolumeHeader h_;
        std::vector<SubVolume> subVolumes_;
        bool swapEndianness_;

        void transformMetaData(MetaDataContainer& mdc, int volumeId);
    };

public:
    ECAT7VolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "ECAT7VolumeReader"; }
    virtual std::string getFormatDescription() const { return "CTI / Siemens ECAT 7 files"; }

    VolumeBase* read(const VolumeURL& origin);

    virtual VolumeList* read(const std::string& url);

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * \param   url         url to load volume from
     * \param   volumeId    id to select the volume, if -1 all volumes will be selected
     **/
    virtual VolumeList* read(const std::string& url, int volumeId);

    std::vector<VolumeURL> listVolumes(const std::string& url) const;

    ECAT7Structure readStructure(const std::string& url) const;

    void printMainHeader(const ECAT7VolumeHeader& h) const;


private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_ECAT7VOLUMEREADER_H
