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

#include "ecat7volumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorswapendianness.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormirror.h"

using tgt::ivec3;

namespace voreen {

const std::string ECAT7VolumeReader::loggerCat_ = "voreen.io.VolumeReader.ECAT7VolumeReader";

void ECAT7VolumeHeader::swapEndianness() {
    sw_version = swapEndian(sw_version);
    system_type = swapEndian(system_type);
    file_type = swapEndian(file_type);

    scan_start_time = swapEndian(scan_start_time);

    isotope_halflife = swapEndian(isotope_halflife);

    gantry_tilt = swapEndian(gantry_tilt);
    gantry_rotation = swapEndian(gantry_rotation);
    bed_elevation = swapEndian(bed_elevation);
    intrinsic_tilt = swapEndian(intrinsic_tilt);
    wobble_speed = swapEndian(wobble_speed);
    transm_source_type = swapEndian(transm_source_type);
    distance_scanned = swapEndian(distance_scanned);
    transaxial_fov = swapEndian(transaxial_fov);
    angular_compression = swapEndian(angular_compression);
    coin_samp_mode = swapEndian(coin_samp_mode);
    axial_samp_mode = swapEndian(axial_samp_mode);
    ecat_calibration_factor = swapEndian(ecat_calibration_factor);
    calibration_units = swapEndian(calibration_units);
    calibration_units_label = swapEndian(calibration_units_label);
    compression_code = swapEndian(compression_code);

    patient_age = swapEndian(patient_age);
    patient_height = swapEndian(patient_height);
    patient_weight = swapEndian(patient_weight);
    patient_birth_date = swapEndian(patient_birth_date);

    acquisition_type = swapEndian(acquisition_type);
    patient_orientation = swapEndian(patient_orientation);

    num_planes = swapEndian(num_planes);
    num_frames = swapEndian(num_frames);
    num_gates = swapEndian(num_gates);
    num_bed_pos = swapEndian(num_bed_pos);
    init_bed_position = swapEndian(init_bed_position);
    for (int i=0; i<15; i++) bed_position[i] = swapEndian(bed_position[i]);
    plane_separation = swapEndian(plane_separation);
    lwr_sctr_thres = swapEndian(lwr_sctr_thres);
    lwr_true_thres = swapEndian(lwr_true_thres);
    upr_true_thres = swapEndian(upr_true_thres);

    acquisition_mode = swapEndian(acquisition_mode);
    bin_size = swapEndian(bin_size);
    branching_fraction = swapEndian(branching_fraction);
    dose_start_time = swapEndian(dose_start_time);
    dosage = swapEndian(dosage);
    well_counter_corr_factor = swapEndian(well_counter_corr_factor);

    septa_state = swapEndian(septa_state);
}

void ECAT7ImageHeader::swapEndianness() {
    data_type = swapEndian(data_type);
    num_dimensions = swapEndian(num_dimensions);
    x_dimension = swapEndian(x_dimension);
    y_dimension = swapEndian(y_dimension);
    z_dimension = swapEndian(z_dimension);
    x_offset = swapEndian(x_offset);
    y_offset = swapEndian(y_offset);
    z_offset = swapEndian(z_offset);
    recon_zoom = swapEndian(recon_zoom);
    scale_factor = swapEndian(scale_factor);
    image_min = swapEndian(image_min);
    image_max = swapEndian(image_max);
    x_pixel_size = swapEndian(x_pixel_size);
    y_pixel_size = swapEndian(y_pixel_size);
    z_pixel_size = swapEndian(z_pixel_size);
    frame_duration = swapEndian(frame_duration);
    frame_start_time = swapEndian(frame_start_time);
    filter_code = swapEndian(filter_code);
    x_resolution = swapEndian(x_resolution);
    y_resolution = swapEndian(y_resolution);
    z_resolution = swapEndian(z_resolution);
    num_r_elements = swapEndian(num_r_elements);
    num_angles = swapEndian(num_angles);
    z_rotation_angle = swapEndian(z_rotation_angle);
    decay_corr_fctr = swapEndian(decay_corr_fctr);
    processing_code = swapEndian(processing_code);
    gate_duration = swapEndian(gate_duration);
    r_wave_offset = swapEndian(r_wave_offset);
    num_accepted_beats = swapEndian(num_accepted_beats);
    filter_cutoff_frequency = swapEndian(filter_cutoff_frequency);
    filter_resolution = swapEndian(filter_resolution);
    filter_ramp_slope = swapEndian(filter_ramp_slope);
    filter_order = swapEndian(filter_order);
    filter_scatter_fraction = swapEndian(filter_scatter_fraction);
    filter_scatter_slope = swapEndian(filter_scatter_slope);

    mt_1_1 = swapEndian(mt_1_1);
    mt_1_2 = swapEndian(mt_1_2);
    mt_1_3 = swapEndian(mt_1_3);
    mt_2_1 = swapEndian(mt_2_1);
    mt_2_2 = swapEndian(mt_2_2);
    mt_2_3 = swapEndian(mt_2_3);
    mt_3_1 = swapEndian(mt_3_1);
    mt_3_2 = swapEndian(mt_3_2);
    mt_3_3 = swapEndian(mt_3_3);
    rfilter_cutoff = swapEndian(rfilter_cutoff);
    rfilter_resolution = swapEndian(rfilter_resolution);
    rfilter_code = swapEndian(rfilter_code);
    rfilter_order = swapEndian(rfilter_order);
    zfilter_cutoff = swapEndian(zfilter_cutoff);
    zfilter_resolution = swapEndian(zfilter_resolution);
    zfilter_code = swapEndian(zfilter_code);
    zfilter_order = swapEndian(zfilter_order);
    mt_1_4 = swapEndian(mt_1_4);
    mt_2_4 = swapEndian(mt_2_4);
    mt_3_4 = swapEndian(mt_3_4);
    scatter_type = swapEndian(scatter_type);
    recon_type = swapEndian(recon_type);
    recon_views = swapEndian(recon_views);
    for (int i=0; i<87; i++) fill_cti[i] = swapEndian(fill_cti[i]);
    for (int i=0; i<49; i++) fill_user[i] = swapEndian(fill_user[i]);
}

//-----------------------------------------------------------------------------

void ECAT7VolumeReader::ECAT7Structure::transformMetaData(MetaDataContainer& mdc, int volumeId) {
    // Handled elsewhere or no metadata:
    //char magic_number[14];
    //char original_file_name[32];
    //int16_t sw_version;
    //int16_t system_type;
    //int16_t file_type;
    //int16_t fill_cti[6];
    //uint16_t num_frames;
    //int16_t compression_code;
    // Used for RealWorldMapping:
    //char data_units[32];
    //float ecat_calibration_factor;
    //int16_t calibration_units;
    //int16_t calibration_units_label;

    // Not sure about these:
    //char serial_number[10];
    //float gantry_tilt;
    //float gantry_rotation;
    //float bed_elevation;
    //float intrinsic_tilt;
    //int16_t wobble_speed;
    //int16_t transm_source_type;
    //float distance_scanned;
    //float transaxial_fov;
    //int16_t angular_compression;
    //int16_t coin_samp_mode;
    //int16_t axial_samp_mode;
    //float ecat_calibration_factor;
    //int16_t calibration_units;
    //int16_t calibration_units_label;
    //int16_t acquisition_type;
    //int16_t patient_orientation;
    //uint16_t num_planes;
    //uint16_t num_gates;
    //uint16_t num_bed_pos;
    //float init_bed_position;
    //float bed_position[15];
    //float plane_separation;
    //int16_t lwr_sctr_thres;
    //int16_t lwr_true_thres;
    //int16_t upr_true_thres;
    //char user_process_code[10];
    //int16_t acquisition_mode;
    //float bin_size;
    //float branching_fraction;
    //float well_counter_corr_factor;
    //int16_t septa_state;

    //TODO:
    //float patient_age;
    //char patient_sex[1];

    mdc.addMetaData("PatientDexterity", new StringMetaData(h_.patient_dexterity)); //char patient_dexterity[1];
    mdc.addMetaData("PatientHeight", new FloatMetaData(h_.patient_height)); //float patient_height;
    if(h_.patient_birth_date != 0)
        mdc.addMetaData("PatientBirthDate", new DateTimeMetaData(DateTime(h_.patient_birth_date))); //uint32_t patient_birth_date;

    if(h_.scan_start_time != 0)
        mdc.addMetaData("AcquisitionDateTime", new DateTimeMetaData(DateTime(h_.scan_start_time)));//uint32_t scan_start_time;
    if(h_.dose_start_time != 0)
        mdc.addMetaData("RadiopharmaceuticalStartTime", new DateTimeMetaData(DateTime(h_.dose_start_time)));//uint32_t dose_start_time;

    mdc.addMetaData("Isotope", new StringMetaData(h_.isotope_name));//char isotope_name[8];
    mdc.addMetaData("IsotopeHalflife", new FloatMetaData(h_.isotope_halflife));//float isotope_halflife;
    mdc.addMetaData("Radiopharmaceutical", new StringMetaData(h_.radiopharmaceutical));//char radiopharmaceutical[32];
    mdc.addMetaData("RadionuclideTotalDose", new FloatMetaData(h_.dosage));//float dosage;
    mdc.addMetaData("OperatorsName", new StringMetaData(h_.operator_name));//char operator_name[32];
    mdc.addMetaData("PerformingPhysiciansName", new StringMetaData(h_.physician_name));//char physician_name[32];
    mdc.addMetaData("PatientName", new StringMetaData(h_.patient_name));//char patient_name[32];
    mdc.addMetaData("PatientID", new StringMetaData(h_.patient_id));//char patient_id[16];
    mdc.addMetaData("PatientWeight", new FloatMetaData(h_.patient_weight));//float patient_weight;

    mdc.addMetaData("InstitutionName", new StringMetaData(h_.facility_name));//char facility_name[20];
    mdc.addMetaData("StudyDescription", new StringMetaData(h_.study_description));//char study_description[32];
    mdc.addMetaData("StudyType", new StringMetaData(h_.study_type));//char study_type[12];

    subVolumes_[volumeId].transformMetaData(mdc);
}

//-----------------------------------------------------------------------------

tgt::svec3 ECAT7VolumeReader::SubVolume::getDimensions() {
    tgt::svec3 dimensions;
    dimensions.x = ih_.x_dimension;
    dimensions.y = ih_.y_dimension;
    dimensions.z = ih_.z_dimension;
    return dimensions;
}

tgt::vec3 ECAT7VolumeReader::SubVolume::getSpacing() {
    tgt::vec3 spacing;
    spacing.x = ih_.x_pixel_size;
    spacing.y = ih_.y_pixel_size;
    spacing.z = ih_.z_pixel_size;
    return spacing;
}

tgt::vec3 ECAT7VolumeReader::SubVolume::getOffset() {
    tgt::vec3 offset;
    offset.x = ih_.x_offset;
    offset.y = ih_.y_offset;
    offset.z = ih_.z_offset;
    return offset;
}

tgt::mat4 ECAT7VolumeReader::SubVolume::getTransformation() {
    tgt::mat4 m = tgt::mat4::identity;

    m[0] = tgt::vec4(ih_.mt_1_1, ih_.mt_1_2, ih_.mt_1_3, 0.f);
    m[1] = tgt::vec4(ih_.mt_2_1, ih_.mt_2_2, ih_.mt_2_3, 0.f);
    m[2] = tgt::vec4(ih_.mt_3_1, ih_.mt_3_2, ih_.mt_3_3, 0.f);

    if(length(m[0]) == 0.f && length(m[1]) == 0.f && length(m[2]) == 0.f)
        m = tgt::mat4::identity;

    return m;
}

int ECAT7VolumeReader::SubVolume::getId() {
    return de_.id_;
}

void ECAT7VolumeReader::SubVolume::transformMetaData(MetaDataContainer& mdc) {
  // Handled elsewhere or no metadata:
  //int16_t data_type;
  //int16_t num_dimensions;
  //uint16_t x_dimension;
  //uint16_t y_dimension;
  //uint16_t z_dimension;
  //float x_offset;
  //float y_offset;
  //float z_offset;
  //int16_t image_min;
  //int16_t image_max;
  //float x_pixel_size;
  //float y_pixel_size;
  //float z_pixel_size;
  //int16_t fill_cti[87];
  //int16_t fill_user[49];
  //float mt_1_1;
  //float mt_1_2;
  //float mt_1_3;
  //float mt_2_1;
  //float mt_2_2;
  //float mt_2_3;
  //float mt_3_1;
  //float mt_3_2;
  //float mt_3_3;
  //float mt_1_4;
  //float mt_2_4;
  //float mt_3_4;
  //float scale_factor;

  // Not sure about these:
  //float x_resolution;
  //float y_resolution;
  //float z_resolution;
  //float recon_zoom;
  //float filter_cutoff_frequency;
  //float filter_resolution;
  //float filter_ramp_slope;
  //int16_t filter_order;
  //float filter_scatter_fraction;
  //float filter_scatter_slope;
  //float rfilter_cutoff;
  //float rfilter_resolution;
  //int16_t rfilter_code;
  //int16_t rfilter_order;
  //float zfilter_cutoff;
  //float zfilter_resolution;
  //int16_t zfilter_code;
  //int16_t zfilter_order;
  //int16_t scatter_type;
  //int16_t recon_type;
  //int16_t recon_views;

  mdc.addMetaData("ActualFrameDuration", new IntMetaData(ih_.frame_duration));//uint32_t frame_duration;
  mdc.addMetaData("FrameTime", new IntMetaData(ih_.frame_start_time));//uint32_t frame_start_time;
  mdc.addMetaData("SeriesDescription", new StringMetaData(ih_.annotation));//char annotation[40];

  //TODO:
  //int16_t filter_code;
  //float num_r_elements;
  //float num_angles;
  //float z_rotation_angle;
  //float decay_corr_fctr;
  //int32_t processing_code;
  //uint32_t gate_duration;
  //int32_t r_wave_offset;
  //int32_t num_accepted_beats;
}

void ECAT7VolumeReader::SubVolume::printSubHeader() const {
    LINFO("SUBHEADER INFO:");
    LINFO(" data_type:               " << ih_.data_type);
    LINFO(" num_dimensions:          " << ih_.num_dimensions);
    LINFO(" x_dimension:             " << ih_.x_dimension);
    LINFO(" y_dimension:             " << ih_.y_dimension);
    LINFO(" z_dimension:             " << ih_.z_dimension);
    LINFO(" x_offset:                " << ih_.x_offset);
    LINFO(" y_offset:                " << ih_.y_offset);
    LINFO(" z_offset:                " << ih_.z_offset);
    LINFO(" recon_zoom:              " << ih_.recon_zoom);
    LINFO(" scale_factor:            " << ih_.scale_factor);
    LINFO(" image_min:               " << ih_.image_min);
    LINFO(" image_max:               " << ih_.image_max);
    LINFO(" x_pixel_size:            " << ih_.x_pixel_size);
    LINFO(" y_pixel_size:            " << ih_.y_pixel_size);
    LINFO(" z_pixel_size:            " << ih_.z_pixel_size);
    LINFO(" frame_duration:          " << ih_.frame_duration);
    LINFO(" frame_start_time:        " << ih_.frame_start_time);
    LINFO(" filter_code:             " << ih_.filter_code);
    LINFO(" x_resolution:            " << ih_.x_resolution);
    LINFO(" y_resolution:            " << ih_.y_resolution);
    LINFO(" z_resolution:            " << ih_.z_resolution);
    LINFO(" num_r_elements:          " << ih_.num_r_elements);
    LINFO(" num_angles:              " << ih_.num_angles);
    LINFO(" z_rotation_angle:        " << ih_.z_rotation_angle);
    LINFO(" decay_corr_fctr:         " << ih_.decay_corr_fctr);
    LINFO(" processing_code:         " << ih_.processing_code);
    LINFO(" gate_duration:           " << ih_.gate_duration);
    LINFO(" r_wave_offset:           " << ih_.r_wave_offset);
    LINFO(" num_accepted_beats:      " << ih_.num_accepted_beats);
    LINFO(" filter_cutoff_frequency: " << ih_.filter_cutoff_frequency);
    LINFO(" filter_resolution:       " << ih_.filter_resolution);
    LINFO(" filter_ramp_slope:       " << ih_.filter_ramp_slope);
    LINFO(" filter_order:            " << ih_.filter_order);
    LINFO(" filter_scatter_fraction: " << ih_.filter_scatter_fraction);
    LINFO(" filter_scatter_slope:    " << ih_.filter_scatter_slope);
    LINFO(" annotation:              " << std::string(ih_.annotation));
    LINFO(" mt_1_1:                  " << ih_.mt_1_1);
    LINFO(" mt_1_2:                  " << ih_.mt_1_2);
    LINFO(" mt_1_3:                  " << ih_.mt_1_3);
    LINFO(" mt_2_1:                  " << ih_.mt_2_1);
    LINFO(" mt_2_2:                  " << ih_.mt_2_2);
    LINFO(" mt_2_3:                  " << ih_.mt_2_3);
    LINFO(" mt_3_1:                  " << ih_.mt_3_1);
    LINFO(" mt_3_2:                  " << ih_.mt_3_2);
    LINFO(" mt_3_3:                  " << ih_.mt_3_3);
    LINFO(" rfilter_cutoff:          " << ih_.rfilter_cutoff);
    LINFO(" rfilter_resolution:      " << ih_.rfilter_resolution);
    LINFO(" rfilter_code:            " << ih_.rfilter_code);
    LINFO(" rfilter_order:           " << ih_.rfilter_order);
    LINFO(" zfilter_cutoff:          " << ih_.zfilter_cutoff);
    LINFO(" zfilter_resolution:      " << ih_.zfilter_resolution);
    LINFO(" zfilter_code:            " << ih_.zfilter_code);
    LINFO(" zfilter_order:           " << ih_.zfilter_order);
    LINFO(" mt_1_4:                  " << ih_.mt_1_4);
    LINFO(" mt_2_4:                  " << ih_.mt_2_4);
    LINFO(" mt_3_4:                  " << ih_.mt_3_4);
    LINFO(" scatter_type:            " << ih_.scatter_type);
    LINFO(" recon_type:              " << ih_.recon_type);
    LINFO(" recon_views:             " << ih_.recon_views);
    LINFO(" fill_cti:                " << ih_.fill_cti);
    LINFO(" int16_fill_user:         " << ih_.fill_user);
}

//-----------------------------------------------------------------------------

void ECAT7DirectoryEntry::swapEndianness() {
    id_ = swapEndian(id_);
    startBlock_ = swapEndian(startBlock_);
    endBlock_ = swapEndian(endBlock_);
    status_ = swapEndian(status_);
}

//-----------------------------------------------------------------------------

ECAT7VolumeReader::ECAT7VolumeReader(ProgressBar* progress)
    : VolumeReader(progress) {
    extensions_.push_back("v");
    protocols_.push_back("ecat7");
}

VolumeBase* ECAT7VolumeReader::read(const VolumeURL& origin) {
    VolumeBase* result = 0;

    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (! tmp.empty())
        volumeId = stoi(tmp);

    VolumeList* collection = read(origin.getPath(), volumeId);

    if (collection && collection->size() == 1) {
        result = collection->first();
    }
    else if (collection && collection->size() > 1) {
        while(!collection->empty()) {
           VolumeBase* vh = collection->first();
           collection->remove(vh);
           delete vh;
        }
        delete collection;
        throw tgt::FileException("Only one volume expected", origin.getPath());
    }

    delete collection;

    return result;
}

VolumeList* ECAT7VolumeReader::read(const std::string &url) {
    VolumeURL origin(url);

    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (! tmp.empty())
        volumeId = stoi(tmp);

    return read(url, volumeId);
}

VolumeList* ECAT7VolumeReader::read(const std::string &url, int volumeId) {
    LINFO("Loading dataset " << url << " vid: " << volumeId);

    VolumeURL origin(url);
    std::string fileName = origin.getPath();
    FILE* fin = fopen(fileName.c_str(), "rb");

    if(!fin) {
        throw tgt::FileNotFoundException("ECAT7: File not found", fileName);
    }

    ECAT7VolumeReader::ECAT7Structure s = readStructure(fileName);

    VolumeList* vc = new VolumeList();

    for(size_t i=0; i<s.subVolumes_.size(); i++) {
        fseek(fin, s.subVolumes_[i].de_.startBlock_ * 512, SEEK_SET);
        fseek(fin, 512, SEEK_CUR); //skip past header (already read)

        tgt::svec3 dimensions = s.subVolumes_[i].getDimensions();
        //tgt::mat4 m = s.subVolumes_[i].getTransformation();
        tgt::vec3 spacing = s.subVolumes_[i].getSpacing();
        tgt::vec3 offset  = s.subVolumes_[i].getOffset();

        if(volumeId != -1) {
            if(volumeId != s.subVolumes_[i].getId())
                continue;
        }

        if(s.subVolumes_[i].ih_.num_dimensions != 3)
            continue;

        if (getProgressBar()) {
            getProgressBar()->setTitle("Loading Volume");
            getProgressBar()->setProgressMessage("Loading volume: " + fileName);
        }

        VolumeRAM* vol;
        RealWorldMapping denormalize;
        if(s.h_.file_type == 6) {
            vol = new VolumeRAM_UInt8(dimensions);
            denormalize = RealWorldMapping::createDenormalizingMapping<uint8_t>();
        }
        else if(s.h_.file_type == 7) {
            vol = new VolumeRAM_UInt16(dimensions);
            denormalize = RealWorldMapping::createDenormalizingMapping<uint16_t>();
        }
        else {
            LERROR("Unknown file format detected.");
            return 0;
        }

        float scale = s.subVolumes_[i].ih_.scale_factor * s.h_.ecat_calibration_factor;
        RealWorldMapping rwm(scale, 0.0f, s.h_.data_units);

        VolumeReader::read(vol, fin);

        // Assume that the pixel size values given in the ecat header are in cm.  Multiply spacing and offset
        // with 0.1 to convert to mm.
        Volume* vh = new Volume(vol, spacing * 0.1f, offset * 0.1f);
        vh->setRealWorldMapping(RealWorldMapping::combine(denormalize, rwm));

        if(s.swapEndianness_)
            VolumeOperatorSwapEndianness::APPLY_OP(vh);

        // TODO: This must depend on some parameter in the headers, figure out which and how
        bool mirrorZ = true;
        if(mirrorZ)
        {
            Volume* mirrored = VolumeOperatorMirrorZ::APPLY_OP(vh);
            delete vh;
            vh = mirrored;
        }

        VolumeURL o("ecat7", fileName);
        o.addSearchParameter("volumeId", itos(s.subVolumes_[i].de_.id_));
        vh->setOrigin(o);

        s.transformMetaData(vh->getMetaDataContainer(), static_cast<int>(i));

        if(length(offset) == 0.f)
            centerVolume(vh);

        vc->add(vh);

        if (getProgressBar())
            getProgressBar()->hide();
    }

    fclose(fin);
    return vc;
}

std::vector<VolumeURL> ECAT7VolumeReader::listVolumes(const std::string& url) const {
    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();
    ECAT7VolumeReader::ECAT7Structure s = readStructure(fileName);

    std::vector<VolumeURL> result;

    for(size_t i=0; i<s.subVolumes_.size(); i++) {
        VolumeURL origin("ecat7", fileName);
        origin.addSearchParameter("volumeId", itos(s.subVolumes_[i].de_.id_));
        origin.getMetaDataContainer().addMetaData("volumeId", new IntMetaData(s.subVolumes_[i].de_.id_));
        origin.getMetaDataContainer().addMetaData("Spacing", new Vec3MetaData(s.subVolumes_[i].getSpacing()));
        origin.getMetaDataContainer().addMetaData("Dimensions", new IVec3MetaData(s.subVolumes_[i].getDimensions()));
        s.transformMetaData(origin.getMetaDataContainer(), static_cast<int>(i));
        result.push_back(origin);
    }

    return result;
}

ECAT7VolumeReader::ECAT7Structure ECAT7VolumeReader::readStructure(const std::string& url) const {

    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    ECAT7Structure s;
    FILE* fin = fopen(fileName.c_str(), "rb");

    if(!fin) {
        throw tgt::FileNotFoundException("ECAT7: File not found", fileName);
    }

    tgtAssert(sizeof(ECAT7VolumeHeader) == 512, "Invalid volume header sizeof!");
    tgtAssert(sizeof(ECAT7ImageHeader) == 512, "Invalid image header sizeof!");

    fread(&s.h_, 512, 1, fin);

    s.swapEndianness_ = true; //TODO: check if necessary

    if(s.swapEndianness_)
        s.h_.swapEndianness();

    std::string magic(s.h_.magic_number, 14);

    if(magic.substr(0, 7) != "MATRIX7")
        throw tgt::CorruptedFileException("error while reading data", fileName);

    if((s.h_.file_type != 6) && (s.h_.file_type != 7)) {
        throw tgt::UnsupportedFormatException("Unknown file_type. Only Volume8 and Volume16 supported.");
    }

    //printMainHeader(s.h_);

    int32_t numMatrices;
    fread(&numMatrices, 4, 1, fin); //TODO?
    fread(&numMatrices, 4, 1, fin); //TODO?
    fread(&numMatrices, 4, 1, fin); //TODO?

    fread(&numMatrices, 4, 1, fin);
    if(s.swapEndianness_)
        numMatrices = swapEndian(numMatrices);

    ECAT7DirectoryEntry dir;
    for(int i=0; i<numMatrices; i++) {
        fread(&dir, sizeof(ECAT7DirectoryEntry), 1, fin);

        if(s.swapEndianness_)
            dir.swapEndianness();

        s.subVolumes_.push_back(SubVolume());
        s.subVolumes_.back().de_ = dir;
        s.subVolumes_.back().de_.startBlock_--;
        s.subVolumes_.back().de_.endBlock_--;
    }

    ECAT7ImageHeader ih;
    for(int i=0; i<numMatrices; i++) {
        fseek(fin, s.subVolumes_[i].de_.startBlock_ * 512, SEEK_SET);

        fread(&ih, 512, 1, fin);

        if(s.swapEndianness_)
            ih.swapEndianness();

        //if(ih.num_dimensions != 3) //TODO
            //continue;

        s.subVolumes_[i].ih_ = ih;
        //s.subVolumes_[i].printSubHeader();
    }

    fclose(fin);
    return s;
}

VolumeReader* ECAT7VolumeReader::create(ProgressBar* progress) const {
    return new ECAT7VolumeReader(progress);
}

void ECAT7VolumeReader::printMainHeader(const ECAT7VolumeHeader& h) const {
    LINFO("HEADER INFO:");
    LINFO("    magic_number              " << h.magic_number);
    LINFO("    original_file_name        " << h.original_file_name);
    LINFO("    sw_version                " << h.sw_version);
    LINFO("    system_type               " << h.system_type);
    LINFO("    file_type                 " << h.file_type);
    LINFO("    serial_number             " << h.serial_number);
    LINFO("    scan_start_time           " << h.scan_start_time);
    LINFO("    isotope_name              " << h.isotope_name);
    LINFO("    isotope_halflife          " << h.isotope_halflife);
    LINFO("    radiopharmaceutical       " << h.radiopharmaceutical);
    LINFO("    gantry_tilt               " << h.gantry_tilt);
    LINFO("    gantry_rotation           " << h.gantry_rotation);
    LINFO("    bed_elevation             " << h.bed_elevation);
    LINFO("    intrinsic_tilt            " << h.intrinsic_tilt);
    LINFO("    wobble_speed              " << h.wobble_speed);
    LINFO("    transm_source_type        " << h.transm_source_type);
    LINFO("    distance_scanned          " << h.distance_scanned);
    LINFO("    transaxial_fov            " << h.transaxial_fov);
    LINFO("    angular_compression       " << h.angular_compression);
    LINFO("    coin_samp_mode            " << h.coin_samp_mode);
    LINFO("    axial_samp_mode           " << h.axial_samp_mode);
    LINFO("    ecat_calibration_factor   " << h.ecat_calibration_factor);
    LINFO("    calibration_units         " << h.calibration_units);
    LINFO("    calibration_units_label   " << h.calibration_units_label);
    LINFO("    compression_code          " << h.compression_code);
    LINFO("    study_type                " << h.study_type);
    LINFO("    patient_id                " << h.patient_id);
    LINFO("    patient_name              " << h.patient_name);
    LINFO("    patient_sex               " << h.patient_sex);
    LINFO("    patient_dexterity         " << h.patient_dexterity);
    LINFO("    patient_age               " << h.patient_age);
    LINFO("    patient_height            " << h.patient_height);
    LINFO("    patient_weight            " << h.patient_weight);
    LINFO("    patient_birth_date        " << h.patient_birth_date);
    LINFO("    physician_name            " << h.physician_name);
    LINFO("    operator_name             " << h.operator_name);
    LINFO("    study_description         " << h.study_description);
    LINFO("    acquisition_type          " << h.acquisition_type);
    LINFO("    patient_orientation       " << h.patient_orientation);
    LINFO("    facility_name             " << h.facility_name);
    LINFO("    num_planes                " << h.num_planes);
    LINFO("    num_frames                " << h.num_frames);
    LINFO("    num_gates                 " << h.num_gates);
    LINFO("    num_bed_pos               " << h.num_bed_pos);
    LINFO("    init_bed_position         " << h.init_bed_position);
    LINFO("    bed_position              " << h.bed_position);
    LINFO("    plane_separation          " << h.plane_separation);
    LINFO("    lwr_sctr_thres            " << h.lwr_sctr_thres);
    LINFO("    lwr_true_thres            " << h.lwr_true_thres);
    LINFO("    upr_true_thres            " << h.upr_true_thres);
    LINFO("    user_process_code         " << h.user_process_code);
    LINFO("    acquisition_mode          " << h.acquisition_mode);
    LINFO("    bin_size                  " << h.bin_size);
    LINFO("    branching_fraction        " << h.branching_fraction);
    LINFO("    dose_start_time           " << h.dose_start_time);
    LINFO("    dosage                    " << h.dosage);
    LINFO("    well_counter_corr_factor  " << h.well_counter_corr_factor);
    LINFO("    data_units                " << h.data_units);
    LINFO("    septa_state               " << h.septa_state);
    LINFO("    fill_cti                  " << h.fill_cti);
}

} // namespace voreen
