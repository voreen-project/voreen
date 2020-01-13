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

#ifndef VRN_MULTISLICEVIEWER_H
#define VRN_MULTISLICEVIEWER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/volume/slice/slicehelper.h"

#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"
#include "voreen/core/interaction/slicecamerainteractionhandler.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"

namespace voreen {

class MultiSliceViewer : public VolumeRenderer {
    /// Determines the current axis-alignment if the displayed slices.
public:
    //TODO: use normal SliceCache
    class SliceCache {
        struct CacheEntry {
            tgt::plane plane_;
            std::string volumeId_; //both md5 hashes
            SliceTexture* slice_;

            CacheEntry(tgt::plane pl, std::string volumeId, SliceTexture* sl) : plane_(pl), volumeId_(volumeId), slice_(sl) {}
        };
        public:
        SliceCache(int cacheSize, float samplingRate);
        ~SliceCache();

        int getCacheSize() const { return cacheSize_; }
        void setCacheSize(int cacheSize);

        float getSamplingRate() const { return samplingRate_; }
        void setSamplingRate(float samplingRate);

        SliceTexture* getSlice(tgt::plane pl, const VolumeBase* vh) const;

        void clear();
    private:
        void setSlice(SliceTexture* m, tgt::plane pl, const VolumeBase* vh) const;

        mutable std::list<CacheEntry> slices_;
        int cacheSize_;
        float samplingRate_;
    };


    //TODO: also in slicerendererbase
    enum TextureMode {
        TEXTURE_2D,
        TEXTURE_3D
    };

    MultiSliceViewer();
    virtual ~MultiSliceViewer();

    virtual std::string getCategory() const { return "Volume"; }
    virtual std::string getClassName() const { return "MultiSliceViewer"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }
    virtual std::string getProcessorInfo() const;

    virtual Processor* create() const { return new MultiSliceViewer(); }

    virtual void initialize();
    virtual void deinitialize();

    virtual void process();
    virtual bool isReady() const;

    virtual void onEvent(tgt::Event* e);

protected:
    int getSliceIndex() const;
    const IntProperty* getSliceIndexProperty() const;
    void setIntProperty(IntProperty& p, int value);

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();
    virtual void setDescriptions() {
        setDescription("Displays 2D slices along one of the three main axis of one or more volumes. Rendering is performed in world space.");
    }

    const VolumePort* getMainInport() const;
    const std::vector<const VolumeBase*> getSecondaryVolumes() const;

    /// Adapts the min/max ranges of the respective properties to the dimensions of the currently connected volume.
    void updateSliceProperties();
    void alignCamera();
    void adjustPropertyVisibility();

    tgt::Shader* sh_;
    tgt::Shader* eepShader_;

    SliceCache sliceCache_;

    IntOptionProperty mainInport_;
    BoolProperty restrictToMainVolume_;
    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    IntProperty xSliceIndexProp_;
    IntProperty ySliceIndexProp_;
    IntProperty zSliceIndexProp_;
    CameraProperty camera_;
    ButtonProperty alignCameraButton_;

    BoolProperty renderCrosshair_;
    ColorProperty crosshairColor_;
    FloatProperty crosshairWidth_;
    FloatProperty crosshairRadius_;
    bool grabbedX_;
    bool grabbedY_;
    tgt::ivec2 lastMousePos_;

    FloatVec3Property plane_;
    FloatProperty planeDist_;

    TransFunc1DKeysProperty transferFunc1_;
    TransFunc1DKeysProperty transferFunc2_;
    TransFunc1DKeysProperty transferFunc3_;
    TransFunc1DKeysProperty transferFunc4_;

    StringOptionProperty blendingMode1_;
    StringOptionProperty blendingMode2_;
    StringOptionProperty blendingMode3_;
    StringOptionProperty blendingMode4_;

    FloatProperty blendingFactor1_;
    FloatProperty blendingFactor2_;
    FloatProperty blendingFactor3_;
    FloatProperty blendingFactor4_;

    IntOptionProperty texFilterMode1_;      ///< filter mode for volume 1
    IntOptionProperty texFilterMode2_;      ///< filter mode for volume 2
    IntOptionProperty texFilterMode3_;      ///< filter mode for volume 3
    IntOptionProperty texFilterMode4_;      ///< filter mode for volume 4

    GLEnumOptionProperty texClampMode1_;    ///< clamp mode for volume 1
    GLEnumOptionProperty texClampMode2_;    ///< clamp mode for volume 2
    GLEnumOptionProperty texClampMode3_;    ///< clamp mode for volume 3
    GLEnumOptionProperty texClampMode4_;    ///< clamp mode for volume 4
    FloatProperty texBorderIntensity_;      ///< border intensity for all volumes

    MWheelNumPropInteractionHandler<int> mwheelCycleHandler_;
    SliceCameraInteractionHandler interactionHandler_;

    OptionProperty<TextureMode> texMode1_;     ///< use 2D slice textures or 3D volume texture?
    OptionProperty<TextureMode> texMode2_;     ///< use 2D slice textures or 3D volume texture?
    OptionProperty<TextureMode> texMode3_;     ///< use 2D slice textures or 3D volume texture?
    OptionProperty<TextureMode> texMode4_;     ///< use 2D slice textures or 3D volume texture?
    FloatProperty samplingRate_;

    IntProperty cacheSize_;

    VolumePort inport1_;
    VolumePort inport2_;
    VolumePort inport3_;
    VolumePort inport4_;
    GeometryPort geomPort_;
    TextPort sliceIndexTextPort_;
    TextPort intensityTextPort_;
    RenderPort outport_;
    RenderPort entryPort_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
