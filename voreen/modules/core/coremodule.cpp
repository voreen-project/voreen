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

#include "coremodule.h"

// input processors
#include "processors/input/geometrysource.h"
#include "processors/input/imagesequencesource.h"
#include "processors/input/imagesource.h"
#include "processors/input/imageselector.h"
#include "processors/input/octreecreator.h"
#include "processors/input/textsource.h"
#include "processors/input/volumelistsource.h"
#include "processors/input/volumesource.h"
#include "processors/input/volumeselector.h"

// output processors
#include "processors/output/canvasrenderer.h"
#include "processors/output/geometrysave.h"
#include "processors/output/imagesequencesave.h"
#include "processors/output/textsave.h"
#include "processors/output/volumesave.h"
#include "processors/output/octreesave.h"
#include "processors/output/volumelistsave.h"

// ports
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"

// properties
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/collectivesettingsproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/planeproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/propertyvector.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/properties/volumeurlproperty.h"
#include "voreen/core/properties/voxeltypeproperty.h"
#include "voreen/core/properties/color/colorswitchproperty.h"

// property link evaluators
#include "voreen/core/properties/link/linkevaluatorboolinvert.h"
#include "voreen/core/properties/link/linkevaluatorboundingbox.h"
#include "voreen/core/properties/link/linkevaluatorcolorswitch.h"
#include "voreen/core/properties/link/linkevaluatorid.h"
#include "voreen/core/properties/link/linkevaluatorinterval.h"
#include "voreen/core/properties/link/linkevaluatormatrixinvert.h"
#include "voreen/core/properties/link/linkevaluatortransferfunctionproperties.h"

// geometry
#include "voreen/core/datastructures/geometry/vertexgeometry.h"
#include "voreen/core/datastructures/geometry/facegeometry.h"
#include "voreen/core/datastructures/geometry/meshgeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"

// meta data
#include "voreen/core/datastructures/meta/metadatabase.h"
#include "voreen/core/datastructures/meta/positionmetadata.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/meta/serializablevectormetadata.h"
#include "voreen/core/datastructures/meta/windowstatemetadata.h"

// octree
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagermmap.h"

// transfer functions
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/utils/transfuncmappingkey.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/transfunc1dgaussian.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/utils/transfuncmappingcurve.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/utils/transfuncprimitive.h"

// volume derived data
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/datastructures/volume/histogram.h"

// ROI
#include "voreen/core/datastructures/roi/roicube.h"
#include "voreen/core/datastructures/roi/roicylinder.h"
#include "voreen/core/datastructures/roi/roisphere.h"
#include "voreen/core/datastructures/roi/roiraster.h"
#include "voreen/core/datastructures/roi/roiunion.h"
#include "voreen/core/datastructures/roi/roisubtract.h"
#include "voreen/core/datastructures/roi/roigraph.h"

// volume i/o
#include "io/datvolumereader.h"
#include "io/datvolumewriter.h"
#include "io/rawvolumereader.h"
#include "io/vvdvolumereader.h"
#include "io/vvdvolumewriter.h"
#include "io/vvodvolumereader.h"
#include "io/vvodvolumewriter.h"

// volume operator
#include "voreen/core/datastructures/volume/operators/volumeoperatorcalcerror.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconnectedcomponents.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorcurvature.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorequalize.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgaussian.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorhalfsample.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorinvert.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorisuniform.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormagnitude.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormedian.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormirror.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormorphology.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatornumsignificant.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorregiongrow.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresizepoweroftwo.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorsecondderivatives.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorsubset.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorswapendianness.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatortranspose.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatoruncertaintymeasure.h"

#include "voreen/core/animation/animation.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string CoreModule::loggerCat_("voreen.CoreModule");

CoreModule::CoreModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Core");
    setGuiName("Core");

    tgtAssert(VoreenApplication::app(), "VoreenApplicaton not instantiated");

    addShaderPath(getModulePath("glsl"));
    addShaderPath(getModulePath("glsl/utils"));
    addShaderPath(getModulePath("glsl/geometry"));

    // processors
    registerSerializableType(new VolumeSource());
    registerSerializableType(new VolumeListSource());
    registerSerializableType(new VolumeSelector());
    registerSerializableType(new ImageSource());
    registerSerializableType(new ImageSequenceSource());
    registerSerializableType(new ImageSelector());
    registerSerializableType(new GeometrySource());
    registerSerializableType(new TextSource());
    registerSerializableType(new OctreeCreator());

    registerSerializableType(new CanvasRenderer());
    registerSerializableType(new VolumeSave());
    registerSerializableType(new OctreeSave());
    registerSerializableType(new VolumeListSave());
    registerSerializableType(new ImageSequenceSave());
    registerSerializableType(new GeometrySave());
    registerSerializableType(new TextSave());

    // ports
    registerSerializableType(new VolumePort(Port::OUTPORT, "dummy"));
    registerSerializableType(new RenderPort(Port::OUTPORT, "dummy"));
    registerSerializableType(new GeometryPort(Port::OUTPORT, "dummy"));
    registerSerializableType(new TextPort(Port::OUTPORT, "dummy"));
    registerSerializableType(new VolumeListPort(Port::OUTPORT, "dummy"));
    registerSerializableType(new ImageSequencePort(Port::OUTPORT, "dummy"));
    registerSerializableType(new LoopPort(Port::OUTPORT, "dummy"));

    // properties
    registerSerializableType(new FloatProperty());
    registerSerializableType(new IntProperty());
    registerSerializableType(new IntVec2Property());
    registerSerializableType(new IntVec3Property());
    registerSerializableType(new IntVec4Property());
    registerSerializableType(new FloatVec2Property());
    registerSerializableType(new FloatVec3Property());
    registerSerializableType(new FloatVec4Property());
    registerSerializableType(new FloatMat2Property());
    registerSerializableType(new FloatMat3Property());
    registerSerializableType(new FloatMat4Property());
    registerSerializableType(new ColorProperty());

    registerSerializableType(new BoolProperty());
    registerSerializableType(new FloatBoundingBoxProperty());
    registerSerializableType(new IntBoundingBoxProperty());
    registerSerializableType(new ButtonProperty());
    registerSerializableType(new CameraProperty());
    registerSerializableType(new FileDialogProperty());
    registerSerializableType(new FontProperty());
    registerSerializableType(dynamic_cast<Property*>(new PropertyVector()));
    registerSerializableType(new ProgressProperty());
    registerSerializableType(new ShaderProperty());
    registerSerializableType(new StringExpressionProperty());
    registerSerializableType(new StringProperty());
    registerSerializableType(new TransFunc1DKeysProperty());
    registerSerializableType(new TransFunc1DGaussianProperty());
    registerSerializableType(new VolumeURLListProperty());
    registerSerializableType(new VolumeURLProperty());
    registerSerializableType(dynamic_cast<Property*>(new VoxelTypeProperty()));

    registerSerializableType(new IntOptionProperty());
    registerSerializableType(new FloatOptionProperty());
    registerSerializableType(new GLEnumOptionProperty());
    registerSerializableType(new StringOptionProperty());
    registerSerializableType(new CollectiveSettingsProperty());
    registerSerializableType(new PlaneProperty());

    // -----------------
    //  Link Evaluators
    // -----------------

    // id
    registerSerializableType(new LinkEvaluatorBoolId());
    registerSerializableType(new LinkEvaluatorColorSwitchId());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxId());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxId());
    registerSerializableType(new LinkEvaluatorIntIntervalId());
    registerSerializableType(new LinkEvaluatorFloatIntervalId());

    registerSerializableType(new LinkEvaluatorIntId());
    registerSerializableType(new LinkEvaluatorFloatId());
    registerSerializableType(new LinkEvaluatorDoubleId());

    registerSerializableType(new LinkEvaluatorIntIdBounds());
    registerSerializableType(new LinkEvaluatorFloatIdBounds());
    registerSerializableType(new LinkEvaluatorDoubleIdBounds());

    registerSerializableType(new LinkEvaluatorRenderSize());      //< specialized ivec2 id link for rendering sizes
    registerSerializableType(new LinkEvaluatorIVec2Id());
    registerSerializableType(new LinkEvaluatorIVec3Id());
    registerSerializableType(new LinkEvaluatorIVec4Id());

    registerSerializableType(new LinkEvaluatorLightSourceId());
    registerSerializableType(new LinkEvaluatorVec2Id());
    registerSerializableType(new LinkEvaluatorVec3Id());
    registerSerializableType(new LinkEvaluatorVec4Id());

    registerSerializableType(new LinkEvaluatorDVec2Id());
    registerSerializableType(new LinkEvaluatorDVec3Id());
    registerSerializableType(new LinkEvaluatorDVec4Id());

    registerSerializableType(new LinkEvaluatorMat2Id());
    registerSerializableType(new LinkEvaluatorMat3Id());
    registerSerializableType(new LinkEvaluatorMat4Id());

    registerSerializableType(new LinkEvaluatorFontId());
    registerSerializableType(new LinkEvaluatorIntListId());

    // id conversion
    registerSerializableType(new LinkEvaluatorDoubleFloatId());
    registerSerializableType(new LinkEvaluatorDoubleIntId());
    registerSerializableType(new LinkEvaluatorDoubleBoolId());
    registerSerializableType(new LinkEvaluatorFloatIntId());
    registerSerializableType(new LinkEvaluatorFloatBoolId());
    registerSerializableType(new LinkEvaluatorIntBoolId());

    registerSerializableType(new LinkEvaluatorDVec2IVec2Id());
    registerSerializableType(new LinkEvaluatorDVec3IVec3Id());
    registerSerializableType(new LinkEvaluatorDVec4IVec4Id());

    registerSerializableType(new LinkEvaluatorDVec2Vec2Id());
    registerSerializableType(new LinkEvaluatorDVec3Vec3Id());
    registerSerializableType(new LinkEvaluatorDVec4Vec4Id());

    registerSerializableType(new LinkEvaluatorVec2IVec2Id());
    registerSerializableType(new LinkEvaluatorVec3IVec3Id());
    registerSerializableType(new LinkEvaluatorVec4IVec4Id());

    // id normalized
    registerSerializableType(new LinkEvaluatorIntIdNormalized());
    registerSerializableType(new LinkEvaluatorFloatIdNormalized());
    registerSerializableType(new LinkEvaluatorDoubleIdNormalized());

    registerSerializableType(new LinkEvaluatorIVec2IdNormalized());
    registerSerializableType(new LinkEvaluatorIVec3IdNormalized());
    registerSerializableType(new LinkEvaluatorIVec4IdNormalized());

    registerSerializableType(new LinkEvaluatorVec2IdNormalized());
    registerSerializableType(new LinkEvaluatorVec3IdNormalized());
    registerSerializableType(new LinkEvaluatorVec4IdNormalized());

    registerSerializableType(new LinkEvaluatorDVec2IdNormalized());
    registerSerializableType(new LinkEvaluatorDVec3IdNormalized());
    registerSerializableType(new LinkEvaluatorDVec4IdNormalized());

    registerSerializableType(new LinkEvaluatorMat2IdNormalized());
    registerSerializableType(new LinkEvaluatorMat3IdNormalized());
    registerSerializableType(new LinkEvaluatorMat4IdNormalized());

    // id normalized conversion
    registerSerializableType(new LinkEvaluatorDoubleFloatIdNormalized());
    registerSerializableType(new LinkEvaluatorDoubleIntIdNormalized());
    registerSerializableType(new LinkEvaluatorFloatIntIdNormalized());
    registerSerializableType(new LinkEvaluatorBoolInvert());

    registerSerializableType(new LinkEvaluatorDVec2IVec2IdNormalized());
    registerSerializableType(new LinkEvaluatorDVec3IVec3IdNormalized());
    registerSerializableType(new LinkEvaluatorDVec4IVec4IdNormalized());

    registerSerializableType(new LinkEvaluatorDVec2Vec2IdNormalized());
    registerSerializableType(new LinkEvaluatorDVec3Vec3IdNormalized());
    registerSerializableType(new LinkEvaluatorDVec4Vec4IdNormalized());

    registerSerializableType(new LinkEvaluatorVec2IVec2IdNormalized());
    registerSerializableType(new LinkEvaluatorVec3IVec3IdNormalized());
    registerSerializableType(new LinkEvaluatorVec4IVec4IdNormalized());

    //registerSerializableType(new LinkEvaluatorMat2Invert());
    registerSerializableType(new LinkEvaluatorMat3Invert());
    registerSerializableType(new LinkEvaluatorMat4Invert());

    registerSerializableType(new LinkEvaluatorTransFunc1DThresholdMin());
    registerSerializableType(new LinkEvaluatorTransFunc1DThresholdMax());
    registerSerializableType(new LinkEvaluatorTransFunc1DDomainMin());
    registerSerializableType(new LinkEvaluatorTransFunc1DDomainMax());
    registerSerializableType(new LinkEvaluatorTransFunc1DGamma());

    //Links between scalars and vectors

// Declare a macro to register all the vec to scalar link classes (value) using OP_VEC_TO_SCALAR_ALL declared in linkevaluatorid.h
#define REGISTER_VEC_TO_SCALAR_LINK(scalarType,vecType,vecDim,vecComponent) \
    registerSerializableType(new VEC_TO_SCALAR_LINK_CLASSNAME(scalarType, vecType, vecDim, vecComponent) ());

    // Use the macros to register 4x4x(2+3+4) = 144 link classes
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,bool,bool)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,bool,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,bool,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,bool,double)

    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,int,bool)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,int,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,int,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,int,double)

    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,float,bool)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,float,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,float,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,float,double)

    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,double,bool)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,double,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,double,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK,double,double)

// The macro is not needed anymore
#undef REGISTER_VEC_TO_SCALAR_LINK

// Declare a macro to register all the vec to scalar link classes (normalized) using OP_VEC_TO_SCALAR_ALL declared in linkevaluatorid.h
#define REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED(scalarType,vecType,vecDim,vecComponent) \
    registerSerializableType(new VEC_TO_SCALAR_LINK_CLASSNAME_NORMALIZED(scalarType, vecType, vecDim, vecComponent) ());

    // Use the macros to register 3x3x(2+3+4) = 81 link classes
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,int,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,int,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,int,double)

    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,float,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,float,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,float,double)

    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,double,int)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,double,float)
    OP_VEC_TO_SCALAR_ALL(REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED,double,double)

// The macro is not needed anymore
#undef REGISTER_VEC_TO_SCALAR_LINK_NORMALIZED

    // interval
    registerSerializableType(new LinkEvaluatorIntIntervalIdWithLimits());
    registerSerializableType(new LinkEvaluatorFloatIntervalIdWithLimits());

    registerSerializableType(new LinkEvaluatorFloatIntervalMinProjection());
    registerSerializableType(new LinkEvaluatorFloatIntervalMaxProjection());
    registerSerializableType(new LinkEvaluatorIntIntervalMinProjection());
    registerSerializableType(new LinkEvaluatorIntIntervalMaxProjection());

    registerSerializableType(new TemplateLinkEvaluatorFloatIntervalIdentityWithBounds());
    registerSerializableType(new TemplateLinkEvaluatorIntIntervalIdentityWithBounds());


    // bounding box
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxIdWithLimits());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxIdWithLimits());

    registerSerializableType(new LinkEvaluatorFloatBoundingBoxMinXComponent());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxMaxXComponent());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxMinYComponent());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxMaxYComponent());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxMinZComponent());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxMaxZComponent());

    registerSerializableType(new LinkEvaluatorFloatBoundingBoxXComponentInterval());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxYComponentInterval());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxZComponentInterval());

    registerSerializableType(new LinkEvaluatorFloatBoundingBoxXComponentIntervalWithBounds());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxYComponentIntervalWithBounds());
    registerSerializableType(new LinkEvaluatorFloatBoundingBoxZComponentIntervalWithBounds());

    registerSerializableType(new LinkEvaluatorIntBoundingBoxMinXComponent());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxMaxXComponent());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxMinYComponent());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxMaxYComponent());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxMinZComponent());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxMaxZComponent());

    registerSerializableType(new LinkEvaluatorIntBoundingBoxXComponentInterval());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxYComponentInterval());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxZComponentInterval());

    registerSerializableType(new LinkEvaluatorIntBoundingBoxXComponentIntervalWithBounds());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxYComponentIntervalWithBounds());
    registerSerializableType(new LinkEvaluatorIntBoundingBoxZComponentIntervalWithBounds());



    // non-numeric links
    registerSerializableType(new LinkEvaluatorStringId());
    registerSerializableType(new LinkEvaluatorIntStringId());
    registerSerializableType(new LinkEvaluatorFloatStringId());
    registerSerializableType(new LinkEvaluatorDoubleStringId());

    registerSerializableType(new LinkEvaluatorShaderId());

    registerSerializableType(new LinkEvaluatorCameraId());
    registerSerializableType(new LinkEvaluatorCameraOrientationId());
    registerSerializableType(new LinkEvaluatorCameraPosId());
    registerSerializableType(new LinkEvaluatorCameraLookId());
    registerSerializableType(new LinkEvaluatorCameraFocusId());
    registerSerializableType(new LinkEvaluatorCameraFrustumId());

    registerSerializableType(new LinkEvaluatorTransFunc1DKeysId());
    registerSerializableType(new LinkEvaluatorTransFunc1DGaussianId());
    registerSerializableType(new LinkEvaluatorTransFunc2DPrimitivesId());
    registerSerializableType(new LinkEvaluatorButtonId());


    registerSerializableType(new LinkEvaluatorColorSwitchBool());
    //----------------------------------------------------------------------
    // link evaluators end

    // geometry
    registerSerializableType(new VertexGeometry());
    registerSerializableType(new FaceGeometry());
    registerSerializableType(new MeshGeometry());
    registerSerializableType(new MeshListGeometry());
    registerSerializableType(new GeometrySequence());
    registerSerializableType(new TriangleMeshGeometrySimple());
    registerSerializableType(new TriangleMeshGeometryNormal());
    registerSerializableType(new TriangleMeshGeometryColor());
    registerSerializableType(new TriangleMeshGeometryTexCoord());
    registerSerializableType(new TriangleMeshGeometryColorNormal());
    registerSerializableType(new TriangleMeshGeometryNormalTexCoord());
    registerSerializableType(new TriangleMeshGeometryColorTexCoord());
    registerSerializableType(new TriangleMeshGeometryColorNormalTexCoord());
    registerSerializableType(new GlMeshGeometryUInt16Simple());
    registerSerializableType(new GlMeshGeometryUInt16Color());
    registerSerializableType(new GlMeshGeometryUInt16Normal());
    registerSerializableType(new GlMeshGeometryUInt16TexCoord());
    registerSerializableType(new GlMeshGeometryUInt16ColorNormal());
    registerSerializableType(new GlMeshGeometryUInt16NormalTexCoord());
    registerSerializableType(new GlMeshGeometryUInt16ColorTexCoord());
    registerSerializableType(new GlMeshGeometryUInt16ColorNormalTexCoord());
    registerSerializableType(new GlMeshGeometryUInt32Simple());
    registerSerializableType(new GlMeshGeometryUInt32Color());
    registerSerializableType(new GlMeshGeometryUInt32Normal());
    registerSerializableType(new GlMeshGeometryUInt32TexCoord());
    registerSerializableType(new GlMeshGeometryUInt32ColorNormal());
    registerSerializableType(new GlMeshGeometryUInt32NormalTexCoord());
    registerSerializableType(new GlMeshGeometryUInt32ColorTexCoord());
    registerSerializableType(new GlMeshGeometryUInt32ColorNormalTexCoord());
    registerSerializableType(new TriangleMeshGeometryUInt16IndexedNormal());
    registerSerializableType(new TriangleMeshGeometryUInt32IndexedNormal());
    registerSerializableType(new TriangleMeshGeometryUInt32IndexedColorNormal());
    registerSerializableType(new TriangleMeshGeometryUInt32IndexedNormalTexCoord());
    registerSerializableType(new PointListGeometryVec3());
    registerSerializableType(new PointSegmentListGeometryVec3());

    // volume derived data
    registerSerializableType(new VolumeMinMax());
    registerSerializableType(new VolumeMinMaxMagnitude());
    registerSerializableType(new VolumeHash());
    registerSerializableType(new VolumePreview());
    registerSerializableType(new VolumeHistogramIntensity());
    registerSerializableType(new VolumeHistogramIntensityGradient());

    // meta data
    registerSerializableType(new BoolMetaData());
    registerSerializableType(new StringMetaData());
    registerSerializableType(new IntMetaData());
    registerSerializableType(new SizeTMetaData());
    registerSerializableType(new FloatMetaData());
    registerSerializableType(new DoubleMetaData());

    registerSerializableType(new Vec2MetaData());
    registerSerializableType(new IVec2MetaData());
    registerSerializableType(new Vec3MetaData());
    registerSerializableType(new IVec3MetaData());
    registerSerializableType(new Vec4MetaData());
    registerSerializableType(new IVec4MetaData());
    registerSerializableType(new Mat2MetaData());
    registerSerializableType(new Mat3MetaData());
    registerSerializableType(new Mat4MetaData());
    registerSerializableType(new DateTimeMetaData());
    registerSerializableType(new PositionMetaData());
    registerSerializableType(new RealWorldMappingMetaData());
    registerSerializableType(new WindowStateMetaData());

    registerSerializableType("SerializableVectorMetaData::Processor", new SerializableVectorMetaData<Processor*>());
    registerSerializableType("SerializableVectorMetaData::StringMetaData", new SerializableVectorMetaData<StringMetaData*>());
    registerSerializableType("SerializableVectorMetaData::MetaDataBase", new SerializableVectorMetaData<MetaDataBase*>());

    registerSerializableType(new VolumeOctree());
    registerSerializableType(new OctreeBrickPoolManagerRAM());
    registerSerializableType(new OctreeBrickPoolManagerDisk(64<<20, 512<<20, ""));
    registerSerializableType(new OctreeBrickPoolManagerMmap("", ""));

    // transfer functions
    registerSerializableType("TransFuncIntensity", new TransFunc1DKeys());
    registerSerializableType("TransFuncGaussian", new TransFunc1DGaussian());
    registerSerializableType("TransFuncIntensityGradient", new TransFunc2DPrimitives());
    registerSerializableType(new TransFuncMappingKey(0.f, tgt::vec4(0.f)));
    registerSerializableType(new TransFuncMappingCurve(0.f, 0.01f, 0.f, tgt::vec4(0.f)));
    registerSerializableType(new TransFuncTriangle());
    registerSerializableType(new TransFuncQuad());
    registerSerializableType(new TransFuncBanana());

    // ROI
    registerSerializableType(new ROICube());
    registerSerializableType(new ROICylinder());
    registerSerializableType(new ROISphere());
    registerSerializableType(new ROIRaster());
    registerSerializableType(new ROIUnion());
    registerSerializableType(new ROISubtract());
    registerSerializableType(new ROIGraph());

    // io
    registerVolumeReader(new DatVolumeReader());
    registerVolumeReader(new RawVolumeReader());
    registerVolumeReader(new VvdVolumeReader());
    registerVolumeReader(new VvodVolumeReader());
    registerVolumeWriter(new DatVolumeWriter());
    registerVolumeWriter(new VvdVolumeWriter());
    registerVolumeWriter(new VvodVolumeWriter());

    // animation (TODO: convert resources into VoreenSerializableObjects)
    if (VoreenApplication::app()) {
        std::vector<SerializableFactory*> animationFactories = Animation::getSerializerFactories();
        for (size_t i=0; i<animationFactories.size(); i++)
            VoreenApplication::app()->registerSerializerFactory(animationFactories.at(i));
    }


    // Instance operator for all scalar types:
    INST_SCALAR_TYPES(VolumeOperatorInvert, VolumeOperatorInvertGeneric)
    //INST_VECTOR_TYPES(VolumeOperatorInvert, VolumeOperatorInvertGeneric)

    INST_SCALAR_TYPES(VolumeOperatorMirrorX, VolumeOperatorMirrorXGeneric)
    INST_VECTOR_TYPES(VolumeOperatorMirrorX, VolumeOperatorMirrorXGeneric)

    INST_SCALAR_TYPES(VolumeOperatorMirrorY, VolumeOperatorMirrorYGeneric)
    INST_VECTOR_TYPES(VolumeOperatorMirrorY, VolumeOperatorMirrorYGeneric)

    INST_SCALAR_TYPES(VolumeOperatorMirrorZ, VolumeOperatorMirrorZGeneric)
    INST_VECTOR_TYPES(VolumeOperatorMirrorZ, VolumeOperatorMirrorZGeneric)

    INST_SCALAR_TYPES(VolumeOperatorTranspose, VolumeOperatorTransposeGeneric)
    INST_VECTOR_TYPES(VolumeOperatorTranspose, VolumeOperatorTransposeGeneric)
    INST_TENSOR_TYPES(VolumeOperatorTranspose, VolumeOperatorTransposeGeneric)

    INST_SCALAR_TYPES(VolumeOperatorSwapEndianness, VolumeOperatorSwapEndiannessGeneric)
    INST_VECTOR_TYPES(VolumeOperatorSwapEndianness, VolumeOperatorSwapEndiannessGeneric)
    INST_TENSOR_TYPES(VolumeOperatorSwapEndianness, VolumeOperatorSwapEndiannessGeneric)

    INST_SCALAR_TYPES(VolumeOperatorMedian, VolumeOperatorMedianGeneric)

    INST_SCALAR_TYPES(VolumeOperatorGaussian, VolumeOperatorGaussianGeneric)

    INST_SCALAR_TYPES(VolumeOperatorCubeDilation, VolumeOperatorCubeDilationGeneric)
    INST_SCALAR_TYPES(VolumeOperatorSphereDilation, VolumeOperatorSphereDilationGeneric)
    INST_SCALAR_TYPES(VolumeOperatorCubeErosion, VolumeOperatorCubeErosionGeneric)
    INST_SCALAR_TYPES(VolumeOperatorSphereErosion, VolumeOperatorSphereErosionGeneric)

    INST_SCALAR_TYPES(VolumeOperatorResample, VolumeOperatorResampleGeneric)
    INST_VECTOR_TYPES(VolumeOperatorResample, VolumeOperatorResampleGeneric)

    INST_SCALAR_TYPES(VolumeOperatorHalfsample, VolumeOperatorHalfsampleGeneric)
    INST_VECTOR_TYPES(VolumeOperatorHalfsample, VolumeOperatorHalfsampleGeneric)

    INST_SCALAR_TYPES(VolumeOperatorRegionGrow, VolumeOperatorRegionGrowGeneric)

    INST_SCALAR_TYPES(VolumeOperatorResize, VolumeOperatorResizeGeneric)
    INST_VECTOR_TYPES(VolumeOperatorResize, VolumeOperatorResizeGeneric)
    INST_TENSOR_TYPES(VolumeOperatorResize, VolumeOperatorResizeGeneric)

    INST_SCALAR_TYPES(VolumeOperatorResizePowerOfTwo, VolumeOperatorResizePowerOfTwoGeneric)
    INST_VECTOR_TYPES(VolumeOperatorResizePowerOfTwo, VolumeOperatorResizePowerOfTwoGeneric)

    INST_SCALAR_TYPES(VolumeOperatorSubset, VolumeOperatorSubsetGeneric)
    INST_VECTOR_TYPES(VolumeOperatorSubset, VolumeOperatorSubsetGeneric)

    INST_SCALAR_TYPES(VolumeOperatorIsUniform, VolumeOperatorIsUniformGeneric)
    INST_VECTOR_TYPES(VolumeOperatorIsUniform, VolumeOperatorIsUniformGeneric)

    INST_SCALAR_TYPES(VolumeOperatorNumSignificant, VolumeOperatorNumSignificantGeneric)
    INST_VECTOR_TYPES(VolumeOperatorNumSignificant, VolumeOperatorNumSignificantGeneric)

    INST_SCALAR_TYPES(VolumeOperatorCalcError, VolumeOperatorCalcErrorGeneric)
    INST_VECTOR_TYPES(VolumeOperatorCalcError, VolumeOperatorCalcErrorGeneric)

    INST_SCALAR_TYPES(VolumeOperatorEqualize, VolumeOperatorEqualizeGeneric)
    //INST_VECTOR_TYPES(VolumeOperatorEqualize, VolumeOperatorEqualizeGeneric)

    INST_SCALAR_TYPES(VolumeOperatorUncertaintyMeasure, VolumeOperatorUncertaintyMeasureGeneric)
    INST_VECTOR_TYPES(VolumeOperatorUncertaintyMeasure, VolumeOperatorUncertaintyMeasureGeneric)

    INST_SCALAR_TYPES(VolumeOperatorConnectedComponentAnalysis, VolumeOperatorConnectedComponentAnalysisGeneric)
}

CoreModule::~CoreModule() {
    Animation::deleteSerializerFactories();
}

} // namespace
