#ifndef VASCRESC_LIB
#define VASCRESC_LIB


#include "ReconstructionLibrary/Core/tubNavContainer.h"
#include "ReconstructionLibrary/MathCore/FiniteDiff.h"
#include "ReconstructionLibrary/MathCore/CoordinateFrame.h"
#include "ReconstructionLibrary/IO/General.h"
#include "ReconstructionLibrary/MathCore/ProjectionOperations.h"
#include "ReconstructionLibrary/Reconstruct2D/CPR.h"
#include "ReconstructionLibrary/Reconstruct2D/Unfolding.h"
#include "ReconstructionLibrary/Reconstruct2D/cprinfo.h"
#include "ReconstructionLibrary/Reconstruct3D/surfacereconstruction.h"

#endif


#ifndef SEMINAR_LIB
#define SEMINAR_LIB

#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageConstIteratorWithIndex.h"

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include <vector>
typedef double T;

//The velocity map is formed by a vector representing the velocity (vx,vy,vz) and time
typedef itk::Image< itk::CovariantVector<T, 4 >,  3u> velocityMap;
typedef itk::Image< itk::CovariantVector<T, 9 >,  3u> TensorMapType;
typedef itk::Image< T,  3u> WallShearStressMap;


typedef signed short InputPixelType;
const unsigned int InputDimension = 3;
typedef itk::Image< InputPixelType, InputDimension > ImageType;

#endif
