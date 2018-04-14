#ifndef SURFACESUBDIVISIONRECONSTRUCTION_H
#define SURFACESUBDIVISIONRECONSTRUCTION_H

/**
 *   The relatively simplest methods for Vascular Reconstruction.
 *   Given Vascular Tree, we generate a simplified mesh & then
 *   apply Mesh Subdivision Methods. The main issue with these kind of
 *   approaches are:
 *       1) In the case of, triangular degeneration... it may be possible that degenerate triangles are
 *          developed. In general we can then apply a decimation method to reduce the number of triangles
 *          & maintain the general method
 *       2) Diminished Radius.. Each subdivision method, in the case of Catmull-Clark subdivision for example
 *          will reduce the radius.
 *
 *    Here we provide different methods to grab the Vascular Tree as defined,
 *    Generate a Mesh from the Vascular Tree.. either as a Simple QuadMesh, or a more Detailed one
 *    And apply subdivision steps.
 *
 *    The possible subdivision schemes are:
 *        * vtkLoopSubdivisionFilter
 *        * vtkLinearSubdivisionFilter
 *        * vtkAdaptiveSubdivisionFilter
 *        * vtkButterflySubdivisionFilter
 *        * Catmull-Clark SubdivisionFilter
 * */


#ifndef __TUBNAVCONTAINER_H
    #include "tubNavContainer.h"
#endif
#include "vtkCatmullClarkFilter.h"

#include <vtkLoopSubdivisionFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkAdaptiveSubdivisionFilter.h>
#include <vtkButterflySubdivisionFilter.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

class SurfaceSubdivisionReconstruction
{
public:
    SurfaceSubdivisionReconstruction();

    vtkSmartPointer<vtkPolyData> GetAllSimpleQuadReconstruction();
    vtkSmartPointer<vtkPolyData> GetReconstructionFromPointCloud();

    void SetVascularTree(std::shared_ptr<tubNav::Container<double>> tree);

private:
    vtkSmartPointer<vtkPolyData> BreakDaughter(vtkSmartPointer<vtkPolyData> parentMesh, vtkSmartPointer<vtkPolyData> daughterMesh, int indexOfPath);

    // Given that a smoothing step is done, we can artificially increase the size of the path by
    // setting percent, the regular one is 1.0
    vtkSmartPointer<vtkPolyData>  GetSinglePathReconstruction(int id, int numSlices, float percent = 1.0);

    // Basically the idea is to create a point cloud... and then apply surface reconstruction from it
    vtkSmartPointer<vtkPoints> GetSinglePathPoints(int id, int numSlices, float percent = 1.0);

    std::shared_ptr<tubNav::Container<double>> vascTree;

    bool ArePointsInPath(int pathToCheck, vtkSmartPointer<vtkPoints> pointsToCheck);

};

#endif // SURFACESUBDIVISIONRECONSTRUCTION_H
