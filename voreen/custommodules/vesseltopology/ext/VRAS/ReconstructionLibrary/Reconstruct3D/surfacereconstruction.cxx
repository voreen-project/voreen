#include "surfacereconstruction.h"
#include <vtkFillHolesFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkClipClosedSurface.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkInformation.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkPlaneCollection.h>
#include <vtkPlane.h>
#include <vtkPolyDataConnectivityFilter.h>
#include "vtkCatmullClarkFilter.h"


SurfaceReconstruction::SurfaceReconstruction(std::shared_ptr<tubNav::Container<double> > tree)
{
    vascTree = tree;
    method = ReconstructionMethod::Subdivision;
    currentReconstruction = nullptr;
    smoothingFilterToApply = nullptr;
    decimationPercentage = 0.5;
    enableDecimation = true;

}

SurfaceReconstruction::SurfaceReconstruction(){
    vascTree = nullptr;
    method = ReconstructionMethod::Subdivision;
    currentReconstruction = nullptr;
    smoothingFilterToApply = nullptr;
    decimationPercentage = 0.5;
    enableDecimation = true;
}


void SurfaceReconstruction::SetVascularTree(std::shared_ptr<tubNav::Container<double> > tree){
    vascTree = tree;
    it.SetContainer(vascTree);
    it.InitTraversal();
    frame = it.GetCoordinateFrame();
    frame->InitTraversal();
}

void SurfaceReconstruction::ChangeReconstructionMethod(ReconstructionMethod newMethod){
    method = newMethod;
}

void SurfaceReconstruction::SetReconstruction(vtkSmartPointer<vtkPolyData> customReconstruction)
{
    method = ReconstructionMethod::Custom;
    _customReconstruction = customReconstruction;
    std::cout << "Initial " << customReconstruction->GetNumberOfPolys() << std::endl;
    currentReconstruction = customReconstruction;
}

void SurfaceReconstruction::SetSmoothingFilter(SmoothingFilter filter, int numberOfSubdivisions){

    if ( filter == SmoothingFilter::LinearSubdivision){
        smoothingFilterToApply = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
        dynamic_cast<vtkLinearSubdivisionFilter *> (smoothingFilterToApply.GetPointer())->SetNumberOfSubdivisions(numberOfSubdivisions);
    }

    if ( filter == SmoothingFilter::LoopSubdivision){
        smoothingFilterToApply = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
        dynamic_cast<vtkLoopSubdivisionFilter *> (smoothingFilterToApply.GetPointer())->SetNumberOfSubdivisions(numberOfSubdivisions);
    }

    if ( filter == SmoothingFilter::ButterflySubdivision){
        smoothingFilterToApply = vtkSmartPointer<vtkButterflySubdivisionFilter>::New();
        dynamic_cast<vtkButterflySubdivisionFilter *> (smoothingFilterToApply.GetPointer())->SetNumberOfSubdivisions(numberOfSubdivisions);
    }

    if ( filter == SmoothingFilter::LaplacianSubdivision){
        smoothingFilterToApply = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        dynamic_cast<vtkSmoothPolyDataFilter *> (smoothingFilterToApply.GetPointer())->SetNumberOfIterations(numberOfSubdivisions);
    }


}


void SurfaceReconstruction::ResetReconstruction(){

    if (method == ReconstructionMethod::Custom){
           currentReconstruction = _customReconstruction;


           ApplyAllClippings();

           //now apply all the after filters ...
           if ( smoothingFilterToApply != nullptr){
               smoothingFilterToApply->SetInputData(currentReconstruction);
               smoothingFilterToApply->Update();
               currentReconstruction = smoothingFilterToApply->GetOutput();
           }

           if (enableDecimation){
               vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
               decimate->SetTargetReduction(decimationPercentage);
               decimate->SetInputData(currentReconstruction);
               decimate->Update();
               currentReconstruction = decimate->GetOutput();
           }




    }
    else {

        std::cout  << "Applying reconstruction ? " << std::endl;
        SurfaceSubdivisionReconstruction reconstruct;
        reconstruct.SetVascularTree(vascTree);
        auto vec = reconstruct.GetAllSimpleQuadReconstruction();
        currentReconstruction =  vec;

        if ( smoothingFilterToApply != nullptr){
            smoothingFilterToApply->SetInputData(currentReconstruction);
            smoothingFilterToApply->Update();
            currentReconstruction = smoothingFilterToApply->GetOutput();
        }


    }
}

vtkSmartPointer<vtkPolyData> SurfaceReconstruction::GetReconstruction(){
       return currentReconstruction;
}

std::shared_ptr<tubNav::CoordinateFrame<double>> SurfaceReconstruction::GetCurrentFrame(){
    return frame;
}


void SurfaceReconstruction::MoveToNextPoint(){
    frame->Next();
    ++it;
    if ( frame->IsAtEnd()){
        frame = it.GetCoordinateFrame();
        frame->InitTraversal();
    }
}

void SurfaceReconstruction::AddNewClipping(tubNav::Point<double> location, tubNav::Point<double> n, bool clipDescendants){
   // adds new clipping at current location...
   Clipping clip;
   clip.location = location;
   clip.normal = n;
   clip.clipDescendants = clipDescendants;
   clippings.push_back(clip);
}

void SurfaceReconstruction::ApplyAllClippings(){
    for(unsigned int i = 0; i < clippings.size();i++)
        ApplyClipping(i);
}

void SurfaceReconstruction::ApplyClipping(int idx){
    Clipping clip = GetClipping(idx);

    std::cout << "Applying clip " << idx << std::endl;
    // Idea No.1 find the radius at the location...

    int pathId = -1;
    tubNav::ContainerIterator<double> clippingIt(vascTree);
    // Assumption... we are only removing from the current path

    clippingIt.InitTraversal();
    auto cframe = clippingIt.GetCoordinateFrame();
    cframe->InitTraversal();

    float  minDistance = 9999999;
    std::cout << "Starting travelsal ? " << clip.location << std::endl;
    while(!clippingIt.IsAtEnd()){
        if ( cframe->IsAtEnd()){
            cframe = clippingIt.GetCoordinateFrame();
            cframe->InitTraversal();
        }

        auto center = cframe->GetCenterPoint();
        double d =  center.distance(clip.location,sqrt);
        if ( d < 0.001){
            pathId = clippingIt.GetCurrentPathId();
            break;
        }

        if ( d < minDistance)
            minDistance = d;

        cframe->Next();
        ++clippingIt;

        std::cout << "." ;
    }
    if ( pathId == -1){
        std::cout << "Couldn't find path with point " << minDistance << std::endl;
        return;
    }
    // Find the cells to be extracted....

    // 1. What cells belong as the descendants or ancestors of the
    //    selection point... first we get the generating points of it....

    // if we are removing descendants... then we move until the end of the path....
    std::cout << "Getting inner locations " << std::endl;
    std::vector< decltype(cframe->GetCenterPoint() )> inner;

    if ( clip.clipDescendants){
        // Let's get all the descendants of the path....
        while(!clippingIt.IsAtEnd() && clippingIt.GetCurrentPathId() == pathId){
            inner.push_back( cframe->GetCenterPoint());
            cframe->Next();
            ++clippingIt;
            if ( cframe->IsAtEnd()){
                cframe = clippingIt.GetCoordinateFrame();
                cframe->InitTraversal();
            }
        }
    }
    // if we are removing ancestors ... then we move until the start of the path....
    else {
        // we have to move in the opposit direction

        while(!clippingIt.IsAtStart() && clippingIt.GetCurrentPathId() == pathId){
            inner.push_back(cframe->GetCenterPoint());
            cframe->Back();
            --clippingIt;

            if (cframe->IsAtStart()){ //?
                cframe = clippingIt.GetCoordinateFrame();
                cframe->SetToEnd();
            }

        }
    }
    // 2. Now we look at the polydata and check whether they are withing the sphere
    //    of influence of the inner points...
    std::cout << "Total " << inner.size() << std::endl;

    vtkSmartPointer<vtkCellArray> originalPolys = currentReconstruction->GetPolys();
    vtkSmartPointer<vtkPoints> points = currentReconstruction->GetPoints();
    int numberOfCells = originalPolys->GetNumberOfCells();
    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);

    std::cout << "Number of cells " << numberOfCells << "number of points " << points->GetNumberOfPoints() << std::endl;

    // Two extractions might be necessary, the first one a working area
    // basically we generate a rectangle from the inner max/min dimensions
    // extract the selection based on the points that exist there
    // and then refine the process


    int numberOfPoints = points->GetNumberOfPoints();

    double p[3];

    for(int j = 0; j < numberOfPoints ; j++){
        points->GetPoint(j, p);
        decltype(clip.location) point(p);
        for(unsigned int h =0; h < inner.size(); h++){
            auto currentPoint = inner.at(h);
            double d = currentPoint.distance(point, sqrt);
            if ( d < currentPoint.getR()*1.5 ){ // add 5% to the radius, to make sure it is inside close
                ids->InsertNextValue(j);
                break;
            }
        }
    }
    // 3. We extract them from the indices
     vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
     selectionNode->SetFieldType(vtkSelectionNode::POINT);
     selectionNode->SetContentType(vtkSelectionNode::INDICES);
     selectionNode->SetSelectionList(ids);
     selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(),1);

     vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
     selection->AddNode(selectionNode);

     vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
     extractSelection->SetInputData(0,currentReconstruction);
     extractSelection->SetInputData(1, selection);
     extractSelection->Update();

     vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
     selected->ShallowCopy(extractSelection->GetOutput());
      vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
      geometryFilter->SetInputData(selected);
      geometryFilter->Update();

     // 4. Close the resulting surface
       auto subarea = geometryFilter->GetOutput();
      vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =    vtkSmartPointer<vtkFillHolesFilter>::New();
      fillHolesFilter->SetInputData(subarea);
      fillHolesFilter->SetHoleSize(1000.0);

      // Make the triangle windong order consistent
      vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
      normals->SetInputConnection(fillHolesFilter->GetOutputPort());
      normals->ConsistencyOn();
      normals->SplittingOff();
      normals->Update();

      // Restore the original normals
      normals->GetOutput()->GetPointData()->SetNormals(subarea->GetPointData()->GetNormals());
      auto cleanSubArea = normals->GetOutput();
     // 5. Apply the clip-closed surface on the extraction
     vtkSmartPointer<vtkPlaneCollection> collection = vtkSmartPointer<vtkPlaneCollection>::New();
     vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
     plane->SetOrigin(clip.location.getX(), clip.location.getY(), clip.location.getZ());

     if ( clip.clipDescendants){ // if we remove descendants, we flip the normal of the location
         plane->SetNormal(clip.normal.getX(), clip.normal.getY(), clip.normal.getZ());
     }
     else {
         plane->SetNormal( -clip.normal.getX(), -clip.normal.getY(),-clip.normal.getZ());
     }

    collection->AddItem(plane);
    vtkSmartPointer<vtkClipClosedSurface> clipper = vtkSmartPointer<vtkClipClosedSurface>::New();
    clipper->SetInputData(cleanSubArea);
    clipper->SetClippingPlanes(collection);
    // a scalar value of "0" indicates an original cell, "1" indicates a new cell on a cut face
    // "2" indicates a new cell on the active plane as set by the set active plane
    clipper->SetScalarModeToLabels();
    clipper->Update();


    // This is the area that we want to keep, of the current path....
    auto toRemove = clipper->GetOutput();
    // 6. now from the original, we find the overlapping cells, and remove them
    // we re-do the selection but now in a different fashion

    vtkSmartPointer<vtkIdTypeArray> ids2 = vtkSmartPointer<vtkIdTypeArray>::New();
    ids2->SetNumberOfComponents(1);

    vtkSmartPointer<vtkPoints> removePoints = toRemove->GetPoints();
    double q[3];

    double epsilon = 0.001;
    for(int j = 0; j < numberOfPoints ; j++){
        points->GetPoint(j, p);

        for(int k =0; k < removePoints->GetNumberOfPoints(); k++){
            removePoints->GetPoint(k, q);

            double distance = 0;
            for(int i = 0; i < 3; i++)
                distance += (p[i]-q[i])* (p[i]-q[i]);
            distance = sqrt(distance);

            if ( distance < epsilon)
            {
                ids2->InsertNextValue(j);
                break;
            }
        }
    }

    std::cout << "Num points " << ids2->GetNumberOfValues() << std::endl;

    //7. Now that we have the points that need to be removed... get the selection
    selectionNode->SetFieldType(vtkSelectionNode::POINT);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids2);
    selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(),1);
    selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1);
    extractSelection->Update();

    vtkSmartPointer<vtkUnstructuredGrid> unselected = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unselected->ShallowCopy(extractSelection->GetOutput());
    geometryFilter->SetInputData(unselected);
    geometryFilter->Update();

    fillHolesFilter->SetInputData(geometryFilter->GetOutput());
    fillHolesFilter->SetHoleSize(1000.0);
    fillHolesFilter->Update();

    // Make the triangle windowing order consistent


    // 5) extract the region closest to the specified point
    //    where the point is the selected plane so we can use that...

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivity = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivity->AddInputData(fillHolesFilter->GetOutput());
    connectivity->SetExtractionModeToClosestPointRegion();
    connectivity->SetClosestPoint(clip.location.getX(), clip.location.getY(), clip.location.getZ());
    //connectivity->SetExtractionModeToLargestRegion();
    connectivity->Update();

    auto largest = connectivity->GetOutput();

    currentReconstruction = largest;
}

void SurfaceReconstruction::RemoveClipping(int idx){

}

Clipping SurfaceReconstruction::GetClipping(int idx){

    return clippings.at(idx);
}


void SurfaceReconstruction::MoveToPreviousPoint(){

    // This one I should check how to go back...
    frame->Back();
    --it;

    // I should still figure this out?

}

void SurfaceReconstruction::MoveToStart(){
    it.InitTraversal();
    frame = it.GetCoordinateFrame();
    frame->InitTraversal();
}

bool SurfaceReconstruction::IsAtEnd(){
    return it.IsAtEnd();
}
