#include "SurfaceSubdivision.h"


#include <vtkCleanPolyData.h>
#include <vtkTriangle.h>
#include <vtkAppendFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkDelaunay3D.h>
#include <vector>
#include <vtkGeometryFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include "../IO/General.h"

SurfaceSubdivisionReconstruction::SurfaceSubdivisionReconstruction()
{
   vascTree = nullptr;
}

void SurfaceSubdivisionReconstruction::SetVascularTree(std::shared_ptr<tubNav::Container<double> > tree){
    vascTree = tree;
}


/**
   First step is the reconstruction of the Vascular Tree
   As D. Meyers et al. Surfaces from Contours. ACM Trans. on Graphics
   There are four fundamental subproblems to deal with when doing surface reconstruction

      1. The corresponce problem. Which contours in one cross-section should be connected
         to which contours in other cross-section.
            -> We do not allow overlapping between the cross-sections, this could be remedied by
               using instead non-linear unfolding... we'll try it out

      2. The tiling problem. How should the pairs of given contours be connected? Which vertices
         and edges should form the triangles.
             -> The minimizing rotatin frame takes care of this, allowing us to connect the
                vertices in continuous locations

      3. The branching problem. How to tile the bifurcations, i.e. cross-sections with a different
         number of contours

      4. Surface-fitting problem. What does the precise geometry look-like? A possible post-processing
         step smoothing the mesh.


*/


vtkSmartPointer<vtkPoints> SurfaceSubdivisionReconstruction::GetSinglePathPoints(int id, int numSlices, float percent){
    // Basically we have defined the paths already, so we construct a generalized cylinder approach.
    // The idea is basically to connect the paths
    tubNav::ContainerIterator<double> it(vascTree);
    it.InitTraversal();
    while ( it.GetCurrentPathId() != id ){
        it.MoveNextPath();
    }
    // We have now found the path that we wanted to generate
    auto frame = it.GetCoordinateFrame();

    std::cout << "Size " << frame->GetSize() << std::endl;
    frame->InitTraversal();
     // Create a polydata to store everything in
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // We can generate the points per each two continuous cross-sections
    // and then call clean poly data

    double stepSize = 360.0/numSlices;

    while (!it.IsAtEnd() && it.GetCurrentPathId() == id ){

        //it.IsPointInsideOtherSegment()
        // Now frame has the points and the generating part... let's see if this works
        float r1 = frame->GetCenterPoint().getR(); // Assume we always left in a correct frame
        auto center = frame->GetCenterPoint();
        auto tangent = frame->GetTangentDirection();


        bool fulfills = false;
        std::vector<decltype(center)> firstRow;

        while(!fulfills){

            firstRow.clear();

            for(int i = 0; i < numSlices+1; i++){
                auto pt = frame->GetPoint(i*stepSize, r1*percent);
                firstRow.push_back(pt);
            }
            center = frame->GetCenterPoint();
            tangent = frame->GetTangentDirection();


            if (true){ //!it.IsPointInsideOtherSegment()){
                fulfills = true;
            }
            else {
                frame->Next();
                ++it;

                if ( frame->IsAtEnd()){
                    frame = it.GetCoordinateFrame();
                    frame->InitTraversal();
                }
            }
            if ( it.GetCurrentPathId() !=id) fulfills = true;
        }
         if ( it.GetCurrentPathId() !=id) break;

         frame->Next();
         ++it;

        // íf it is the first row ever, then we generate the closing
        if ( points->GetNumberOfPoints() == 0 || (frame->IsAtEnd() && it.GetNextPathId() != id) ){

            points->InsertNextPoint( center.getX(), center.getY(),center.getZ());

            for(int i = 0; i < numSlices+1; i++){
                auto pt0 = firstRow.at(i);
                // This are the points in the direction to angle i, lets just add points
                auto dir = pt0 - center;
                // int id1 =
                dir.normalize(sqrt);

                double distance = sqrt(dir.normSquared());
                double step = distance*0.02;
                double start = step;

                while( start <= distance){
                    auto newPt = center +  dir*start;
                    points->InsertNextPoint( newPt.getX(), newPt.getY(), newPt.getZ());
                    start += step;
                }
            }
        }



        if ( frame->IsAtEnd()){
            frame = it.GetCoordinateFrame();
            frame->InitTraversal();
        }

        if ( it.GetCurrentPathId() !=id) break;

        std::vector<decltype(center)> secondRow;
        decltype(center) lastCenter;

        fulfills = false;

        while(!fulfills){

            secondRow.clear();

            int pointsBelow = 0;

            float r2 = frame->GetCenterPoint().getR();
            for(int i = 0; i < numSlices+1; i++){
                auto pt = frame->GetPoint(i*stepSize, r2*percent);
                if ( tubNav::IsPointBelow(tangent, center, pt) ) pointsBelow++;
                secondRow.push_back(pt);
            }

            if ( pointsBelow == 0 ){ // && !it.IsPointInsideOtherSegment()){
                fulfills = true;
                // if it doesn't full fill, then it's not drawn
                lastCenter = frame->GetCenterPoint();
            }
            else {
                frame->Next();
                ++it;

                if ( frame->IsAtEnd()){
                    frame = it.GetCoordinateFrame();
                    frame->InitTraversal();
                }
            }
            if ( it.GetCurrentPathId() !=id) fulfills = true;
        }

        if ( it.GetCurrentPathId() !=id) {

            points->InsertNextPoint( center.getX(), center.getY(),center.getZ());
            for(int i = 0; i < numSlices; i++){

                auto pt0 = firstRow.at(i);
                // This are the points in the direction to angle i, lets just add points
                auto dir = pt0 - center;
                // int id1 =
                dir.normalize(sqrt);

                double distance = sqrt(dir.normSquared());
                double step = distance*0.02;
                double start = step;

                while( start <= distance){
                    auto newPt = center +  dir*start;
                    points->InsertNextPoint( newPt.getX(), newPt.getY(), newPt.getZ());
                    start += step;
                }
            }
            break;
        }

        for(int i = 0; i < numSlices; i++){
            auto pt0 = firstRow.at(i);
            auto pt1 = firstRow.at(i+1);
            auto pt2 = secondRow.at(i);
            auto pt3 = secondRow.at(i+1);

            //Now that we have the actual location, we add the two triangles
            points->InsertNextPoint( pt0.getX(), pt0.getY(), pt0.getZ());
            points->InsertNextPoint( pt1.getX(), pt1.getY(), pt1.getZ());
            points->InsertNextPoint( pt2.getX(), pt2.getY(), pt2.getZ());
            points->InsertNextPoint( pt3.getX(), pt3.getY(), pt3.getZ());
        }

        // fullFills and didn't break
        if (  !(!it.IsAtEnd() && it.GetCurrentPathId() == id) ){

            std::cout << "Actually reaches this point?" << std::endl;
            int id0 = points->InsertNextPoint( lastCenter.getX(), lastCenter.getY(), lastCenter.getZ());
            for(int i = 0; i < numSlices; i++){
                auto pt0 = secondRow.at(i);
                auto pt1 = secondRow.at(i+1);
                int id1 = points->InsertNextPoint( pt0.getX(), pt0.getY(), pt0.getZ());
                int id2 = points->InsertNextPoint( pt1.getX(), pt1.getY(), pt1.getZ());

                vtkSmartPointer<vtkTriangle> triangle1 =  vtkSmartPointer<vtkTriangle>::New();

                triangle1->GetPointIds()->SetId(0,id0);
                triangle1->GetPointIds()->SetId(1,id2);
                triangle1->GetPointIds()->SetId(2,id1);
            }
        }
    }

    return points;
}

vtkSmartPointer<vtkPolyData>  SurfaceSubdivisionReconstruction::GetSinglePathReconstruction(int id, int numSlices, float percent){
   // Basically we have defined the paths already, so we construct a generalized cylinder approach.
   // The idea is basically to connect the paths
   tubNav::ContainerIterator<double> it(vascTree);
   it.InitTraversal();
   while ( it.GetCurrentPathId() != id ){
       it.MoveNextPath();
   }
   // We have now found the path that we wanted to generate
   auto frame = it.GetCoordinateFrame();
   frame->InitTraversal();
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polydata =  vtkSmartPointer<vtkPolyData>::New();

    // Add the points and quads to the dataset
   vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
   // We can generate the points per each two continuous cross-sections
   // and then call clean poly data

   double stepSize = 360.0/numSlices;


   while (!it.IsAtEnd() && it.GetCurrentPathId() == id ){

       //it.IsPointInsideOtherSegment()
       // Now frame has the points and the generating part... let's see if this works
       float r1 = frame->GetCenterPoint().getR(); // Assume we always left in a correct frame
       auto center = frame->GetCenterPoint();
       auto tangent = frame->GetTangentDirection();


       bool fulfills = false;
       std::vector<decltype(center)> firstRow;

       while(!fulfills){

           firstRow.clear();

           for(int i = 0; i < numSlices+1; i++){
               auto pt = frame->GetPoint(i*stepSize, r1*percent);
               firstRow.push_back(pt);
           }
           center = frame->GetCenterPoint();
           tangent = frame->GetTangentDirection();


           if (true){ //!it.IsPointInsideOtherSegment()){
               fulfills = true;
           }
           else {
               frame->Next();
               ++it;

               if ( frame->IsAtEnd()){
                   frame = it.GetCoordinateFrame();
                   frame->InitTraversal();
               }
           }
           if ( it.GetCurrentPathId() !=id) fulfills = true;
       }
        if ( it.GetCurrentPathId() !=id) break;

        frame->Next();
        ++it;

       // íf it is the first row ever, then we generate the closing
       if ( quads->GetNumberOfCells() == 0 || (frame->IsAtEnd() && it.GetNextPathId() != id) ){

           int id0 = points->InsertNextPoint( center.getX(), center.getY(),center.getZ());
           for(int i = 0; i < numSlices; i++){
               auto pt0 = firstRow.at(i);
               auto pt1 = firstRow.at(i+1);
               int id1 = points->InsertNextPoint( pt0.getX(), pt0.getY(), pt0.getZ());
               int id2 = points->InsertNextPoint( pt1.getX(), pt1.getY(), pt1.getZ());

               vtkSmartPointer<vtkTriangle> triangle1 =  vtkSmartPointer<vtkTriangle>::New();

               triangle1->GetPointIds()->SetId(0,id0);
               triangle1->GetPointIds()->SetId(1,id2);
               triangle1->GetPointIds()->SetId(2,id1);
               quads->InsertNextCell(triangle1);
           }
       }



       if ( frame->IsAtEnd()){
           frame = it.GetCoordinateFrame();
           frame->InitTraversal();
       }

       if ( it.GetCurrentPathId() !=id) break;

       std::vector<decltype(center)> secondRow;
       decltype(center) lastCenter;

       fulfills = false;

       while(!fulfills){

           secondRow.clear();

           int pointsBelow = 0;

           float r2 = frame->GetCenterPoint().getR();
           for(int i = 0; i < numSlices+1; i++){
               auto pt = frame->GetPoint(i*stepSize, r2*percent);
               if ( tubNav::IsPointBelow(tangent, center, pt) ) pointsBelow++;
               secondRow.push_back(pt);
           }

           if ( pointsBelow == 0 ){ // && !it.IsPointInsideOtherSegment()){
               fulfills = true;
               // if it doesn't full fill, then it's not drawn
               lastCenter = frame->GetCenterPoint();
           }
           else {              
               frame->Next();
               ++it;

               if ( frame->IsAtEnd()){
                   frame = it.GetCoordinateFrame();
                   frame->InitTraversal();
               }
           }
           if ( it.GetCurrentPathId() !=id) fulfills = true;
       }

       if ( it.GetCurrentPathId() !=id) {

           int id0 = points->InsertNextPoint( center.getX(), center.getY(),center.getZ());
           for(int i = 0; i < numSlices; i++){
               auto pt0 = firstRow.at(i);
               auto pt1 = firstRow.at(i+1);
               int id1 = points->InsertNextPoint( pt0.getX(), pt0.getY(), pt0.getZ());
               int id2 = points->InsertNextPoint( pt1.getX(), pt1.getY(), pt1.getZ());

               vtkSmartPointer<vtkTriangle> triangle1 =  vtkSmartPointer<vtkTriangle>::New();

               triangle1->GetPointIds()->SetId(0,id0);
               triangle1->GetPointIds()->SetId(1,id2);
               triangle1->GetPointIds()->SetId(2,id1);
               quads->InsertNextCell(triangle1);
           }


           break;
       }

       for(int i = 0; i < numSlices; i++){
           auto pt0 = firstRow.at(i);
           auto pt1 = firstRow.at(i+1);
           auto pt2 = secondRow.at(i);
           auto pt3 = secondRow.at(i+1);

           //Now that we have the actual location, we add the two triangles
           int id0 = points->InsertNextPoint( pt0.getX(), pt0.getY(), pt0.getZ());
           int id1 = points->InsertNextPoint( pt1.getX(), pt1.getY(), pt1.getZ());
           int id2 = points->InsertNextPoint( pt2.getX(), pt2.getY(), pt2.getZ());
           int id3 = points->InsertNextPoint( pt3.getX(), pt3.getY(), pt3.getZ());

           // Now that we have the ids and the points we add the two triangles
           vtkSmartPointer<vtkTriangle> triangle1 =  vtkSmartPointer<vtkTriangle>::New();
           vtkSmartPointer<vtkTriangle> triangle2 =  vtkSmartPointer<vtkTriangle>::New();


           triangle1->GetPointIds()->SetId(0,id0);
           triangle1->GetPointIds()->SetId(1,id1);
           triangle1->GetPointIds()->SetId(2,id3);

           triangle2->GetPointIds()->SetId(0,id0);
           triangle2->GetPointIds()->SetId(1,id3);
           triangle2->GetPointIds()->SetId(2,id2);

           quads->InsertNextCell(triangle1);
           quads->InsertNextCell(triangle2);
       }

       // fullFills and didn't break
       if (  !(!it.IsAtEnd() && it.GetCurrentPathId() == id) ){
           std::cout << "Actually reaches this point?" << std::endl;

           int id0 = points->InsertNextPoint( lastCenter.getX(), lastCenter.getY(), lastCenter.getZ());
           for(int i = 0; i < numSlices; i++){
               auto pt0 = secondRow.at(i);
               auto pt1 = secondRow.at(i+1);
               int id1 = points->InsertNextPoint( pt0.getX(), pt0.getY(), pt0.getZ());
               int id2 = points->InsertNextPoint( pt1.getX(), pt1.getY(), pt1.getZ());

               vtkSmartPointer<vtkTriangle> triangle1 =  vtkSmartPointer<vtkTriangle>::New();

               triangle1->GetPointIds()->SetId(0,id0);
               triangle1->GetPointIds()->SetId(1,id2);
               triangle1->GetPointIds()->SetId(2,id1);
               quads->InsertNextCell(triangle1);
           }
       }
   }

   polydata->SetPoints(points);
   polydata->SetPolys(quads);


   vtkSmartPointer<vtkCleanPolyData> cleanPolyData =  vtkSmartPointer<vtkCleanPolyData>::New();
   cleanPolyData->SetInputData(polydata);
   cleanPolyData->Update();


   vtkSmartPointer<vtkPolyData> fnlResult = cleanPolyData->GetOutput();

   return fnlResult;
}

vtkSmartPointer<vtkPolyData> SurfaceSubdivisionReconstruction::BreakDaughter(vtkSmartPointer<vtkPolyData> parentMesh, vtkSmartPointer<vtkPolyData> daughterMesh, int indexOfPath){
    // Use to check whether it is inside or not


    vtkSmartPointer<vtkCellArray> polys = daughterMesh->GetPolys();

    vtkSmartPointer<vtkCellArray> newPolys =  vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPoints> points = daughterMesh->GetPoints();
    int numberOfCells = polys->GetNumberOfCells();

    polys->InitTraversal();

    double tempPoint[3];
    for(int j = 0; j < numberOfCells ; j++){
        vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();

        polys->GetNextCell(list);

        vtkSmartPointer<vtkPoints> testPoints = vtkSmartPointer<vtkPoints>::New();

        for(int i = 0; i < list->GetNumberOfIds(); i++){
           points->GetPoint(list->GetId(i), tempPoint);
           testPoints->InsertNextPoint(tempPoint);
        }
        // test points have the triangles ...

        bool toRemove = false;
        for(int k = 0; k < indexOfPath; k++){
            if ( ArePointsInPath(k, testPoints)){
                //cellsToRemove.push_back(j);
                toRemove = true;
                break;
            }
        }

        if (!toRemove){
            vtkSmartPointer<vtkTriangle> triangle1 =  vtkSmartPointer<vtkTriangle>::New();

            triangle1->GetPointIds()->SetId(0,list->GetId(0));
            triangle1->GetPointIds()->SetId(1,list->GetId(1));
            triangle1->GetPointIds()->SetId(2,list->GetId(2));
            newPolys->InsertNextCell(triangle1);
        }

    }


    vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();
    newPolyData->SetPoints(points);
    newPolyData->SetPolys(newPolys);
    newPolyData->Modified();

    vtkSmartPointer<vtkCleanPolyData> cleanPolyData =  vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputData(newPolyData);
    cleanPolyData->Update();


    auto result = cleanPolyData->GetOutput();


    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivity = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivity->AddInputData(result);
    connectivity->SetExtractionModeToLargestRegion();
    connectivity->Update();

    auto largest = connectivity->GetOutput();

    std::cout << result->GetNumberOfPoints() << " / " << points->GetNumberOfPoints() << std::endl;

   return largest;
}

bool SurfaceSubdivisionReconstruction::ArePointsInPath(int pathToCheck, vtkSmartPointer<vtkPoints> pointsToCheck){

    tubNav::ContainerIterator<double> it(vascTree);
    it.InitTraversal();
    while ( it.GetCurrentPathId() != pathToCheck ){
        it.MoveNextPath();
    }
    // We have now found the path that we wanted to generate
    auto frame = it.GetCoordinateFrame();
    frame->InitTraversal();

    bool inside[pointsToCheck->GetNumberOfPoints()];
    for(int i = 0; i < pointsToCheck->GetNumberOfPoints(); i++)
        inside[i] = false;



    while (!it.IsAtEnd() && it.GetCurrentPathId() == pathToCheck ){
        if ( frame->IsAtEnd()){
            frame = it.GetCoordinateFrame();
            frame->InitTraversal();
        }
        //it.IsPointInsideOtherSegment()
        // Now frame has the points and the generating part... let's see if this works
        auto point = frame->GetCenterPoint();

        // now we check whether the points are within the radius ...

        for(int k = 0; k < pointsToCheck->GetNumberOfPoints(); k++){

            double pt[3];
            pointsToCheck->GetPoint(k, pt);

            double distance =  sqrt( pow(point.getX() - pt[0],2.0) + pow(point.getY() - pt[1],2.0) + pow(point.getZ() - pt[2],2.0));

            if ( distance <= point.getR()*0.95){
                inside[k] = true;
            }
        }

        frame->Next();
        ++it;
    }

    int totalInside = 0;
    for(int i = 0; i < pointsToCheck->GetNumberOfPoints(); i++)
          if (inside[i]) totalInside++;

    return totalInside == pointsToCheck->GetNumberOfPoints();
}



vtkSmartPointer<vtkPolyData>   SurfaceSubdivisionReconstruction::GetAllSimpleQuadReconstruction(){

    vtkSmartPointer<vtkPolyData> overallMesh = GetSinglePathReconstruction(0,40);

    std::cout  << "Get Number of Paths " << vascTree->GetNumberOfPaths() << std::endl;
    for(int i = 1; i <= 3/*vascTree->GetNumberOfPaths()*/; i++){

        std::cout << "***********************************" << std::endl;

        vtkSmartPointer<vtkPolyData> DaughterMesh = GetSinglePathReconstruction(i,40);

        if  ( DaughterMesh->GetNumberOfPolys() < 40 ) continue;

        std::cout << i << " :: " << DaughterMesh->GetNumberOfPolys() << std::endl;
        std::cout << overallMesh->GetNumberOfPolys() << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOp = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
        booleanOp->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
        booleanOp->SetInputData(0, overallMesh);
        booleanOp->SetInputData(1, DaughterMesh);
        booleanOp->Update();

        std::cout << "A" << std::endl;
        auto res = booleanOp->GetOutput();

        std::cout << "res ? " << res->GetNumberOfPolys() << std::endl;
        vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        normals->SetInputData(res);
        normals->ConsistencyOn();
        normals->SplittingOff();
        normals->Update();

        std::cout << "B" << std::endl;

        std::cout << normals->GetOutput()->GetNumberOfPolys() << std::endl;
        vtkSmartPointer<vtkCleanPolyData> cleanPolyData =  vtkSmartPointer<vtkCleanPolyData>::New();
        cleanPolyData->SetInputData(normals->GetOutput());
        cleanPolyData->Update();



        std::cout << "C" << std::endl;

        auto resCleaned = cleanPolyData->GetOutput();
        std::cout << resCleaned->GetNumberOfPolys() << std::endl;

        std::cout << "D " << std::endl;
        overallMesh->DeepCopy(resCleaned);
    //    auto result = cleanPolyData->GetOutput();
    }


    return overallMesh;
}

vtkSmartPointer<vtkPolyData> SurfaceSubdivisionReconstruction::GetReconstructionFromPointCloud(){

    vtkSmartPointer<vtkPoints> allPoints = vtkSmartPointer<vtkPoints>::New();
    double tmpPoint[3];
    for(int i = 0; i < vascTree->GetNumberOfPaths(); i++){
        auto pointsFromPath = GetSinglePathPoints(i, 40);

        for(int j = 0; j < pointsFromPath->GetNumberOfPoints(); j++){
            pointsFromPath->GetPoint(j, tmpPoint);
            allPoints->InsertNextPoint(tmpPoint);
        }
    }

    // Neiteher vtkSurfaceReconstructionFilter no Delaunay3D worked to generate the points from
    // this point cloud...
    vtkSmartPointer<vtkPolyData> tmpPoly = vtkSmartPointer<vtkPolyData>::New();
    tmpPoly->SetPoints(allPoints);

    return tmpPoly;
}


