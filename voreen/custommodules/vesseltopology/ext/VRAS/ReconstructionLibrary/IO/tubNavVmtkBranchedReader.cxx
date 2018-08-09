
#include "tubNavVmtkBranchedReader.h"

#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>


namespace tubNav {

   bool BranchedCenterlinesReader::IsBranchingData(vtkSmartPointer<vtkPolyData> polydata){
       vtkSmartPointer<vtkIntArray> groupData = vtkIntArray::SafeDownCast(polydata->GetCellData()->GetArray("GroupIds"));
       bool hasBranchingArray = false;
       if (groupData)
           hasBranchingArray = true;
       return hasBranchingArray;

   }


   void BranchedCenterlinesReader::RemoveDuplicates(int firstId[], vtkSmartPointer<vtkPoints> points, const double epsilon){
       const int numPoints = points->GetNumberOfPoints();
       // Given a point i, find which one is the first id
       for(int i = 0; i < numPoints; i++ ){
            Point<double> current = points->GetPoint(i);
            firstId[i] = numPoints + 2; // Just a failsafe, it should cause a crash if
            // it doesn't get a location...
            for(int j = 0; j < numPoints; j++){
                Point<double> other = points->GetPoint(j);

                if ( current.distance(other, sqrt) < epsilon )
                {
                    firstId[i] = j;
                    break;
                }
            }

            // It may be possible that the point are mapped to the same id but are just a bit farther
            // than each other, so  it will create a chain to itself
            while ( firstId[i] != firstId[firstId[i]])
               firstId[i] = firstId[firstId[i]];
        }

   }


   int BranchedCenterlinesReader::GetNextInnerPoints(vtkSmartPointer<vtkIdList> list,  std::unordered_map<int, std::set<int>* >& nextPoint, int firstId[]){
       int lastId = -1;
       for(int j = 0; j < list->GetNumberOfIds(); j++ ){

                       if ( j < list->GetNumberOfIds() -1){
                           vtkIdType ptId = list->GetId(j);
                           vtkIdType ptIdNxt = list->GetId(j+1);

                           int actualId = firstId[ptId];
                           int actualNextId = firstId[ptIdNxt];

                           if ( actualId != actualNextId){
                               // If the actual point is not found then we add it,
                               // otherwise we just get it
                               if ( nextPoint.find(actualId) == nextPoint.end()){
                                   //It would be better to make it as a single step perhaps
                                   std::set<int>* tmp = new std::set<int>();
                                   std::pair<int,std::set<int>* >  keyValuePair(actualId, tmp);
                                   nextPoint.insert(keyValuePair);
                               }
                               nextPoint[actualId]->insert(actualNextId);
                           }
                       }
                       else if ( j +1 == list->GetNumberOfIds() ){
                           lastId = list->GetId(j);
                       }
       }
       return lastId;
   }


   void BranchedCenterlinesReader::SetConnectionLastPoint(const int lastId, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkPolyData> polyData,
                                                          std::unordered_map<int, std::set<int>* >& nextPoint, int firstId[] ){
       double minDistance = 1; // It should be below 1
       // The last point in the current line. Find the closest point
       int actualId = firstId[lastId];
       double pt[3];
       points->GetPoint(actualId, pt);
       Point<double> endPoint = pt;
       int continuation = -1;

       vtkSmartPointer<vtkCellArray> others = vtkSmartPointer<vtkCellArray>::New();
       others->DeepCopy(polyData->GetLines());

       //vtkSmartPointer<vtkCellArray> others = polyData->GetLines();
       others->InitTraversal();
       for(int j = 0; j < polyData->GetNumberOfLines(); j++){
           vtkSmartPointer<vtkIdList> list2 = vtkSmartPointer<vtkIdList>::New();

           others->GetNextCell(list2);
           int o = 0;
           vtkIdType ptId = list2->GetId(o); // the first
           while (firstId[ptId] == actualId ){
               o++;
               ptId = list2->GetId(o);
           }

           Point<double> startPoint = points->GetPoint(firstId[ptId]);
           if ( startPoint.distance(endPoint, sqrt) < minDistance )
           {
               minDistance = startPoint.distance(endPoint, sqrt);
               continuation = firstId[ptId];
           }
       }

       if ( continuation != -1){
           std::set<int>* tmp = new std::set<int>();
           std::pair<int,std::set<int>* >  keyValuePair(actualId, tmp);
           nextPoint.insert(keyValuePair);
           nextPoint[actualId]->insert(continuation);
       }
   }

   void BranchedCenterlinesReader::CreateGraph(vtkSmartPointer<vtkPolyData> polyData,
                                               std::vector< std::shared_ptr< Segment<double> > > &segments,
                                               std::vector< std::pair<int, int> >& connections){
        // We have the polydata.
        vtkSmartPointer<vtkPoints> points = polyData->GetPoints();

         // Need to find the roots, bifurcations
         // First the points are actually repeated several times. It's easier if we compare the ids of the points
         const int numPoints = points->GetNumberOfPoints();
         int firstId[numPoints];
         RemoveDuplicates(firstId, points, 0.2);

        // Now we have all the points, without repetition and in a single coordinate locations..
        // The lines in the branching vmtk are separated in the branching section and the segments
        // but they are not necessarily divided in the branching point...
        std::unordered_map<int, std::set<int>* > nextPoint;
        vtkSmartPointer<vtkCellArray> lines = polyData->GetLines();

        lines->InitTraversal();
        for(int i = 0; i < polyData->GetNumberOfLines(); i++){
            vtkSmartPointer<vtkIdList> list= vtkSmartPointer<vtkIdList>::New();
            lines->GetNextCell(list);
            int numPointsInLine = list->GetNumberOfIds();
            if (numPointsInLine <= 0) continue;
            int lastId = GetNextInnerPoints(list, nextPoint, firstId);
            // Now for the case of the last point
            SetConnectionLastPoint(lastId, points, polyData, nextPoint, firstId);
        }
        // Remember to clear the memory of the points....



        //std::vector< Segment<double> > segments;
        //std::vector< std::pair<int, int> > connections;

        std::vector<int> visited;
        connections.push_back( std::make_pair(-1, 0 ));

        CreateSegments(segments, 0, nextPoint, firstId, polyData, connections, visited);


        std::cout << "Total segments " << segments.size() << std::endl;
        CleanMapMemory(nextPoint);
   }


   void BranchedCenterlinesReader::CleanMapMemory( std::unordered_map<int, std::set<int>* >& nextPoint){

       for(auto it : nextPoint){
           delete it.second;
       }
   }


   int BranchedCenterlinesReader::CreateSegments(std::vector< std::shared_ptr< Segment<double> > > &segments, int currentId,
                                                   std::unordered_map<int, std::set<int>* >& nextPoint, int firstId[],
                                                   vtkSmartPointer<vtkPolyData> polydata,  std::vector< std::pair<int, int> >& connections , std::vector<int> &visited)
   {
           // Create a new segment

          if (  std::find(visited.begin(), visited.end(), currentId) != visited.end())
              return -1;


          visited.push_back(currentId);

          std::shared_ptr< Segment<double> > segment = std::make_shared<Segment<double>>();
          vtkSmartPointer<vtkDoubleArray> radiusData = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("MaximumInscribedSphereRadius"));
          double pt[3], r;
          polydata->GetPoints()->GetPoint(currentId, pt);
          r = radiusData->GetTuple(currentId)[0];

          Point<double> cpt(pt[0], pt[1], pt[2], r);

          segment->AddPoint(cpt);

          while( nextPoint.find(currentId) != nextPoint.end() && nextPoint[firstId[currentId]]->size() == 1  ){

                  std::set<int>::iterator first = nextPoint[firstId[currentId]]->begin();
                  int nxtId = firstId[*first];

                  polydata->GetPoints()->GetPoint(nxtId, pt);
                  r = radiusData->GetTuple(nxtId)[0];
                  Point<double> npt(pt[0], pt[1], pt[2], r);
                  segment->AddPoint(npt);
                  currentId = nxtId;
          }


          segments.push_back(segment);
          int loc = segments.size() -1;

          if (  nextPoint.find(currentId) != nextPoint.end() && nextPoint[firstId[currentId]]->size() >1 ){

              std::set<int>::iterator furcations = nextPoint[firstId[currentId]]->begin();

              for( ; furcations != nextPoint[firstId[currentId]]->end(); furcations++ )
              {
                  int furcation = *furcations;
                  int o = CreateSegments(segments, furcation, nextPoint, firstId, polydata, connections, visited);

                  if ( o != -1){
                      connections.push_back(  std::make_pair(loc, o));
                  }
              }
          }
          return loc;
   }



}
