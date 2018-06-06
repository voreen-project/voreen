#ifndef __TUBNAV_VMTKBRANCHEDCENTERLINES_READER
#define __TUBNAV_VMTKBRANCHEDCENTERLINES_READER

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <memory>
#include "../Core/tubNavContainer.h"
#include <set>
#include <unordered_map>

namespace tubNav {

      class BranchedCenterlinesReader {

           public:
               bool IsBranchingData(vtkSmartPointer<vtkPolyData> polydata);

               void CreateGraph(vtkSmartPointer<vtkPolyData> polyData, std::vector<std::shared_ptr<Segment<double> > > &segments, std::vector<std::pair<int, int> > &connections);

               void RemoveDuplicates(int firstId[], vtkSmartPointer<vtkPoints> points, const double epsilon = 0.2);

               int GetNextInnerPoints(vtkSmartPointer<vtkIdList> list,  std::unordered_map<int, std::set<int>* >& nextPoint, int firstId[]);
               void SetConnectionLastPoint(const int lastId, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkPolyData> polyData,
                                            std::unordered_map<int, std::set<int>* >& nextPoint, int firstId[]);


               void CleanMapMemory( std::unordered_map<int, std::set<int>* >& nextPoint);


               int CreateSegments(std::vector< std::shared_ptr< Segment<double> > > &segments, int currentId, std::unordered_map<int, std::set<int>* >& nextPoint, int firstId[],
                                   vtkSmartPointer<vtkPolyData> polydata, std::vector< std::pair<int, int> >& connections, std::vector<int>& visited);



               template<typename T,
                        typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
               std::shared_ptr< Container<T> > CreateContainer(vtkSmartPointer<vtkPolyData> polydata ){

                   // The polydata is saved already with a double precision, we don't need that to
                   // specify the type to read the data, but once we read it, we need to convert it to the user
                   // specified datatype
                   std::vector< std::shared_ptr< Segment<double> > > segments;
                   std::vector<std::pair<int, int> > connections;
                   CreateGraph(polydata, segments, connections);
                   //OverlapsOfSegments(segments);

                   std::shared_ptr< Container<T> > container(new Container<T>());
                   container->CreateContainer(segments, connections);

                   return container;
               }

      };
}

#endif
