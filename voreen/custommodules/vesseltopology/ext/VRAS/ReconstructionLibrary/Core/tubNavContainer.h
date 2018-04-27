#ifndef __TUBNAVCONTAINER_H
#define __TUBNAVCONTAINER_H

/*
   Vascular Reconstruction Library.
   Vascular Container
   This Container holds a vascular graph.
   Here we define the two main representations of the bifurcations for the centerlines.
   Need to find a better definition
   1. As in topological thinning
   2. Fast marching




   A graph implementation.
   Each Node is either a root, bifurcation or leaf location.
   The combination between two nodes is a segment
   Each node can have multiple parents or multiple children.


   // A DFS implementation needs to be used to calculate the
*/

#include <type_traits>
#include <vector>
#include <memory>
#include <atomic>
#include <unordered_map>
#include <queue>
#include <stack>
#include <map>
#include <set>

#include "tubNavSegment.h"
#include "../MathCore/CoordinateFrame.h"


namespace tubNav {


   // A node is either the root, bifurcation or leaf.
   // defined by a begin value and a identifier,
   // the begin
   template<typename T,
   typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
   struct Node {

       Node(){
           // Imposible location
           begin = -99999999;

       }
       // A list of parents and a list of children
       int begin; // In the connections, which is the starting point value
       // -1 for root

       Point<T> stlocation;

       std::set<std::shared_ptr<Node> > parents;
       std::set<std::shared_ptr<Node> > children;
       std::set<int> childrenConnections;

    };


   // A key is a connection, from-to. Id1 is the starting point
   // Id2 is the ending point. The root is -1, so the combination of
   // a starting and end point is a segment.
   struct Key {
       int id1;
       int id2;
       double geodesicDistance;
       double volume;
       double ratioVolumeDistance;

       bool operator==(const Key &other) const
         { return (id1 == other.id1
                   && id2 == other.id2);
         }
   };

   struct KeyHasher
   {
     std::size_t operator()(const Key& k) const
     {
       using std::size_t;
       using std::hash;

       return ((hash<int>()(k.id1)
                ^ (hash<int>()(k.id2) << 1)) >> 1);
     }
   };


   static bool ContainerComparisonHelper(const std::pair<int,float> &a,const std::pair<int,float> &b)
   {
          return a.second>b.second;
   }


   template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
   class ContainerIterator ;


    template<typename T,
    typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
    class Container {

        friend class ContainerIterator<T>;

    public:
        Container(): root(nullptr), size(0){

        }  

        void CreateContainer(std::vector< std::shared_ptr< Segment<double> > > &segments,
                             std::vector<std::pair<int, int> >& connections){

            // The segments are the collection of points.. in between "node" locations
            // The connections are basically where this segment start and end...
            // since we start from -1 as the root node, the segment is
            int numNodes = 0;
            for(auto connection: connections){
                Key newKey;
                newKey.id1 = connection.first;
                newKey.id2 = connection.second;
                std::shared_ptr<Segment<double>> segmentToAdd = segments.at(connection.second);
                std::shared_ptr<Segment<T> > newSegment = segmentToAdd; //std::make_shared<Segment<T>>(segmentToAdd);

                //************
                // Calculate some helpful information of the segment, the geodesic distance
                // The volume and the ratio
                SegmentIterator<double> sit(newSegment);
                sit.InitTraversal();
                T volume = 0;
                T geodesicDistance = 0;
                bool enters = true;
                while((!sit.IsNextEnd()) || enters){ // This is done for segments that are only 2 points long.. it enters at least once

                    Point<T> p1 = sit.GetPoint();
                    Point<T> p2 = sit.GetNextPoint();

                    // We use the frustum of a cone formula to get the volume
                    geodesicDistance += p1.distance(p2, sqrt);
                    volume += p1.frustumVolume(p2, sqrt);
                    ++sit;
                    enters = false;
                }
                newKey.geodesicDistance = geodesicDistance;
                newKey.volume = volume;
                newKey.ratioVolumeDistance = volume/geodesicDistance;

                // add information of a key, save information of the
                // segment
                _keyMap[std::make_pair(connection.first, connection.second)] = newKey;
                //*************
                auto keyValuePair = std::make_pair(newKey, newSegment);
                _segments.insert(keyValuePair);
                // If the begin is not found
                if ( _nodes.find(connection.first) == _nodes.end()){
                    std::shared_ptr<Node<T>> newNode( new Node<T>());
                    newNode->stlocation = segmentToAdd->GetStart();
                    // each node is defined by its starting connection
                    newNode->begin = connection.first;
                    // and the children are where it connects to
                    newNode->childrenConnections.insert(connection.second);
                    auto keyValuePair = std::make_pair(newNode->begin, std::move(newNode));
                    _nodes.insert(keyValuePair);
                    numNodes++;
                }
                else {
                    // same starting point connection,
                    _nodes[connection.first]->childrenConnections.insert(connection.second);
                }

                if (_nodes.find(connection.second) == _nodes.end()){
                    //In this case the node is a leaf,
                    std::shared_ptr<Node<T>> newNode( new Node<T>());
                    std::shared_ptr<Segment<double>> segmentToAdd = segments.at(connection.second);
                    newNode->stlocation = segmentToAdd->GetEnd();
                    newNode->begin = connection.second;
                    auto keyValuePair = std::make_pair(newNode->begin, std::move(newNode));
                    _nodes.insert(keyValuePair);
                    _nodes[connection.second]->parents.insert(_nodes[connection.first]);
                    numNodes++;
                }

            }
            // Now we have the nodes, with the children connections...
            for(auto connection: connections){
                // Create the parents or children
                // Adding parents, ending connections
                if ( _nodes.find(connection.second) != _nodes.end()){
                    _nodes[connection.second]->parents.insert(_nodes[connection.first]);
                }

                // In the other direction as well
                if ( _nodes.find(connection.first) != _nodes.end()){
                    // Now we add the children
                    auto currentNode = _nodes[connection.first];
                    for(auto child: currentNode->childrenConnections){
                        if ( _nodes.find(child) != _nodes.end()){
                            currentNode->children.insert(_nodes[child]);
                        }
                    }
                }
            }


            root = _nodes[-1];

            std::map< std::pair<int,int>, bool > seenSegments;
            std::vector<int> path(1,-1);

            _dfsPaths.clear();
            size = 0;
            DFS(-1, path, seenSegments);
            SortPaths();

            std::cout << "After Sorting " << std::endl;
            int p = 0;
            for(auto path: _dfsPaths){
                _startOfPathIdx.push_back(size);
                size += path.size();
                //
                double vol = 0;
                std::cout << p << " >> " ;

                for(unsigned int i = 0; i < path.size(); i++){
                    _segmentToPathIdx.push_back(p);
                    auto key = _keyMap[ path.at(i)];
                    vol +=  key.ratioVolumeDistance;
                    std::cout << "[" << path.at(i).first << "," << path.at(i).second << "]   ";

                }
                std::cout << " :: vol " << vol << std::endl;

                //
                p++;
            }
            GenerateFrames(connections);
            OverlapsOfSegments();
        }

        void CheckSegment(std::shared_ptr<Segment<double> > segment1, std::shared_ptr<Segment<double> > segment2 ){
            segment1->SetToStart();

            for(unsigned int k = 0; k < segment1->GetSize(); k++){
                auto point = segment1->GetCurrentPoint();

                segment2->SetToStart();
                for(unsigned int h = 0; h < segment2->GetSize(); h++ ){
                     auto otherPoint = segment2->GetCurrentPoint();
                     if ( point.distance(otherPoint, sqrt)  <= otherPoint.getR() ){
                         segment1->SetCurrentPointInsideProperty(true);
                         break;
                     }
                     segment2->MoveForward();
                }
                segment1->MoveForward();
            }
        }

        void OverlapsOfSegments(){
            // The paths are already sorted...
            // so for each path, let's look at the previous paths and see whether it is inside of the previous ones
            for(unsigned int i = 0; i < _dfsPaths.size(); i++){
                //
                auto path = _dfsPaths.at(i);

                for(unsigned int j =  0; j < i; j++){
                    // look at the previous paths
                    auto prevPath = _dfsPaths.at(j);
                    //std::vector< std::vector< std::pair<int,int> >> _dfsPaths;
                    for(unsigned int k = 0; k < path.size(); k++){
                        auto keyOfSegment = _keyMap[ path.at(k)];
                        auto segmentToCheck = _segments[keyOfSegment];

                        for(unsigned int h = 0; h < prevPath.size(); h++){
                            auto keyOfSegmentAgainst = _keyMap[ prevPath.at(k)];
                            auto segmentToCheckAgainst = _segments[keyOfSegmentAgainst];

                            CheckSegment(segmentToCheck, segmentToCheckAgainst);
                        }
                    }
                }
            }
        }

        void SortPaths(){
            int p = 0;
            std::vector< std::pair<int,float> > pathsVolumes;
            for(auto path: _dfsPaths){
                double vol = 0;
                std::cout << p << " >> " ;
                for(unsigned int i = 0; i < path.size(); i++){

                    auto key = _keyMap[ path.at(i)];

                    std::cout << "[" << path.at(i).first << "," << path.at(i).second << "(  " << key.geodesicDistance << "," << key.volume  << ") ] ";
                    vol +=  key.geodesicDistance; //key.ratioVolumeDistance;
                }

                std::cout << " :: vol " << vol << std::endl;
                pathsVolumes.push_back(std::make_pair(p, vol));
                p++;
            }
            std::sort(pathsVolumes.begin(), pathsVolumes.end(), ContainerComparisonHelper);
            //    std::vector< std::vector< std::pair<int,int> >> _dfsPaths;
            std::vector< std::vector< std::pair<int,int> >> tmpPaths;

            for(unsigned int i = 0; i < pathsVolumes.size(); i++){
                tmpPaths.push_back( _dfsPaths.at(pathsVolumes.at(i).first));
            }
            _dfsPaths.clear();
            for(auto path: tmpPaths)
                _dfsPaths.push_back(path);

        }

        int GetNumberOfPaths(){
            return _dfsPaths.size();
        }

        void GetPath(int idx, bool* success,  std::shared_ptr< Segment<T> > fnlSegment){
            if ( idx < 0 || idx >= static_cast<int>(_dfsPaths.size()) ){
                *success = false;
                return;
            }
            auto selectedPath = _dfsPaths.at(idx);
            // each element is a segment ...

            for(auto segmentKey: selectedPath){
                auto currentSegment = _segments[_keyMap[segmentKey]];
                std::vector<Point<double>> points;
                currentSegment->GetAllPoints(points);
                for(auto point: points)
                  fnlSegment->AddPoint(point);
            }

            *success = true;
            return;
        }


        //TODO - Handle loops
        void DFS(int loc, std::vector<int>&  currentPath, std::map< std::pair<int,int>, bool >& seenSegments){
            // From loc

            if ( _nodes.find(loc) != _nodes.end() ){
                auto currentNode = _nodes[loc];
                std::vector< std::pair<int,float> > connections;

                for(auto connection: currentNode->childrenConnections){                
                    auto key = _keyMap[std::make_pair(loc, connection) ];
                    connections.push_back( std::make_pair(connection, key.ratioVolumeDistance) );
                }

                std::sort(connections.begin(), connections.end(), ContainerComparisonHelper);

                for(auto connection: connections){
                    auto key = _keyMap[std::make_pair(loc, connection.first) ];
                    currentPath.push_back(key.id2);
                    DFS( key.id2, currentPath, seenSegments);
                    currentPath.pop_back();
                }

                // If the node is a leaf
                if ( connections.empty()){
                    std::vector< std::pair<int, int> > finalPath;
                    for(unsigned int i = 0; i < currentPath.size() -1; i++){
                        auto pair = std::make_pair(currentPath.at(i), currentPath.at(i+1));
                        if ( seenSegments.find(pair)  == seenSegments.end() ){
                            seenSegments[pair] = true;
                            finalPath.emplace_back(pair);
                        }
                    }

                    if (finalPath.size() > 0)
                        _dfsPaths.push_back(finalPath);

                }
            }

        }

        std::shared_ptr< Segment<T> > GetSegment(int idx){
            auto key = GetSegmentKey(idx);
            return _segments[key];
        }

        std::shared_ptr< CoordinateFrame<T> > GetCoordinateFrame(int idx){
            auto key = GetSegmentKey(idx);
            return _frames[key];
        }

        Key GetSegmentKey(int idx){
            int p = 0;
            //std::cout << "Getting segment key? " << idx << " .. ";
            std::pair<int, int> selectedPath;
            for(auto path: _dfsPaths){
                for(unsigned int i = 0; i < path.size(); i++){
                    if (p == idx){
                         selectedPath = path.at(i);
                         //std::cout << "selected Path? " << path.at(i).first << " , " << path.at(i).second << " :: ";
                    }
                    p++;
                }
            }
            auto key = _keyMap[selectedPath];       
            //std::cout << " .. [" << key.id1 << "," << key.id2 << "]" << std::endl;
            return key;
        }

        int GetSegmentIndex(Key key){
            int p = 0;
            int index = 0;
            for(auto path: _dfsPaths){
                for(auto segment: path){

                    if (segment.first == key.id1 && segment.second == key.id2 ){
                        index = p;
                    }
                    p++;
                }
            }
            return index;
        }

        std::size_t GetNumberOfDaughters(int idx){
            auto key = GetSegmentKey(idx);
            auto node = _nodes[key.id2];
            std::size_t numDaughters = node->childrenConnections.size();
            return numDaughters;
        }




     private:
         std::shared_ptr<Node<T>> root;
         std::unordered_map<int, std::shared_ptr<Node<T>> > _nodes;
         std::unordered_map<Key, std::shared_ptr< Segment<T> >, KeyHasher > _segments;
         std::unordered_map<Key, std::shared_ptr< CoordinateFrame<T> > , KeyHasher> _frames;

         std::map< std::pair<int,int>, Key> _keyMap;
         std::vector< std::vector< std::pair<int,int> >> _dfsPaths;
         std::vector<int> _segmentToPathIdx;
         std::vector<int> _startOfPathIdx;
         std::vector<int> _endOfPathIdx;

         std::size_t size;


         void GenerateFrames( std::vector<std::pair<int, int> >& connections){
              auto key = GetSegmentKey(0); // Root key
              std::shared_ptr< CoordinateFrame<T> > firstFrame = std::make_shared<CoordinateFrame<T>>();
              firstFrame->setNavigationCurve(_segments[key], true);

              auto keyValuePair = std::make_pair(key, firstFrame);
              _frames.insert(keyValuePair);

              // For every key, find the previous vector & generate the
              // coordinate frame according to the root. That way the whole path
              // is always from the root with minimal rotation & coming back as well
              std::queue<Key> connecting;
              connecting.push(key);

              while(!connecting.empty()){
                  auto currentKey = connecting.front();
                  connecting.pop();
                  Point<T> lastVect = _frames[currentKey]->GetLastReferenceVector();

                  for(auto connection: connections){
                      if ( connection.first == currentKey.id2){
                          auto nextKey = _keyMap[connection];

                          std::shared_ptr< CoordinateFrame<T> > nextFrame = std::make_shared<CoordinateFrame<T>>();
                          nextFrame->SetReferenceVector(lastVect);
                          nextFrame->setNavigationCurve(_segments[nextKey]);
                          auto keyValuePair = std::make_pair(nextKey, nextFrame);
                          _frames.insert(keyValuePair);
                          connecting.push(nextKey);
                      }
                  }
              }
         }


    };




    template<typename T, typename U >
    class ContainerIterator : public std::iterator<std::random_access_iterator_tag, Container<T> >{

         public:
            ContainerIterator():_container(nullptr), _segmentIterator(nullptr), _pathIndex(0), _currentSegment(0),_daughterIndex(0)  {
                _nextSegmentId = 0;
                _prevSegmentId = 0;
                _overallPosition = 0;

            }

            void SetContainer(std::shared_ptr<Container<T> > container){
                 _container = container;
                 _currentSegment = 0; _daughterIndex = 0; _pathIndex = 0;
                 _segmentIterator = SegmentIterator<T>(nullptr);
                 _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
                 _nextSegmentId = _currentSegment+1;
                 _prevSegmentId = 0;
                 _overallPosition = 0;
            }

           ContainerIterator(std::shared_ptr<Container<T> > container):_container(container), _segmentIterator(nullptr), _pathIndex(0), _currentSegment(0),_daughterIndex(0)  {
               _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
               _nextSegmentId = _currentSegment+1;
               _prevSegmentId = 0;
               _overallPosition = 0;
           }

           ContainerIterator& operator++(){
               // If we get the end, we pass to the next segment
               if (_segmentIterator.IsAtEnd()){
                   MoveNextSegment();
               }
               else {
                    ++_segmentIterator;
               }

               _overallPosition++;
               _daughterIndex = 0;
               return *this;
           }


           ContainerIterator& operator--(){

               if (_segmentIterator.IsAtStart()){
                   MovePrevSegment();
               }
               else {
                   --_segmentIterator;
               }

               _overallPosition--;
               _daughterIndex = 0;
               return *this;}


           Point<T> GetPoint(){ return _segmentIterator.GetPoint(); }
           std::shared_ptr< Segment<T>> GetSegment(){  return _container->GetSegment(_currentSegment); }
           std::shared_ptr< CoordinateFrame<T>> GetCoordinateFrame(){  return _container->GetCoordinateFrame(_currentSegment); }

           int GetCurrentPathId(){ return _pathIndex;}
           int GetNextPathId(){
               if ( _currentSegment +1 >= _container->_segmentToPathIdx.size())
                   return -1;

               return _container->_segmentToPathIdx.at(_currentSegment+1);
           }

           void MoveNextSegment(){

             if ( _nextSegmentId < _container->size){
                _currentSegment = _nextSegmentId ;
                DefinePrevSegment();
                _nextSegmentId = _currentSegment +1;
                _daughterIndex = 0;
                _pathIndex = _container->_segmentToPathIdx.at(_currentSegment);
                // std::cout << "Segment? " << _currentSegment <<    ", " << _container->GetSegmentKey(_currentSegment).id1 << " , " << _container->GetSegmentKey(_currentSegment).id2 << std::endl;
                _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
             }
           }

           bool IsAtEndSegment(){ return (_nextSegmentId) == _container->size;}
           bool IsAtEnd(){ return IsAtEndSegment() && _segmentIterator.IsAtEnd() ;}
           bool IsBeforeBifurcation(){ return _segmentIterator.IsNextEnd(); }

           bool IsBifurcation(){ return _segmentIterator.IsAtEnd() ; }

           std::size_t GetNumberOfDaughters(){
               return _container->GetNumberOfDaughters(_currentSegment);
           }

           void ChangeToNextDaugther(){
              //TODO-here the overall position should be changed
              //     to an appropiate value
              auto key = _container->GetSegmentKey(_currentSegment);
              auto node = _container->_nodes[key.id2];
              std::vector< std::pair<int, float> > childrenRatios;
              int idx = 0;
              for(auto child: node->children){
                    auto keyPair = std::make_pair(node->begin, child->begin);
                    auto childKey = _container->_keyMap[keyPair];
                    childrenRatios.push_back( std::make_pair( idx,childKey.ratioVolumeDistance ));
                    idx++;
              }
              std::sort(childrenRatios.begin(), childrenRatios.end(), ContainerComparisonHelper);
              _daughterIndex = (_daughterIndex+1)% node->children.size();
              idx = 0;
              for(auto child: node->children){
                  if ( childrenRatios.at(idx).first == _daughterIndex ){
                      auto keyPair = std::make_pair(node->begin, child->begin);
                      auto childKey = _container->_keyMap[keyPair];
                      _nextSegmentId = _container->GetSegmentIndex(childKey);
                  }
                  idx++;
              }
           }

           bool IsPointInsideOtherSegment(){ return _segmentIterator.IsInside(); }

           bool IsRoot(){ return _pathIndex == 0 && _currentSegment == 0;}

           bool IsLeaf(){
               if(IsAtEndSegment() && _segmentIterator.IsNextEnd()) return true;

               if (IsBifurcation()){
                   std::size_t nextPathIndex = _container->_segmentToPathIdx.at(_currentSegment+1);
                   if (_pathIndex != nextPathIndex)
                       return true;
               }
               return false;
           }

           bool IsAtStart(){
               return _overallPosition == 0;

           }
           void InitTraversal(){ _pathIndex = 0; _currentSegment = 0;  _daughterIndex = 0;
                                 _prevSegmentId = 0; _nextSegmentId = _currentSegment +1;
                                 _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
                                 _overallPosition = 0;
                               }



           void DefinePrevSegment(){
               auto key = _container->GetSegmentKey(_currentSegment);
               auto node = _container->_nodes[key.id1];

               std::vector< std::pair<int, float> > parentRatios;
               int idx = 0;
               for(auto parent: node->parents ){
                   auto keyPair = std::make_pair(parent->begin, node->begin);
                   auto parentKey = _container->_keyMap[keyPair];
                   parentRatios.push_back( std::make_pair( idx, parentKey.ratioVolumeDistance ));
                   idx++;
               }
               std::sort(parentRatios.begin(), parentRatios.end(), ContainerComparisonHelper);

               idx = 0;
               for(auto parent: node->parents ){ //
                   if ( parentRatios.at(0).first == idx ){
                       auto keyPair = std::make_pair(parent->begin, node->begin);
                       auto parentKey = _container->_keyMap[keyPair];
                       _prevSegmentId = _container->GetSegmentIndex(parentKey);
                       break;
                   }
                   idx++;
               }
           }

           void MovePrevSegment(){
               if (_currentSegment > 0){
                  _nextSegmentId = _currentSegment;
                  _currentSegment = _prevSegmentId;
                  DefinePrevSegment();
                  _pathIndex = _container->_segmentToPathIdx.at(_currentSegment);

                  _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
                  _daughterIndex = 0;
               }
           }

           void MoveNextPath(){
               //TODO-here the overall position should be changed
               //     to an appropiate value
               if (_pathIndex +1 < _container->_dfsPaths.size() ){
                   _pathIndex += 1;
                   _currentSegment = _container->_startOfPathIdx.at(_pathIndex);
                   _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
                   _daughterIndex = 0;
               }
           }
           void MovePrevPath(){
               //TODO-here the overall position should be changed
               //     to an appropiate value
               if (_pathIndex > 0){
                   _pathIndex -= 1;
                   _currentSegment = _container->_startOfPathIdx.at(_pathIndex);
                   _segmentIterator.SetCurve(_container->GetSegment(_currentSegment));
                   _daughterIndex = 0;
               }
           }

           void GetSize(){
                return _container->size;
           }
         private:
           std::shared_ptr<Container<T>> _container;
           SegmentIterator<T> _segmentIterator;

           std::size_t _pathIndex;
           std::size_t _currentSegment;
           std::size_t _nextSegmentId;
           int _prevSegmentId;
           int _daughterIndex;
           int _overallPosition;

    };
}



#endif
