#ifndef __TUBNAVPOINT_H
#define __TUBNAVPOINT_H

/*
   Base Centerline Point, can be better described as a
   node. It needs a the physical location & the radius.
   
   The tube could be either 2D or 3D, if 2D then we default 
   the Z construction to 0 
*/

#include <type_traits>
#include <iostream>
#include <limits>


namespace tubNav {

   // The only types permitted can have varying precision,
   // but all of them should support arithmetic opeations

   template<typename T,
            typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
   class Point {

       public:
          // Regular 3D construction
          Point() :_x(0), _y(0), _z(0), _r(0){

          }

          Point(const T x,const T y,const T r) : _x(x), _y(y), _z(0), _r(r)
          {

          }

          Point(const T pt[]): _x(pt[0]), _y(pt[1]),_z(pt[2]), _r(0){

          }
 
          //2D Construction
          Point(const T x,const T y,const T z,const T r) : _x(x), _y(y), _z(z),_r(r) {

          }

          // Copy constructor
          Point( const Point& rhs): _x( rhs.getX()), _y(rhs.getY()), _z(rhs.getZ()),_r(rhs.getR())   {

          }

          // Move constructor
          Point( Point&& rhs): _x( rhs.getX()), _y(rhs.getY()), _z(rhs.getZ()),_r(rhs.getR())  {

          }

          // Default assignation
          Point& operator=(const Point&& rhs) {
              _x = rhs.getX();
              _y = rhs.getY();
              _z = rhs.getZ();
              _r = rhs.getR();
              return *this;
          }

          Point& operator=(const Point& rhs) {
              _x = rhs.getX();
              _y = rhs.getY();
              _z = rhs.getZ();
              _r = rhs.getR();
              return *this;
          }


          bool operator==(const Point& o) const{
              int tot = 0;
              for(int i = 0; i < 4; i++){

                  if ( o.get(i) > get(i)  && o.get(i) - get(i) < std::numeric_limits<T>::epsilon() ){ // o is greater then the difference should be positive
                      tot++;
                  }
                  else if ( get(i) >= o.get(i) &&  get(i) - o.get(i) < std::numeric_limits<T>::epsilon() ){  // get(i) is greater or equal
                      tot++;
                  }
              }
              return tot == 4;
          }


          bool operator!=(const Point& o) const{
              int tot = 0;
              for(int i = 0; i < 4; i++){

                  if ( o.get(i) > get(i)  && o.get(i) - get(i) < std::numeric_limits<T>::epsilon() ){ // o is greater then the difference should be positive
                      tot++;
                  }
                  else if ( get(i) >= o.get(i) &&  get(i) - o.get(i) < std::numeric_limits<T>::epsilon() ){  // get(i) is greater or equal
                      tot++;
                  }
              }
              return tot != 4;
          }

          // Operator + and - may be useful for vector calculations
          // The radius is not useful here so we can initalize to 0

          // Operator+, given that Return Value optimization exist
          // Then we just generate the sum of the locations. 
          virtual Point operator+(Point&& o) const noexcept {
              Point newPoint(_x + o.getX(),
                                       _y + o.getY(),
                                       _z + o.getZ(),
                                       0);

              return newPoint;
          }
          
          virtual Point operator+(Point& o) const noexcept {
              Point newPoint(_x + o.getX(),
                                       _y + o.getY(),
                                       _z + o.getZ(),
                                       0);

              return newPoint;
          }

          virtual Point operator-(Point<T>&& o) const noexcept {
              Point newPoint(_x - o.getX(),
                             _y - o.getY(),
                             _z - o.getZ(),
                                  0);
              return newPoint;
          }

          virtual Point operator-(Point<T>& o) const noexcept {
              Point newPoint(_x - o.getX(),
                             _y - o.getY(),
                             _z - o.getZ(),
                                       0);
              return newPoint;
          }
          
          virtual Point operator-(const Point<T>& o) const noexcept {
              Point newPoint(_x - o.getX(),
                             _y - o.getY(),
                             _z - o.getZ(),
                                       0);
              return newPoint;
          }
          virtual Point operator*(const T& value) const {
              Point newPoint(value * _x,
                             value * _y,
                             value * _z,
                             0);
              return newPoint;
          }

          virtual Point operator /(const T& value) const {
              Point newPoint(_x / value,
                             _y / value,
                             _z / value,
                             0);
              return newPoint;
          }



          virtual T operator[]( const int idx) const {
              return get(idx);
          }
          virtual T dot(Point&& o) const noexcept {
              T dot = _x * o.getX() +  _y * o.getY() +  _z * o.getZ();
              return dot ;
          }

          virtual T dot(Point& o) const noexcept {
              T dot = _x * o.getX() +  _y * o.getY() +  _z * o.getZ();
              return dot;
          }
          
         virtual Point cross(Point&& o) const noexcept {
              Point newPoint(_y * o.getZ() - _z* o.getY(),
                             _z * o.getX() - _x* o.getZ(),
                             _x * o.getY() - _y* o.getX(),
                             0);
              return newPoint;
          }

          virtual Point cross(Point& o) const noexcept {
               Point newPoint(_y * o.getZ() - _z* o.getY(),
                              _z * o.getX() - _x* o.getZ(),
                              _x * o.getY() - _y* o.getX(),
                              0);
               return newPoint;
           }

          virtual void normalize(T (*sqrt)(T)) {
              T size = sqrt(normSquared());
              _x /= size;
              _y /= size;
              _z /= size;
           }
         
          // The Sqrt function has several operators which may cause a change of precision due
          // to implicit conversions, therefore only supply the normSquared
          virtual T normSquared() const noexcept {
             return _x*_x + _y*_y + _z*_z;
          }
         
          // Other points may be generated as long as they support
          // arithmetic operations
          virtual ~Point() = default;

          friend std::ostream& operator<<( std::ostream& o, const Point& t ) {
             std::cout <<  "Coors: [" << t.getX() << "," << t.getY() << "," << t.getZ() << "]" << " R " << t.getR(); 
             return o;
          }

          bool IsZeroth(){ // Check if it unitialized, by looking at the comparison with the zeroth
                 Point pt;
                 return equals(pt);
          }

          bool equals(const double o[]){
              int tot = 0;

              for(int i = 0; i < 3; i++){

                  if ( o[i] > get(i)  && o[i] - get(i) < std::numeric_limits<T>::epsilon() ){ // o is greater then the difference should be positive
                      tot++;
                  }
                  else if ( get(i) >= o[i] && get(i) - o[i] < std::numeric_limits<T>::epsilon() ){  // get(i) is greater or equal
                      tot++;
                  }
              }
              return (tot == 3);
          }

          bool equals(const Point o){
              int tot = 0;

              for(int i = 0; i < 4; i++){

                  if ( o.get(i) > get(i)  && o.get(i) - get(i) < std::numeric_limits<T>::epsilon() ){ // o is greater then the difference should be positive
                      tot++;
                  }
                  else if ( get(i) >= o.get(i) && get(i) - o.get(i) < std::numeric_limits<T>::epsilon() ){  // get(i) is greater or equal
                      tot++;
                  }
              }
              return (tot == 4);
          }


          T distance(const Point o, T (*sqrt)(T)){
              Point newPoint(_x - o.getX(),
                             _y - o.getY(),
                             _z - o.getZ(),
                             0);


              return sqrt(newPoint.normSquared());
          }

          T frustumVolume(const Point o, T (*sqrt)(T)){
              T h = distance(o, sqrt)*3.14159;
              T volume = (o.getR()*o.getR() + _r*o.getR() + _r*_r )*h;
              volume /= 3.0;
              return volume;
          }


          const T getX() const { return _x;}
          const T getY() const { return _y;}
          const T getZ() const { return _z;}
          const T getR() const { return _r;}

          const T get(int dim) const{
              T value = std::numeric_limits<T>::max();

              switch (dim) {
              case 0: value = _x; break;
              case 1: value = _y; break;
              case 2: value = _z; break;
              case 3: value = _r; break;

              default:
                  break;
              }
              return value;
          }
       private: 
         T _x;
         T _y;
         T _z;
         T _r;
   };

}

#endif
