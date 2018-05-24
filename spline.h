/*
  Library for 1-d splines
  Copyright Ryan Michael
  Licensed under the LGPLv3 
*/
#ifndef spline_h
#define spline_h

#include <math.h>

#define Hermite 10
#define Catmull 11

template<class T>
class Spline
{
  public:
    Spline( void );
    Spline( T x[], T y[], int numPoints, int degree = 1 );
    Spline( T x[], T y[], T m[], int numPoints );
    T value( T x );
    void setPoints( T x[], T y[], int numPoints );
    void setPoints( T x[], T y[], T m[], int numPoints );
    void setDegree( int degree );
    
  private:
    T calc( T, int);
    T* _x;
    T* _y;
    T* _m;
    int _degree;
    int _length;
    int _prev_point;
    
    T hermite( T t, T p0, T p1, T m0, T m1, T x0, T x1 );
    inline T hermite_00( T t );
    inline T hermite_10( T t );
    inline T hermite_01( T t );
    inline T hermite_11( T t );
    T catmull_tangent( int i );
};

template <typename T>
Spline<T>::Spline(void) {
  _prev_point = 0;
}

template <typename T>
Spline<T>::Spline( T x[], T y[], int numPoints, int degree )
{
  setPoints(x, y, numPoints);
  setDegree(degree);
  _prev_point = 0;
}

template <typename T>
Spline<T>::Spline( T x[], T y[], T m[], int numPoints )
{
  setPoints(x,y,m,numPoints);
  setDegree(Hermite);
  _prev_point = 0;
}

template <typename T>
void Spline<T>::setPoints( T x[], T y[], int numPoints ) {
  _x = x;
  _y = y;
  _length = numPoints;
}

template <typename T>
void Spline<T>::setPoints( T x[], T y[], T m[], int numPoints ) {
  _x = x;
  _y = y;
  _m = m;
  _length = numPoints;
}

template <typename T>
void Spline<T>::setDegree( int degree ){
  _degree = degree;
}

template <typename T>
T Spline<T>::value( T x )
{
  if( _x[0] > x ) { 
    return _y[0]; 
  }
  else if ( _x[_length-1] < x ) { 
    return _y[_length-1]; 
  }
  else {
    for(int i = 0; i < _length; i++ )
    {
      int index = ( i + _prev_point ) % _length;
      
      if( _x[index] == x ) {
        _prev_point = index;
        return _y[index];
      } else if( (_x[index] < x) && (x < _x[index+1]) ) {
        _prev_point = index;
        return calc( x, index );
      }
    }    
  }
}

template <typename T>
T Spline<T>::calc( T x, int i )
{
  switch( _degree ) {
    case 0:
      return _y[i];
    case 1:
      if( _x[i] == _x[i+1] ) {
        // Avoids division by 0
        return _y[i];
      } else {
        return _y[i] + (_y[i+1] - _y[i]) * ( x - _x[i]) / ( _x[i+1] - _x[i] );      
      }
    case Hermite:
      return hermite( ((x-_x[i]) / (_x[i+1]-_x[i])), _y[i], _y[i+1], _m[i], _m[i+1], _x[i], _x[i+1] );
    case Catmull:
      if( i == 0 ) {
        // x prior to spline start - first point used to determine tangent
        return _y[1];
      } else if( i == _length-2 ) {
        // x after spline end - last point used to determine tangent
        return _y[_length-2];
      } else {
        T t = (x-_x[i]) / (_x[i+1]-_x[i]);
        T m0 = (i==0 ? 0 : catmull_tangent(i) );
        T m1 = (i==_length-1 ? 0 : catmull_tangent(i+1) );
        return hermite( t, _y[i], _y[i+1], m0, m1, _x[i], _x[i+1]);        
      }
  }
}

template <typename T>
T Spline<T>::hermite( T t, T p0, T p1, T m0, T m1, T x0, T x1 ) {
  return (hermite_00(t)*p0) + (hermite_10(t)*(x1-x0)*m0) + (hermite_01(t)*p1) + (hermite_11(t)*(x1-x0)*m1);
}
template <>
inline double Spline<double>::hermite_00( double t ) { return (2*pow(t,3)) - (3*pow(t,2)) + 1;}
template <>
inline double Spline<double>::hermite_10( double t ) { return pow(t,3) - (2*pow(t,2)) + t; }
template <>
inline double Spline<double>::hermite_01( double t ) { return (3*pow(t,2)) - (2*pow(t,3)); }
template <>
inline double Spline<double>::hermite_11( double t ) { return pow(t,3) - pow(t,2); }

template <>
inline float Spline<float>::hermite_00( float t ) { return (2*powf(t,3)) - (3*powf(t,2)) + 1;}
template <>
inline float Spline<float>::hermite_10( float t ) { return powf(t,3) - (2*powf(t,2)) + t; }
template <>
inline float Spline<float>::hermite_01( float t ) { return (3*powf(t,2)) - (2*powf(t,3)); }
template <>
inline float Spline<float>::hermite_11( float t ) { return powf(t,3) - powf(t,2); }

template <typename T>
T Spline<T>::catmull_tangent( int i ) 
{ 
  if( _x[i+1] == _x[i-1] ) {
    // Avoids division by 0
    return 0;
  } else {  
    return (_y[i+1] - _y[i-1]) / (_x[i+1] - _x[i-1]);    
  } 
}

#endif