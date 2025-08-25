
//  a, b, c      points on the unit 2-sphere
//  for best numerics, the shortest side should be opposite c
//
//  returns *signed* area of the triangle.
//  positive area is counterclockwise, when viewed from the outside of the sphere

extern  double  area_spherical_triangle( const double a[3], const double b[3], const double c[3] );


//  compute signed area of a planar polygon
//  counter-clockwise is positive

extern  double  area_polygon( const double x[], const double y[], int n );


//  compute 3x3 determinant
extern  double  det3x3( const double a[3], const double b[3], const double c[3] );
