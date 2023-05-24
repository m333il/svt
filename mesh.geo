// Gmsh project created on Wed May 24 21:50:53 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {3, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {1, 4};
//+
Line(4) = {1, 2};
//+
Line(5) = {3, 1};
//+
Curve Loop(1) = {4, -1, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -3, -5};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {5} = 10 Using Bump 1;
//+
Transfinite Curve {1} = 10 Using Progression 3;
//+
Transfinite Curve {4} = 10 Using Progression 3;
//+
Transfinite Curve {3} = 10 Using Progression 3;
//+
Transfinite Curve {2} = 10 Using Progression 3;
