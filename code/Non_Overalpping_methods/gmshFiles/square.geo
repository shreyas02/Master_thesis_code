// Gmsh project created on Wed Oct 23 16:03:23 2024

//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {05, 0, 0, 1.0};
Point(3) = {05, 05, 0, 1.0};
Point(4) = {0, 05, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {4, 3, 2, 1} = 21 Using Progression 1;
//+
Transfinite Surface {1} = {1, 4, 3, 2};
//+
Recombine Surface {1};
//+
Physical Curve("sides", 5) = {4, 3, 2, 1};
//+
Physical Point("sides", 6) = {4, 3, 2, 1};
//+
Physical Surface("domain", 7) = {1};
