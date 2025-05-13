// Gmsh project created on Wed Oct 23 16:03:23 2024

//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {05, 0, 0, 1.0};
Point(3) = {05, 05, 0, 1.0};
Point(4) = {0, 05, 0, 1.0};
Point(5) = {3, 0, 0, 1.0};
Point(6) = {3, 5, 0, 1.0};
Line(1) = {4, 6};
Line(2) = {6, 3};
Line(3) = {3, 2};
Line(4) = {2, 5};
Line(5) = {5, 1};
Line(6) = {1, 4};
Line(7) = {6, 5};

//+
Curve Loop(1) = {-7, -1, -6, -5};
Plane Surface(1) = {1};
Curve Loop(2) = {7, -4, -3, -2};
Plane Surface(2) = {2};

//+
Transfinite Curve {6,7,3} = 21 Using Progression 1;
Transfinite Curve {5, 1} = 13 Using Progression 1; 
Transfinite Curve {2, 4} = 9 Using Progression 1;
Transfinite Surface {1} = {1, 5, 6, 4};
Transfinite Surface {2} = {5, 2, 3, 6};
Recombine Surface {1, 2};

//+
Physical Surface("LeftDomain", 8) = {1};
Physical Surface("RightDomain", 9) = {2};
Physical Curve("interface", 10) = {7};
Physical Curve("RightBoundary", 11) = {4, 3, 2};
Physical Curve("LeftBoundary", 12) = {1, 6, 5};
Physical Point("TopPoint", 13) = {6};
Physical Point("BottomPoint", 14) = {5};
Physical Point("points", 15) = {1, 2, 3, 4};