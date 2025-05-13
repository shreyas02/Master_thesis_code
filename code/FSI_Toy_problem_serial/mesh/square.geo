// Gmsh project created on Wed Oct 23 16:03:23 2024

//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {06, 0, 0, 1.0};
Point(3) = {06, 01, 0, 1.0};
Point(4) = {0, 01, 0, 1.0};
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
Transfinite Curve {3, 1} = 41 Using Progression 1;
//+
Transfinite Curve {4, 2} = 6 Using Progression 1;
//+
Transfinite Surface {1} = {1, 4, 3, 2};
//+
Recombine Surface {1};
//+
Physical Point("Left", 10) = {1 , 4};
//+
Physical Point("Right", 11) = {2, 3};
//+
Physical Curve("Right", 12) = {2};
//+
Physical Curve("Left", 13) = {4};
//+
Physical Curve("Bottom", 14) = {1};
//+
Physical Curve("Top", 15) = {3};
//+
Physical Point("BCTop", 17) = {4, 3};
//+
Physical Surface("Domain", 9) = {1};
//+
Physical Point("Bottom", 18) = {1 , 2};
//+