Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 0, 1, 1};
Point(4) = {0, 0, 1, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};

Transfinite Surface {6};
Transfinite Line {1,3} = 24 Using Progression 1;
Transfinite Line {2,4} = 2 Using Progression 1;

Recombine Surface{6};
Physical Line("fixAll") = {2,4};
Physical Surface("Membrane") = {6};
