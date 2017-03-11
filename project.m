%% For Bunny:
%   set the number of sites for the Voronoi diagram:
k = 100;
%   load the origial mesh from off file:
bunny = Mesh('../meshes/bunny2_smooth.off');
%bunny.showMesh();
bunny_bm = BeamMesh(bunny,k);

%   Calculate CVT w/ respect to mean curvature:
bunny_cvt = CentroidalVoronoiTesselation(bunny, bunny_bm.sites);
%   show CVT results:
bunny_cvt.showResults();

%% For Car:
%   set the number of sites for the Voronoi diagram:
k = 50;
%   load the origial mesh from off file:
car = Mesh('../meshes/octane_smooth.off');
%car.showMesh();
car_bm = BeamMesh(car,k);

%   Calculate CVT w/ respect to mean curvature:
car_cvt = CentroidalVoronoiTesselation(car, car_bm.sites);
%   show CVT results:
car_cvt.showResults();

%% For Companion Cube:
%   set the number of sites for the Voronoi diagram:
k = 75;
%   load the origial mesh from off file:
cube = Mesh('../meshes/companionCubeRemashed.off');
%car.showMesh();
cube_bm = BeamMesh(cube,k);

%   Calculate CVT w/ respect to mean curvature:
cube_cvt = CentroidalVoronoiTesselation(cube, cube_bm.sites);
%   show CVT results:
cube_cvt.showResults();

%% For Cow:
%   set the number of sites for the Voronoi diagram:
k = 150;
%   load the origial mesh from off file:
cow = Mesh('../meshes/spot_s4.off');
%car.showMesh();
cow_bm = BeamMesh(cow,k);

%   Calculate CVT w/ respect to mean curvature:
cow_cvt = CentroidalVoronoiTesselation(cow, cow_bm.sites);
%   show CVT results:
cow_cvt.showResults();

%% For moomoo:
%   set the number of sites for the Voronoi diagram:
k = 50;
%   load the origial mesh from off file:
moo = Mesh('../meshes/bunny2.off');
%car.showMesh();
moo_bm = BeamMesh(moo,k);

%   Calculate CVT w/ respect to mean curvature:
moo_cvt = CentroidalVoronoiTesselation(moo, moo_bm.sites);
%   show CVT results:
moo_cvt.showResults();




