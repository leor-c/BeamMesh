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
%bunny_cvt.showResults();

%% For Dragon:
%   set the number of sites for the Voronoi diagram:
k = 80;
%   load the origial mesh from off file:
dragon = Mesh('../meshes/dragon_50kV.off');
%car.showMesh();
dragon_bm = BeamMesh(dragon,k);

%   Calculate CVT w/ respect to mean curvature:
dragon_cvt = CentroidalVoronoiTesselation(dragon, dragon_bm.sites);
%   show CVT results:
dragon_cvt.showResults();

%% For Companion Cube:
%   set the number of sites for the Voronoi diagram:
k = 120;
%   load the origial mesh from off file:
cube = Mesh('../meshes/companionCubeRemashed.off');
%car.showMesh();
cube_bm = BeamMesh(cube,k);

%   Calculate CVT w/ respect to mean curvature:
cube_cvt = CentroidalVoronoiTesselation(cube, cube_bm.sites);
%   show CVT results:
%cube_cvt.showResults();

%% For Cow:
%   set the number of sites for the Voronoi diagram:
k = 100;
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
k = 60;
%   load the origial mesh from off file:
moo = Mesh('../meshes/bunny2.off');
%car.showMesh();
moo_bm = BeamMesh(moo,k);

%   Calculate CVT w/ respect to mean curvature:
moo_cvt = CentroidalVoronoiTesselation(moo, moo_bm.sites);
%   show CVT results:
moo_cvt.showResults();

%% For Chamelion:
%   set the number of sites for the Voronoi diagram:
k = 50;
%   load the origial mesh from off file:
chamelion = Mesh('../meshes/chamelion_light.off');
%car.showMesh();
chamelion_bm = BeamMesh(chamelion,k);

%   Calculate CVT w/ respect to mean curvature:
chamelion_cvt = CentroidalVoronoiTesselation(chamelion, chamelion_bm.sites);
%   show CVT results:
chamelion_cvt.showResults();

%%  Spaceship2

k=50;

spaceship = Mesh('../meshes/spaceship2.off');

spaceship_bm = BeamMesh(spaceship,k);

spaceship_cvt = CentroidalVoronoiTesselation(spaceship, spaceship_bm.sites);

spaceship_cvt.showResults();




