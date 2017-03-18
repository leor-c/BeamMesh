%% For Smooth Bunny:
%   set the number of sites for the Voronoi diagram:
k = 100;
%   load the origial mesh from off file:
bunnySmooth = Mesh('../meshes/bunny2_smooth.off');
%bunny.showMesh();
bunnySmooth_bm = BeamMesh(bunnySmooth,k);

%   Calculate CVT w/ respect to mean curvature:
bunnySmooth_cvt = CentroidalVoronoiTesselation(bunnySmooth, bunnySmooth_bm.sites);
%   show CVT results:
bunnySmooth_cvt.showResults();

if(bunnySmooth_cvt.validateCVT())
    bunnySmooth_bm.importDataFromCVT(bunnySmooth_cvt);
    bunnySmooth_bm.convergeConstraints(false);
    %bunny_bm.convergeConstraints();
    %bunnySmooth_bm.showBeamMesh();
    
    bunnySmooth_bm.computeRadii();
    bunnySmooth_bm.computeThickBeamMesh();
end

%% For Regular Bunny - Report Code:
%   load the origial mesh from off file:
bunny = Mesh('../meshes/bunny2.off');
%bunny.showMesh();
%   set the number of sites for the Voronoi diagram:
k = 60;
bunny_bm = BeamMesh(bunny,k);

%   Calculate CVT w/ respect to mean curvature:
bunny_cvt = CentroidalVoronoiTesselation(bunny, bunny_bm.sites);
%   show CVT results:
bunny_cvt.showResults();

if(bunny_cvt.validateCVT())
    bunny_bm.importDataFromCVT(bunny_cvt);
    bunny_bm.convergeConstraints();
    %bunny_bm.convergeConstraints(false);
    bunny_bm.showBeamMesh();
    
    bunny_bm.computeRadii();
    bunny_bm.computeThickBeamMesh();
end

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
if(dragon_cvt.validateCVT())
    dragon_bm.importDataFromCVT(dragon_cvt);
    dragon_bm.convergeConstraints(false);
    dragon_bm.showBeamMesh();
    dragon_bm.computeRadii();
    dragon_bm.computeThickBeamMesh();
end

%% For Companion Cube:
%   set the number of sites for the Voronoi diagram:
k = 100;
%   load the origial mesh from off file:
cube = Mesh('../meshes/companionCubeRemashed.off');
%car.showMesh();
cube_bm = BeamMesh(cube,k);

%   Calculate CVT w/ respect to mean curvature:
cube_cvt = CentroidalVoronoiTesselation(cube, cube_bm.sites);
%   show CVT results:
cube_cvt.showResults();

if(cube_cvt.validateCVT())
    cube_bm.importDataFromCVT(cube_cvt);
    cube_bm.convergeConstraints(false);
    cube_bm.showBeamMesh();
    cube_bm.computeRadii();
    cube_bm.computeThickBeamMesh();
end

%% For Cow:
%   set the number of sites for the Voronoi diagram:
k = 70;
%   load the origial mesh from off file:
cow = Mesh('../meshes/spot_s4.off');
%car.showMesh();
cow_bm = BeamMesh(cow,k);

%   Calculate CVT w/ respect to mean curvature:
cow_cvt = CentroidalVoronoiTesselation(cow, cow_bm.sites);
%   show CVT results:
cow_cvt.showResults();

if(cow_cvt.validateCVT())
    cow_bm.importDataFromCVT(cow_cvt);
    cow_bm.convergeConstraints(false);
    cow_bm.showBeamMesh();
    
    cow_bm.computeRadii();
    cow_bm.computeThickBeamMesh();
end

%% For moomoo:
%   set the number of sites for the Voronoi diagram:
k = 100;
%   load the origial mesh from off file:
moo = Mesh('../meshes/moomoo_s3.off');
%car.showMesh();
moo_bm = BeamMesh(moo,k);

%   Calculate CVT w/ respect to mean curvature:
moo_cvt = CentroidalVoronoiTesselation(moo, moo_bm.sites);
%   show CVT results:
moo_cvt.showResults();

if(moo_cvt.validateCVT())
    moo_bm.importDataFromCVT(moo_cvt);
    moo_bm.convergeConstraints(false);
    moo_bm.showBeamMesh();
    moo_bm.computeRadii();
    moo_bm.computeThickBeamMesh();
end

%% For Chamelion:
%   set the number of sites for the Voronoi diagram:
k = 80;
%   load the origial mesh from off file:
chamelion = Mesh('../meshes/chamelion_light.off');
%car.showMesh();
chamelion_bm = BeamMesh(chamelion,k);

%   Calculate CVT w/ respect to mean curvature:
chamelion_cvt = CentroidalVoronoiTesselation(chamelion, chamelion_bm.sites);
%   show CVT results:
chamelion_cvt.showResults();


%%  Cat mesh:

k = 120;

cat = Mesh('../meshes/cat_s3.off');
cat_bm = BeamMesh(cat,k);
cat_cvt = CentroidalVoronoiTesselation(cat, cat_bm.sites);
cat_cvt.showResults();
cat_bm.importDataFromCVT(cat_cvt);
cat_bm.convergeConstraints(false);
cat_bm.computeRadii();
cat_bm.computeThickBeamMesh();



