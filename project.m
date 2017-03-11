%% For Bunny:
%   set the number of sites for the Voronoi diagram:
k = 40;
%   load the origial mesh from off file:
bunny = Mesh('../meshes/bunny2_smooth.off');
bunny.showMesh();
bunny_bm = BeamMesh(bunny,k);

%   Calculate CVT w/ respect to mean curvature:
bunny_cvt = CentroidalVoronoiTesselation(bunny, bunny_bm.sites);
%   show CVT results:
bunny_cvt.showResults();

%% For Car:
%   set the number of sites for the Voronoi diagram:
k = 60;
%   load the origial mesh from off file:
car = Mesh('../meshes/companionCubeRemashed.off');
%car.showMesh();
car_bm = BeamMesh(car,k);

%   Calculate CVT w/ respect to mean curvature:
car_cvt = CentroidalVoronoiTesselation(car, car_bm.sites);
%   show CVT results:
car.showMesh(car_cvt.cells);
alpha(0.85);
hold on;
[V_v,V_e] = car_cvt.findVoronoiVertices();
scatter3(car_cvt.sites(:,1),car_cvt.sites(:,2),car_cvt.sites(:,3),100,'filled','yellow');
scatter3(V_e(:,1),V_e(:,2),V_e(:,3),25,'filled','cyan');
scatter3(V_v(:,1),V_v(:,2),V_v(:,3),50,'filled','red');

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
k = 100;
%   load the origial mesh from off file:
moo = Mesh('../meshes/bunny2.off');
%car.showMesh();
moo_bm = BeamMesh(moo,k);

%   Calculate CVT w/ respect to mean curvature:
moo_cvt = CentroidalVoronoiTesselation(moo, moo_bm.sites);
%   show CVT results:
moo_cvt.showResults();




