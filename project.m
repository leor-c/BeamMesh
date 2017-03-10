%% For the project:
% first init the mesh, and all necessary condition:
mesh = Mesh('../hw2/bunny2.off');

%   number of sites for Voronoi diagram:
k = 30;



bunny_bm = BeamMesh('../hw2/bunny2.off',k);

bunny_cvt = CentroidalVoronoiTesselation(bunny_bm, bunny_bm.sites);
[voronoiDi, ~] = bunny_cvt.run();
%   show CVT results:
bunny_bm.showMesh(ones(bunny_bm.dimensions(2),1));
bunny_bm.showMesh(bunny_cvt.cells);
alpha(0.95);
hold on;
[V_v,V_e] = bunny_cvt.findVoronoiVertices();
scatter3(bunny_cvt.sites(:,1),bunny_cvt.sites(:,2),bunny_cvt.sites(:,3),100,'filled','yellow');
%scatter3(V_e(:,1),V_e(:,2),V_e(:,3),25,'filled','cyan');
scatter3(V_v(:,1),V_v(:,2),V_v(:,3),50,'filled','red');


%% For Bunny:
%   set the number of sites for the Voronoi diagram:
k = 40;
%   load the origial mesh from off file:
bunny = Mesh('../meshes/bunny2_smooth.off');
%bunny.showMesh();
bunny_bm = BeamMesh(bunny,k);

%   Calculate CVT w/ respect to mean curvature:
bunny_cvt = CentroidalVoronoiTesselation(bunny, bunny_bm.sites);
%   show CVT results:
bunny.showMesh(bunny_cvt.cells);
alpha(0.85);
hold on;
[V_v,V_e] = bunny_cvt.findVoronoiVertices();
scatter3(bunny_cvt.sites(:,1),bunny_cvt.sites(:,2),bunny_cvt.sites(:,3),100,'filled','yellow');
scatter3(V_e(:,1),V_e(:,2),V_e(:,3),25,'filled','cyan');
scatter3(V_v(:,1),V_v(:,2),V_v(:,3),50,'filled','red');
%% For Car:
%   set the number of sites for the Voronoi diagram:
k = 70;
%   load the origial mesh from off file:
car = Mesh('../meshes/Octane_Fixed.off');
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
k = 40;
%   load the origial mesh from off file:
cow = Mesh('../meshes/spot_s4.off');
%car.showMesh();
cow_bm = BeamMesh(cow,k);

%   Calculate CVT w/ respect to mean curvature:
cow_cvt = CentroidalVoronoiTesselation(cow, cow_bm.sites);
%   show CVT results:
cow.showMesh(cow_cvt.cells);
alpha(0.85);
hold on;
[V_v,V_e] = cow_cvt.findVoronoiVertices();
scatter3(cow_cvt.sites(:,1),cow_cvt.sites(:,2),cow_cvt.sites(:,3),100,'filled','yellow');
scatter3(V_e(:,1),V_e(:,2),V_e(:,3),25,'filled','cyan');
scatter3(V_v(:,1),V_v(:,2),V_v(:,3),50,'filled','red');

