%%  Report Code:
%   Bunny Mesh - Initialization and first look:
%   load the origial mesh from off file:
bunny = Mesh('../meshes/bunny2.off');
%   show the bunny with uniform color:
bunny.showMesh(ones(bunny.dimensions(2),1));
%%  Sampling Initial sites:
%   set the number of sites for the Voronoi diagram:
k = 60;
bunny_bm = BeamMesh(bunny,k);
sampled_sites = bunny_bm.sites;

%%  Generating the CVT:
%   Calculate CVT w/ respect to mean curvature:
bunny_cvt = CentroidalVoronoiTesselation(bunny, bunny_bm.sites);

%   show initial sites:
bunny_cvt.showInitialSites();

%   show CVT results:
bunny_cvt.showResults();

%%  Show CVT polygon resulting from the CVT process:

bunny_cvt.showPolygon();

%%  Generating the BeamMesh:

if(bunny_cvt.validateCVT())
    bunny_bm.importDataFromCVT(bunny_cvt);
    %   show before convergence:
    bunny_bm.showBeamMesh();
    
    %bunny_bm.convergeConstraints();
    bunny_bm.convergeConstraints(false);
    bunny_bm.showBeamMesh();
    
    bunny_bm.computeRadii();
    bunny_bm.computeThickBeamMesh();
end

