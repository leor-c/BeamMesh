classdef BeamMesh < handle
    %   Class for Beam Meshes. 
    %   * sites is the points x_i that determines the Voronoi cell.
    %   * VD_vertices - vertices of the voronoi diagram (not neccessarily
    %   needed)
    %   * vertices_plus/minus is the v_i+, v_i- coordinates of the vertices of all
    %   quads.
    properties
        vertices_plus
        vertices_minus
        faces
        verticesNormals
        numberOfSites
        metricTensors
        sites
        VD_vertices
        VD_adjacencyMatrix
    end
    
    methods
        function obj = BeamMesh(shape, k)
            %   Initialize a Beam Mesh. 
            %   Params: 
            %       -   shape: an input mesh of the shape to approximate
            %       -   k: number of sites \ cells
            
            
            
            obj.numberOfSites = k;
            
            %   init metric tensors:
            metricTensors = zeros(3,3,k);
            for i=1:k
                metricTensors(:,:,i) = eye(3); 
            end
            
            obj.sampleInitialSites(shape);
            
            %   Compute the Centroidal Voronoi Tesselation for the mesh.
            %obj.CVT(shape);
            
            
            
        end
        
        function sampleInitialSites(obj, shape)
            H = shape.getMeanCurvature();
            % TRY SMOOTHING!:
            %L = shape.getLaplacian();
            %H = smooth(H,20);
            %H = smooth(H,20);
            [~, H_F] = shape.vertex_to_face_interp(H);
            %shape.showMesh(H_F);

            %   Use mean curvature as density function, to randomly sample points:
            %   compute the cumulative distribution function
            F_X = H_F;
            sum = 0;
            for i=1:length(H_F)
                sum = sum + H_F(i);
                F_X(i) = sum;
            end
            F_X = F_X ./ sum;

            %   Use random variable transformation: use F_X to find the inverse by
            %   finding for a value the last index that F_X <= x, so the uniform
            %   variable U (random_sample) has mean curvature distribution.
            uniform_random_sample = rand(obj.numberOfSites, 1);
            hat = zeros(shape.dimensions(2), 1);
            obj.sites = zeros(obj.numberOfSites, 3);
            for i=1:obj.numberOfSites
                face = find(F_X <= uniform_random_sample(i), 1, 'last');
                hat(face) = 1;
                %   choose barycenter of chosen faces as x_i:
                for j=1:3
                    obj.sites(i,:) = obj.sites(i,:) + shape.vertices(shape.faces(face,j),:)./3;
                end
            end
        end
        
        function CVT(obj, shape)
            %   Compute the Centroidal Voronoi Tesselation for the mesh.
            %   (Voronoi diagram).
            cvt = CentroidalVoronoiTesselation(shape, obj.sites);
            obj.importDataFromCVT(cvt);
        end
        
        function importDataFromCVT(obj, cvtObj)
            %   This enables you to generate beam mesh from an existing
            %   CVT instance (class CentroidalVoronoiTesselation).
            obj.sites = cvtObj.sites;
            
            obj.VD_adjacencyMatrix = cvtObj.voronoiAdjMatrix;
            obj.VD_vertices = 1000.*cvtObj.voronoiVertices;
            
            obj.vertices_plus = obj.VD_vertices + cvtObj.voronoiVertexNormals;
            obj.vertices_minus = obj.VD_vertices - cvtObj.voronoiVertexNormals;
        end
        
        function generateBeamVertices(obj)
            %   As described, for each vertex of the Voronoi Diagram, 
            %   calculate the normal of the vertex, then set v_i+-
            %   accordingly. Normals are calculated using sum of normals to
            %   adjacent faces (cross products of adjacent incident edge
            %   vectors).
            TODO implement this function
        end
        
        function projectBeamConstraints(obj)
            
        end
        
        function planarity_projection = getPlanarityProjection(obj)
            %   Calculate the projection for the planarity constraint.
            %   For each row in adj. matrix, there will be exactly 3
            %   nonzero values. because adj. matrix is symmetric, let's
            %   operate only on upper triange. for each edge (that define a
            %   beam) for all beam edges, project it.
            [n,~] = size(obj.VD_vertices);
            %   for each edge in VD, there are 4 edges in the beam mesh.
            %   number of edges is (3/2 * 2)|VD_vertices| = 3*|VD_vertices|
            planarity_constraints = 12*n;
            planarity_projection = zeros(planarity_constraints,1);
            
            for i=1:n
                % for vertex i in VD_vertices:
                %   get neighbors to v_i from upper triangle of adj. mat.
                neighbors = find(obj.VD_adjacencyMatrix(i,i:end));
                for j = neighbors
                    TODO project edges of beam defined by v_i and v_j
                end
            end
        end
        
        function showBeamMesh(obj)
            %   Visualize the beam mesh:
            [V_size,~] = size(obj.VD_adjacencyMatrix);
            
            figure;
            
            for i=1:V_size
                vertices = find(obj.VD_adjacencyMatrix(i,:));
                for j=vertices
                    X = [obj.vertices_minus(i,1) obj.vertices_plus(i,1) ...
                        obj.vertices_plus(j,1) obj.vertices_minus(j,1)];
                    Y = [obj.vertices_minus(i,2) obj.vertices_plus(i,2) ...
                        obj.vertices_plus(j,2) obj.vertices_minus(j,2)];
                    Z = [obj.vertices_minus(i,3) obj.vertices_plus(i,3) ...
                        obj.vertices_plus(j,3) obj.vertices_minus(j,3)];
                    fill3(X,Y,Z, 'g');
                    hold on;
                end
            end
            
        end
        
        
    end
end