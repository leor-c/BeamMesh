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
        beamHeight
        heightTolerance
        minEdgeLengthThreshold
    end
    
    methods
        function obj = BeamMesh(shape, k, beamHeight, heightTolerance, edgeLenThreshold)
            %   Initialize a Beam Mesh. 
            %   Params: 
            %       -   shape: an input mesh of the shape to approximate
            %       -   k: number of sites \ cells
            %       -   beamHeight: the desired height for each beam. will
            %       be set to default value if not passed.
            %       -   heightTolerance: tolerance for each beam how far
            %       it can be from the target height. will be set to
            %       default value if not passed
            
            
            
            obj.numberOfSites = k;
            
            if (nargin <= 2)
                obj.beamHeight = 2;
                obj.heightTolerance = 0.02;
            elseif (nargin == 3)
                obj.beamHeight = beamHeight;
                obj.heightTolerance = 0.02;
            else
                obj.beamHeight = beamHeight;
                obj.heightTolerance = heightTolerance;
            end
            if (nargin < 5)
                obj.minEdgeLengthThreshold = 2;
            else
                obj.minEdgeLengthThreshold = edgeLenThreshold;
            end
            
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
            if (~cvtObj.validateCVT())
                disp('CVT is BAD, please build another one.');
                return;
            end
            obj.sites = cvtObj.sites;
            
            obj.VD_adjacencyMatrix = cvtObj.voronoiAdjMatrix;
            obj.VD_vertices = cvtObj.voronoiVertices;
            obj.verticesNormals = cvtObj.voronoiVertexNormals;
            
            %normailze the mesh coordinates so they are around 100:
            targetMean = 100;
            meandist = norm(mean(obj.VD_vertices));
            c = targetMean/meandist;
            obj.VD_vertices = c .* obj.VD_vertices;
            
            obj.vertices_plus = obj.VD_vertices + cvtObj.voronoiVertexNormals;
            obj.vertices_minus = obj.VD_vertices - cvtObj.voronoiVertexNormals;
        end
        
        function BeamEdgesMat = getBeamEdges(obj,i,j)
            %   calculate the beam generated from vertices i,j (+-).
            %   direction is from minus to plus.
            %   output matrix is composed of 4 column vectors, each
            %   represents an edge of the beam.
            %   column 1 is e_nrml_i, column 2 is e_nrml_j,
            %   column 3 is e_ofst_pl, column 4 is e_ofst_mi
            
            e_nrml_i = (obj.vertices_plus(i,:) - obj.vertices_minus(i,:))';
            e_nrml_j = (obj.vertices_plus(j,:) - obj.vertices_minus(j,:))';
            e_ofst_pl = (obj.vertices_plus(j,:) - obj.vertices_plus(i,:))';
            e_ofst_mi = (obj.vertices_minus(j,:) - obj.vertices_minus(i,:))';
            
            BeamEdgesMat = [e_nrml_i e_nrml_j e_ofst_pl e_ofst_mi];
            
        end
        
        function [v_i, v_j, e_ij] = getMidVertices(obj, i, j)
            %   get the vertices in the middle of the normal edges
            v_i = (obj.vertices_plus(i,:) + obj.vertices_minus(i,:))./2;
            v_j = (obj.vertices_plus(j,:) + obj.vertices_minus(j,:))./2;
            e_ij = (v_j - v_i)';
        end
        
        function [n_i, n_j] = getBeamNormals(obj, i, j)
            %   get normals to mid vertices v_i, v_j as described in
            %   article:
            [v_i, v_j, ~] = obj.getMidVertices(i, j);
            n_i = (obj.vertices_plus(i,:) - v_i)';
            n_j = (obj.vertices_plus(j,:) - v_j)';
        end
        
        function distance = projectBeamConstraints(obj)
            %   Projects all the constraints needed for the mesh as
            %   described in the article. I will implement it as I
            %   understand though, and build a matrix A so the global step
            %   is ||AV-p||2^2. The matrix A is constructed so each row
            %   will represent an edge (because we're projecting edges).
            %   For all projections Im assumming that A is constructed the
            %   following way: each row is |V+| + |V-| long, and first |V|
            %   indices are for v+ vertices and the last |V| are for v-.
            %   each row represents an edge of the projection.
            
            %   For this to work, it is neccessary that all vertices
            %   invovled will be centered at the mean of their positions.
            %   in this application, all vertices involved in all
            %   projections so I'll just center them once:
            %   Actually, it is neccessary for LS to be 0 if projection =
            %   vertices positions, but here Im working with edges, so no
            %   need to center vertices as far as I understand.
            [n,~] = size(obj.VD_vertices);
            %   by definition in the article:
            %{
            N = kron( ( speye(2*n) - (1/(2*n)).*ones((2*n),(2*n)) ), speye(3) );
            V = [obj.vertices_plus; obj.vertices_minus];
            V = V';
            centeredV = N*(V(:))';
            %}
            
            [A_planarity, p_planarity] = obj.getPlanarityProjection();
            [A_height, p_height] = obj.getHeightProjection();
            [A_parallel, p_parallel] = obj.getParallelityProjection();
            [A_offset, p_offset] = obj.getOffsetDirProjection();
            [A_length, p_length] = obj.getMinLengthProjection();
            
            A = [A_planarity; A_height; A_parallel; A_offset; A_length];
            p = [p_planarity; p_height; p_parallel; p_offset; p_length];
            
            newV = kron(speye(3), A) \ p(:);
            newV = reshape(newV, 2*n, 3);
            
            diff_plus = newV(1:n,:) - obj.vertices_plus;
            diff_minus = newV((n+1):end,:) - obj.vertices_minus;
            distance = trace(diff_plus*diff_plus') + trace(diff_minus*diff_minus');
            
            obj.vertices_plus = newV(1:n,:);
            obj.vertices_minus = newV((n+1):end,:);
        end
        
        function convergeConstraints(obj)
            for i=1:100
                d = obj. projectBeamConstraints();
                disp(['current distance: ' num2str(d)]);
            end
        end
        
        function [A, planarity_projection] = getPlanarityProjection(obj)
            %   Calculate the projection for the planarity constraint.
            %   For each row in adj. matrix, there will be exactly 3
            %   nonzero values. because adj. matrix is symmetric, let's
            %   operate only on lower triangle. for each edge (that defines
            %   a beam) for all beam edges, project it.
            [n,~] = size(obj.VD_vertices);
            %   for each edge in VD, there are 4 edges in the beam mesh.
            %   number of edges is (3/2)|VD_vertices|
            planarity_constraints = 6*n;
            planarity_projection = zeros(planarity_constraints, 3);
            
            A_i = repmat(1:planarity_constraints, 2, 1);
            A_i = A_i(:);
            A_j = zeros(2*planarity_constraints, 1);
            A_val = [ones(1, planarity_constraints); (-1).*ones(1, planarity_constraints)];
            A_val = A_val(:);
            
            current = 1;
            
            for i=1:n
                % for vertex i in VD_vertices:
                %   get neighbors to v_i from lower triangle of adj. mat.
                neighbors = find(obj.VD_adjacencyMatrix(i,1:i));
                
                for j = neighbors
                    [n_i, n_j] = obj.getBeamNormals(i, j);
                    [~, ~, e_ij] = obj.getMidVertices(i, j);
                    n_ij = cross( (n_i + n_j), e_ij);
                    
                    %   assumming the order of the edges is e_nrml_i,
                    %   e_nrml_j, e_ofst_pl, e_ofst_mi: I will assemble the
                    %   required A matrix rows for the projection:
                    %   the code below means: in the 'current' row, and
                    %   i-th column (for example) put the value (1 or -1).
                    %   e_nrml_i:
                    current_j_idx = 8*((current-1)/4);
                    A_j(current_j_idx + 1) = i; % v_i+
                    A_j(current_j_idx + 2) = n + i; % v_i-
                    %   e_nrml_j:
                    A_j(current_j_idx + 3) = j; % v_j+
                    A_j(current_j_idx + 4) = n + j; % v_j-
                    %   e_ofst_pl:
                    A_j(current_j_idx + 5) = j; % v_j+
                    A_j(current_j_idx + 6) = i; % v_i+
                    %   e_ofst_mi:
                    A_j(current_j_idx + 7) = n + j; % v_j-
                    A_j(current_j_idx + 8) = n + i; % v_i-
                    
                    BeamEdgesMat = obj.getBeamEdges(i,j);
                    for edge_idx=1:4
                        edge = BeamEdgesMat(:, edge_idx);
                        c = (edge' * n_ij)./(n_ij' * n_ij);
                        planarity_projection(current,:) = (edge - c.*n_ij)';
                        
                        current = current + 1;
                    end
                end
            end
            
            A = sparse(A_i(:), A_j, A_val(:), 6*n, 2*n);
        end
        
        function [A_height, height_projection] = getHeightProjection(obj)
            %   I will use the 1st option suggested: use one target height
            %   and a tolerance variable.
            %   each beam is going to yield 2 constraints
            [n,~] = size(obj.VD_vertices);
            %   there are 2|E| = 2*(3/2)|V| = 3|V| constraints
            height_constrinats = 3*n;
            height_projection = zeros(height_constrinats, 3);
            
            A_i = repmat(1:height_constrinats, 2, 1);
            A_i = A_i(:);
            A_j = zeros(2*height_constrinats, 1);
            A_val = [ones(1, height_constrinats); (-1).*ones(1, height_constrinats)];
            A_val = A_val(:);
            
            current = 1;
            
            for i=1:n
                %   get neighbors to v_i from lower triangle of adj. mat.
                neighbors = find(obj.VD_adjacencyMatrix(i,1:i));
                
                for j = neighbors
                    [n_i, n_j] = obj.getBeamNormals(i, j);
                    [~, ~, e_ij] = obj.getMidVertices(i, j);
                    
                    %   assumming the order of the edges is e_nrml_i,
                    %   e_nrml_j, e_ofst_pl, e_ofst_mi: I will assemble the
                    %   required A matrix rows for the projection:
                    %   the code below means: in the 'current' row, and
                    %   i-th column (for example) put the value (1 or -1).
                    %   e_nrml_i:
                    current_j_idx = 4*((current-1)/2);
                    A_j(current_j_idx + 1) = i; % v_i+
                    A_j(current_j_idx + 2) = n + i; % v_i-
                    %   e_nrml_j:
                    A_j(current_j_idx + 3) = j; % v_j+
                    A_j(current_j_idx + 4) = n + j; % v_j-
                    
                    nrmls = [n_i n_j];
                    
                    BeamEdgesMat = obj.getBeamEdges(i,j);
                    for idx=1:2
                        c = (nrmls(:,idx)' * e_ij)./(e_ij' * e_ij);
                        h_ij = 2*(nrmls(:,idx) - c.*e_ij);
                        h_ij_len = norm(h_ij);
                        if (norm(h_ij_len - obj.beamHeight) > obj.heightTolerance)
                            if (h_ij_len > obj.beamHeight)
                                constraint = ((obj.beamHeight+obj.heightTolerance)/h_ij_len).*BeamEdgesMat(:,idx);
                            else
                                constraint = ((obj.beamHeight-obj.heightTolerance)/h_ij_len).*BeamEdgesMat(:,idx);
                            end
                        else
                            constraint = BeamEdgesMat(:,idx);
                        end
                        height_projection(current,:) = constraint;
                        current = current + 1;
                    end
                end
            end
            A_height = sparse(A_i(:), A_j, A_val(:), 3*n, 2*n);
        end
        
        function [A_parallel, p_parallel] = getParallelityProjection(obj)
            [n,~] = size(obj.VD_vertices);
            %   there are 2|E| = 2*(3/2)|V| = 3|V| constraints
            parallel_constrinats = 3*n;
            p_parallel = zeros(parallel_constrinats, 3);
            
            A_i = repmat(1:parallel_constrinats, 2, 1);
            A_i = A_i(:);
            A_j = zeros(2*parallel_constrinats, 1);
            A_val = [ones(1, parallel_constrinats); (-1).*ones(1, parallel_constrinats)];
            A_val = A_val(:);
            
            current = 1;
            
            for i=1:n
                %   get neighbors to v_i from lower triangle of adj. mat.
                neighbors = find(obj.VD_adjacencyMatrix(i,1:i));
                
                for j = neighbors
                    [~, ~, e_ij] = obj.getMidVertices(i, j);
                    BeamEdgesMat = obj.getBeamEdges(i,j);
                    %   assumming the order of the edges is e_nrml_i,
                    %   e_nrml_j, e_ofst_pl, e_ofst_mi: I will assemble the
                    %   required A matrix rows for the projection:
                    current_j_idx = 4*((current-1)/2);
                    A_j(current_j_idx + 1) = j; % v_j+
                    A_j(current_j_idx + 2) = i; % v_i+
                    %   e_ofst_mi:
                    A_j(current_j_idx + 3) = n + j; % v_j-
                    A_j(current_j_idx + 4) = n + i; % v_i-
                    
                    for idx=3:4
                        c = (e_ij' * BeamEdgesMat(:,idx))./(e_ij' * e_ij);
                        p_parallel(current,:) = c.*e_ij;
                        current = current + 1;
                    end
                end
            end
            A_parallel = sparse(A_i(:), A_j, A_val(:), 3*n, 2*n);
        end
        
        function [A_offset, p_offset] = getOffsetDirProjection(obj)
            [n,~] = size(obj.VD_vertices);
            %   there are 2|E| = 2*(3/2)|V| = 3|V| constraints
            offset_constrinats = 3*n;
            p_offset = zeros(offset_constrinats, 3);
            
            A_i = repmat(1:offset_constrinats, 2, 1);
            A_i = A_i(:);
            A_j = zeros(2*offset_constrinats, 1);
            A_val = [ones(1, offset_constrinats); (-1).*ones(1, offset_constrinats)];
            A_val = A_val(:);
            
            current = 1;
            
            for i=1:n
                %   get neighbors to v_i from lower triangle of adj. mat.
                neighbors = find(obj.VD_adjacencyMatrix(i,1:i));
                
                for j = neighbors
                    BeamEdgesMat = obj.getBeamEdges(i,j);
                    %   assumming the order of the edges is e_nrml_i,
                    %   e_nrml_j, e_ofst_pl, e_ofst_mi: I will assemble the
                    %   required A matrix rows for the projection:
                    current_j_idx = 4*((current-1)/2);
                    A_j(current_j_idx + 1) = i; % v_i+
                    A_j(current_j_idx + 2) = n + i; % v_i-
                    %   e_ofst_mi:
                    A_j(current_j_idx + 3) = j; % v_j+
                    A_j(current_j_idx + 4) = n + j; % v_j-
                    
                    dir_i = (BeamEdgesMat(:,1)'*(obj.verticesNormals(i,:)'./norm(obj.verticesNormals(i,:)))).* obj.verticesNormals(i,:);
                    p_offset(current,:) = dir_i;
                    dir_j = (BeamEdgesMat(:,2)'*(obj.verticesNormals(j,:)'./norm(obj.verticesNormals(j,:)))).* obj.verticesNormals(j,:);
                    p_offset(current+1,:) = dir_j;
                    current = current + 2;
                end
            end
            A_offset = sparse(A_i(:), A_j, A_val(:), 3*n, 2*n);
        end
        
        function [A_length, p_length] = getMinLengthProjection(obj)
            [n,~] = size(obj.VD_vertices);
            %   there are 2|E| = 2*(3/2)|V| = 3|V| constraints
            length_constrinats = 3*n;
            p_length = zeros(length_constrinats, 3);
            
            A_i = repmat(1:length_constrinats, 2, 1);
            A_i = A_i(:);
            A_j = zeros(2*length_constrinats, 1);
            A_val = [ones(1, length_constrinats); (-1).*ones(1, length_constrinats)];
            A_val = A_val(:);
            
            current = 1;
            
            for i=1:n
                %   get neighbors to v_i from lower triangle of adj. mat.
                neighbors = find(obj.VD_adjacencyMatrix(i,1:i));
                
                for j = neighbors
                    BeamEdgesMat = obj.getBeamEdges(i,j);
                    %   assumming the order of the edges is e_nrml_i,
                    %   e_nrml_j, e_ofst_pl, e_ofst_mi: I will assemble the
                    %   required A matrix rows for the projection:
                    current_j_idx = 4*((current-1)/2);
                    A_j(current_j_idx + 1) = j; % v_j+
                    A_j(current_j_idx + 2) = i; % v_i+
                    %   e_ofst_mi:
                    A_j(current_j_idx + 3) = n + j; % v_j-
                    A_j(current_j_idx + 4) = n + i; % v_i-
                    
                    e_pl_len = norm(BeamEdgesMat(:,3));
                    if (e_pl_len < obj.minEdgeLengthThreshold)
                        p_length(current,:) = (obj.minEdgeLengthThreshold / e_pl_len).* (BeamEdgesMat(:,3)');
                    else
                        p_length(current,:) = BeamEdgesMat(:,3)';
                    end
                    
                    e_mi_len = norm(BeamEdgesMat(:,4));
                    if (e_mi_len < obj.minEdgeLengthThreshold)
                        p_length(current+1,:) = (obj.minEdgeLengthThreshold / e_mi_len).* (BeamEdgesMat(:,4)');
                    else
                        p_length(current+1,:) = BeamEdgesMat(:,4)';
                    end
                    current = current + 2;
                end
            end
            A_length = sparse(A_i(:), A_j, A_val(:), 3*n, 2*n);
        end
        
        function showBeamMesh(obj)
            %   Visualize the beam mesh:
            [V_size,~] = size(obj.VD_adjacencyMatrix);
            
            figure;
            axis equal;
            for i=1:V_size
                vertices = find(obj.VD_adjacencyMatrix(i,1:i));
                if(isempty(vertices))
                    continue;
                end
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