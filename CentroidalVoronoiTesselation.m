classdef CentroidalVoronoiTesselation < handle
    % Centroidal Voronoi Tesselation for a 3d mesh
    %  
    
    properties
        mesh
        sites
        numberOfSites
        cells
        metrics
        metricTensors
        voronoiVertices
        voronoiEdges
    end
    
    methods
        function obj = CentroidalVoronoiTesselation(mesh, initial_x)
            %   C'tor for the class. 
            %   params:
            %       * mesh: the mesh (of class Mesh) of the shape
            %       * initial_x: the initial sites to begin with. This
            %       should be a matrix of size k X 3 when k is number
            %       of initial sites.
            %       * metric_tensors: the metric of each site. default will
            %       be the identity matrix I. This should be a matrix of
            %       size k X 3 X 3 when k is number of initial sites.
            obj.mesh = mesh;
            obj.sites = initial_x;
            obj.numberOfSites = length(obj.sites);
            %   the 'cells' will be an array of size |V| such that every
            %   entry will represent the index of the voronoi cell of the
            %   vertex:
            obj.cells = zeros(obj.mesh.dimensions(1), 1);
            
            %   Initial value of Metric Tensor Matrices - Identity matrix I
            obj.metrics = zeros(3,3,obj.numberOfSites);
            obj.metricTensors = cell(obj.numberOfSites,1);
            for i=1:obj.numberOfSites
                obj.metrics(:,:,i) = eye(3);
                obj.metricTensors{i} = speye(3);
            end
            
            obj.runIterations();
            
        end
        
        function iteration_dist = calculateCellsAndSites(obj)
            %   calculate the cells based on current values of sites.
            %   during that proccess, calculate the mean value of each new
            %   cell, and at the end update the sites to that mean.
            
            iteration_dist = 0;
            
            %   efficient version:
            %   calculate in a matrix p-x_i for all p and all x_i such that
            %   each row is (x,y,z) triplets of p-x_i for all i. and
            %   different rows are different p. then, stack them all to one
            %   vector.
            V_size = obj.mesh.dimensions(1);
            
            p = obj.mesh.vertices';
            p = repmat(p, obj.numberOfSites, 1);
            p = p(:);
            allsites = obj.sites';
            allsites = allsites(:);
            allsites = repmat(allsites, V_size, 1);
            distanceVec = p-allsites;
            
            %   make a matrix whose diagonal is composed out of the metric
            %   tensors M_i, and repeat it |V| times so it will match the
            %   vector p-x_i we made.
            M = blkdiag(obj.metricTensors{:});
            %M_cell = repmat({M}, V_size, 1);
            %M_all = blkdiag(M_cell{:});
            M_all = kron(speye(V_size),M);
            
            res = M_all * distanceVec;
            res_size = length(res);
            
            % put mid-result to a diagonal matrix:
            diagRes = sparse(1:res_size, 1:res_size, res);
            % now let's transofrm it to a matrix of size (3k|V| X k|V|) so each vector
            % has 3 ones for the current triplet. so after multiplication
            % we'll have k distances for each point p. then we'll get the
            % minimum.
            kronMat = kron(speye(obj.numberOfSites * V_size), sparse([1;1;1]));
            
            res = diagRes*kronMat;
            
            %finish computation:
            res = distanceVec' * res;
            res = reshape(res,obj.numberOfSites, V_size);
            [~,new_cells] = min(res);
            new_cells = new_cells';
            
            %   now calc new sites: 
            %   for this, let's create a matrix so vector k of it will contain
            %   1 in indexes of p that clsoset to site k and 0 for others.
            %   then multiply transpose of p matrix by it. it will result
            %   in 3 X k matrix where each column is sum of all points in
            %   its cell.
            summationMatrix = sparse(1:V_size, new_cells, ones(V_size,1));
            cellSum = obj.mesh.vertices' * summationMatrix;
            cellSize = sum(summationMatrix);
            for i=1:obj.numberOfSites
                if cellSize ~= 0
                    cellSum(:,i) = cellSum(:,i) ./ cellSize(i);
                    iteration_dist = iteration_dist + norm(obj.sites(i,:) - cellSum(:,i)');
                end
            end
            
            new_sites = cellSum';
            
            %   assign the new cells and sites to the object:
            obj.cells = new_cells;
            obj.sites = new_sites;
        end
        
        function calculateMetrics(obj)
            %   call this function only after at least one iteration has
            %   been done! This function calculates a new metric tensors
            %   based on the cells that has been calculated earlier.
            new_metrics = zeros(3, 3, obj.numberOfSites);
            
            %   compute C_i:
            for i=1:length(obj.cells)
                %   calculate C_i, using the sum of all points in its cell.
                %   then, calculate eigen decomposition, and M_i.
                cell_idx = obj.cells(i);
                p = obj.mesh.vertices(i,:);
                site = obj.sites(cell_idx,:);
                new_metrics(:,:,cell_idx) = new_metrics(:,:,cell_idx) + ( (p - site)' * (p - site) );
            end
            
            for k=1:obj.numberOfSites
                %   compute M_i using C_i:
                C_i = new_metrics(:,:,k);
                [Q,S,Q_T] = svd(C_i);
                S_inv = diag(1 ./ diag(S));
                new_metrics(:,:,k) = S(3,3) .* (Q*S_inv*Q_T');
                obj.metricTensors{k} = new_metrics(:,:,k);
            end
            
            obj.metrics = new_metrics;
        end
        
        function runIterations(obj)
            %   do first iteration, then compute metrics, and then do the
            %   rest of the iterations until convergence.
            w = waitbar(0,'Initializing...');
            
            obj.calculateCellsAndSites();
            obj.calculateMetrics();
            
            numOfIterations = 10;
            
            waitbar(0,w,'current iteration');
            
            i=1;
            while i < numOfIterations
                waitbar(i/numOfIterations);
                %profile on;
                iteration_dist = obj.calculateCellsAndSites();
                obj.calculateMetrics();
                %profile off;
                %profile viewer;
                
                disp(['iteration ' num2str(i) ...
                ' distance from centroid = ' num2str(iteration_dist)]);
                
                if iteration_dist == 0
                    break;
                end
                i = i+1;
                if i == numOfIterations
                    prompt = 'Do you want more? Y/N [Y]: ';
                    str = input(prompt,'s');
                    if isempty(str)
                        str = 'Y';
                    end
                    if str == 'Y'
                        numOfIterations = numOfIterations + 100;
                    end
                end
            end
            
            %DT = delaunayTriangulation(obj.sites);
            %[VoronoiDiagram, VoronoiFaces] = voronoiDiagram(DT);
            [obj.voronoiVertices, obj.voronoiEdges] = obj.findVoronoiVertices();
            
            close(w);
        end
        
        function [V,V2] = findVoronoiVertices(obj)
            N_F = obj.mesh.dimensions(2);
            V = zeros(N_F,3);
            V2 = zeros(N_F,3);
            N_voronoi = 0;
            n_v2 = 0;
            for i=1:N_F
                %   iterate over all faces, and find those who have
                %   vertices from 3 different cells => vertices of the
                %   Voronoi diagram will be the center of those faces.
                %    This can be done using the adjacency matrix.
                neighbors = obj.mesh.faces(i,:);
                neighbors_cells = unique(obj.cells(neighbors));
                if (length(neighbors_cells) >= 3)
                    %   this is a vertex of the voronoy diagram:
                    N_voronoi = N_voronoi + 1;
                    V(N_voronoi,:) = obj.mesh.getFaceBarycentricPoint(i,[1/3,1/3,1/3]);
                end
                if (length(neighbors_cells) >= 2)
                    %   voronoi edge
                    n_v2 = n_v2 + 1;
                    V2(n_v2,:) = obj.mesh.getFaceBarycentricPoint(i,[1/3,1/3,1/3]);
                end
            end
            
            % remove redundant last lines that are 0s:
            V = V(1:N_voronoi,:);
            V2 = V2(1:n_v2,:);
        end
        
        function faces = findVoronoiFaces(obj) 
            % the method is to find the closest point to the site in its
            % cell, then using the normal and
        end
        
    end
    
end

