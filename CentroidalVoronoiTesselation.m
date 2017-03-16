classdef CentroidalVoronoiTesselation < handle
    % Centroidal Voronoi Tesselation for a 3d mesh
    %  I denote K as the number of sites of the CVT.
    
    properties
        mesh
        sites
        numberOfSites
        cells
        metrics
        metricTensors
        voronoiVertices
        voronoiVertexNormals
        %   matrix where each row i is the 3 neighbor cells of vertex i.
        voronoiVertexNeighborCells
        
        %   a martix with row for each cell, and each row has 1 in indices
        %   of vertices that assocciated with the cell:
        voronoiCellVertices
        
        voronoiAdjMatrix
        
        %   3d coordinates in cetroid of faces that are on edge of voronoi
        %   diagram.
        voronoiEdgeTriangles
        
        %   a matrix of size |V| x K where each column j represents which
        %   points p are assocciated with site j.
        summationMatrix
        
        %   a matrix of size K x |V|, where each column i is distances of
        %   point p_i from each site. so d_i,j is the distance of vertex j
        %   from site i. distances are measured as described in article:
        distancesMatrix
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
            
            success = obj.runIterations();
            if (~success)
                return;
            end
            
            disp('Extracting Cells...');
            obj.getVoronoiCells();
            obj.findVoronoiVertices();
            obj.findVoronoiFaces();
            
            
            isValid = obj.validateCVT();
            if(~isValid)
                disp('Trying to repair CVT...');
                obj.repairCVT();
                if(~obj.validateCVT())
                    disp('CVT is BAD!!! run again with different conditions.');
                end
            end
        end
        
        function iteration_dist = calculateCellsAndSites(obj)
            %   calculate the cells based on current values of sites.
            %   during that proccess, calculate the mean value of each new
            %   cell, and at the end update the sites to that mean.
            
            iteration_dist = 0;
            
            %   efficient version:
            %   calculate in a matrix p-x_i for all p and all x_i such that
            %   each row is (x,y,z) triplets of p-x_i for all i. and
            %   different rows are different p. 
            V_size = obj.mesh.dimensions(1);
            
            p = obj.mesh.vertices';
            p = repmat(p, obj.numberOfSites, 1);
            allsites = obj.sites';
            allsites = allsites(:);
            allsites = repmat(allsites, 1, V_size);
            distanceVec = p-allsites;
            
            %   make a matrix whose diagonal is composed out of the metric
            %   tensors M_i:
            M = blkdiag(obj.metricTensors{:});
            
            %   calculate M*(p-xi). this is half the way. then stack all
            %   columns:
            res = M * distanceVec;
            res = res(:);
            
            % now let's transofrm it to a matrix of size (3k|V| X k|V|) so each vector
            % has 3 entries for the current triplet. so after multiplication
            % we'll have k distances for each point p. then we'll get the
            % minimum.
            triMat = repmat(1:(obj.numberOfSites * V_size),3,1);
            rowNum = (3 * obj.numberOfSites * V_size);
            res = sparse(1:rowNum, triMat(:), res, rowNum, obj.numberOfSites * V_size);
            
            %finish computation:
            res = distanceVec(:)' * res;
            res = reshape(res,obj.numberOfSites, V_size);
            obj.distancesMatrix = res;
            [~,new_cells] = min(res);
            new_cells = new_cells';
            
            %   now calc new sites: 
            %   for this, let's create a matrix so vector k of it will contain
            %   1 in indexes of p that clsoset to site k and 0 for others.
            %   then multiply transpose of p matrix by it. it will result
            %   in 3 X k matrix where each column is sum of all points in
            %   its cell.
            obj.summationMatrix = sparse(1:V_size, new_cells, ones(V_size,1));
            cellSum = obj.mesh.vertices' * obj.summationMatrix;
            cellSize = sum(obj.summationMatrix);
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
            %   summationMatrix (|V| x k) contain in the i-th column 1
            %   in indices of points that belong to cell i.
            %V_size = obj.mesh.dimensions(1);
            %obj.summationMatrix = sparse(1:V_size, obj.cells, ones(V_size,1));
            %   calc matrix whith |V| columns, where each col. i is the
            %   site coordinates of the site that p_i is in.
            p_sites = obj.sites' * obj.summationMatrix';
            distanceMat = obj.mesh.vertices' - p_sites;
            
            for k=1:obj.numberOfSites
                idx = find(obj.summationMatrix(:,k));
                if (isempty(idx))
                    new_metrics(:,:,k) = speye(3);
                    continue;
                end
                cellPointsDiff = distanceMat(:, idx);
                new_metrics(:,:,k) = cellPointsDiff * cellPointsDiff';
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
        
        function success = runIterations(obj)
            %   do first iteration, then compute metrics, and then do the
            %   rest of the iterations until convergence.
            success = true;
            w = waitbar(0,'Initializing...');
            
            %obj.calculateCellsAndSites();
            %obj.calculateMetrics();
            
            numOfIterations = 10;
            
            waitbar(0,w,'current iteration');
            
            i=1;
            while i < numOfIterations
                waitbar(i/numOfIterations);
                
                iteration_dist = obj.calculateCellsAndSites();
                obj.calculateMetrics();
                
                disp(['iteration ' num2str(i) ...
                ' distance from centroid = ' num2str(iteration_dist)]);
                
                if iteration_dist == 0
                    if (i == 1)
                        close(w);
                        disp('bad initial sites, please try again.');
                        success = false;
                        return;
                    end
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
            
            close(w);
        end
        
        function [V,V2] = findVoronoiVertices(obj)
            N_F = obj.mesh.dimensions(2);
            V = zeros(N_F,3);
            V_neighbors = zeros(N_F,3);
            V2 = zeros(N_F,3);
            N_voronoi = 0;
            n_v2 = 0;
            
            facesNormals = obj.mesh.getFacesNormals();
            obj.voronoiVertexNormals = zeros(N_F,3);
            
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
                    V_neighbors(N_voronoi,:) = neighbors_cells;
                    obj.voronoiVertexNormals(N_voronoi,:) = facesNormals(i,:);
                end
                if (length(neighbors_cells) == 2)
                    %   voronoi edge
                    n_v2 = n_v2 + 1;
                    V2(n_v2,:) = obj.mesh.getFaceBarycentricPoint(i,[1/3,1/3,1/3]);
                end
            end
            
            % remove redundant last lines that are 0s:
            V = V(1:N_voronoi,:);
            obj.voronoiVertices = V;
            obj.voronoiVertexNeighborCells = V_neighbors(1:N_voronoi,:);
            obj.voronoiVertexNormals = obj.voronoiVertexNormals(1:N_voronoi,:);
            V2 = V2(1:n_v2,:);
            obj.voronoiEdgeTriangles = V2;
        end
        
        function faces = findVoronoiFaces(obj) 
            % the method is to find the closest point to the site in its
            % cell, then using the normal and the cell's site, use them as
            % reference point and define order clockwise of the vertices.
            %   first, let's generate for each cell a list of it's
            %   vertices:
            [numOfV,~] = size(obj.voronoiVertices);
            trimat = repmat(1:numOfV,3,1);
            
            %   matrix where each row i is 1 in indices of vetrices
            %   assocciated with cell i, and 0 otherwise.
            obj.voronoiCellVertices = sparse(obj.voronoiVertexNeighborCells', ...
                trimat, ones(3*numOfV,1));
            
            obj.voronoiAdjMatrix = sparse(numOfV,numOfV);
            for k=1:obj.numberOfSites
                %   iterate over couples of vertices of site k, and for 2
                %   vertices who share 2 common cells, put an edge between
                %   them. (each cell is at most 10 vertices from what I
                %   experienced so shouldnt be a problem).
                vertices = find(obj.voronoiCellVertices(k,:));
                n = length(vertices);
                
                for i=1:(n-1)
                    for j=i+1:n
                        common = unique([obj.voronoiVertexNeighborCells(vertices(i),:), ...
                            obj.voronoiVertexNeighborCells(vertices(j),:)]);
                        if (length(common) == 4)
                            %   set an edge in adj Matrix:
                            obj.voronoiAdjMatrix(vertices(i), vertices(j)) = 1;
                            obj.voronoiAdjMatrix(vertices(j), vertices(i)) = 1;
                        end
                    end
                end
                
            end
            
        end
        
        function faces = findVoronoiFaces_omri(obj)
            % ### This is not working well!
            %   Try using omri's idea - using angles!
            faces = cell(obj.numberOfSites);
            
            % the method is to find the closest point to the site in its
            % cell, then using the normal and the cell's site, use them as
            % reference point and define order clockwise of the vertices.
            %   first, let's generate for each cell a list of it's
            %   vertices:
            [numOfV,~] = size(obj.voronoiVertices);
            trimat = repmat(1:numOfV,3,1);
            
            %   matrix where each row i is 1 in indices of vetrices
            %   assocciated with cell i, and 0 otherwise.
            obj.voronoiCellVertices = sparse(obj.voronoiVertexNeighborCells', ...
                trimat, ones(3*numOfV,1));
            
            obj.voronoiAdjMatrix = sparse(numOfV,numOfV);
            
            for k=1:obj.numberOfSites
                %   use first vertex as reference point (direction on
                %   plane). so for first vertex the angel is 0. so
                %   calculate angles of vertiex 2-n:
                vertices = find(obj.voronoiCellVertices(k,:));
                n = length(vertices);
                angles = zeros(n,1);
                v_reference = (obj.voronoiVertices(vertices(1),:) - obj.sites(k,:))';
                for v_idx=2:n
                    v_current = (obj.voronoiVertices(vertices(v_idx),:) - obj.sites(k,:))';
                    angles(v_idx) = obj.getAngleBetweenVectors(v_reference, v_current);
                end
                
                Order = sortrows([angles vertices(:)]);
                faces{k} = Order(:,2);
                for v_idx = 1:(n-1)
                    %   set an edge in adj Matrix:
                    nxt = v_idx+1;
                    obj.voronoiAdjMatrix(Order(v_idx, 2), Order(nxt, 2)) = 1;
                    obj.voronoiAdjMatrix(Order(nxt, 2), Order(v_idx,2)) = 1;
                end
                %   now for first and last:
                obj.voronoiAdjMatrix(Order(n, 2), Order(1, 2)) = 1;
                obj.voronoiAdjMatrix(Order(1, 2), Order(n,2)) = 1;
            end
        end
        
        function angle = getAngleBetweenVectors(v1, v2)
            n1 = norm(v1);
            n2 = norm(v2);
            if (n1 == 0 && n2 == 0)
                angle = 0;
                return;
            end
            angle = acos( dot(v1,v2) ./ (n1.*n2) );
        end
        
        function N = getFacesNormals(obj)
            %   For each voronoi cell calculate it's normal.
            %   for now, I'll do it by calculating the mean of all normals
            %   in the cell.
            N_V = obj.mesh.getVertexNormals()';
            V_size = obj.mesh.dimensions(1);
            obj.summationMatrix = sparse(1:V_size, obj.cells, ones(V_size,1));
            
            N = N_V * obj.summationMatrix;
            cellSize = sum(obj.summationMatrix);
            N = N ./ repmat(cellSize,3,1);
        end
        
        function [d,r] = getRepresentativeVertices(obj)
            %   compute for each cell its representative vertex, that is,
            %   the vertex whose closest to the site of the cell. also,
            %   include the distance from x_i.
            [d,r] = min(obj.distancesMatrix,[],2);
            unique_r = unique(r);
            if (length(unique_r) < length(r))
                disp('discarding sites!');
                r = unique_r;
            end
        end
        
        function cells = getVoronoiCells(obj)
            %   Compute cells according to 2nd article. (dijkstra like). I
            %   will ignore k-d tree for simplicity. I will implement a
            %   slight changed version of the algorithm: I will assume the
            %   distances from real sites were computed in
            %   obj.distancesMatrix. I will add vertices to queue only if a
            %   neighbor was assigned to a cell. This will ensure
            %   connectivity.
            [d,r] = obj.getRepresentativeVertices();
            
            %   3-col matrix. 1st col is distance from site, 2nd col is
            %   index of point p, and 3-rd is the site from which p is
            %   measured.
            Q = zeros((2*obj.mesh.dimensions(3)),3);
            Q(1:obj.numberOfSites,:) = [d r (1:obj.numberOfSites)'];
            Q(1:obj.numberOfSites,:) = sortrows(Q(1:obj.numberOfSites,:));
            occupation = obj.numberOfSites;
            start = 1;
            
            cells = zeros(obj.mesh.dimensions(1),1);
            
            [r,c] = find(obj.mesh.AdjMatrix);
            
            while (~isempty(Q))
                %   pop from queue:
                current = Q(start,:);
                start = start + 1;
                %occupation = occupation - 1;
                
                if(start >= occupation)
                    break;
                end
                
                %   assign to cell (if hasn't been assigned yet):
                if (cells(current(2)) == 0)
                    cells(current(2)) = current(3);
                    
                    %   add neighbors:
                    %neighbors = find(obj.mesh.AdjMatrix(current(2),:));
                    %efficient try:
                    neighbors = find(r == current(2));
                    neighbors = c(neighbors)';
                    for neighbor=neighbors
                        %   insert to Q if not assigned to cell yet:
                        if (cells(neighbor) == 0)
                            dist = obj.distancesMatrix(current(3), neighbor); 
                            row = [dist neighbor current(3)];
                            
                            place = find(Q(start:occupation,1) > dist, 1);
                            place = place + start -1;
                            occupation = occupation + 1;
                            Q(occupation,:) = row;
                            
                            if (~isempty(place))
                                Q(start:occupation,:) = [Q(start:(place-1),:);...
                                                        Q(occupation,:); ...
                                                        Q(place:(occupation-1),:)];
                            end
                        end
                    end
                end
                
            end
            
            obj.cells = cells;
        end
        
        function discardInvalidCells(obj)
            %   try removing cells with less than 3 voronoi vertices, and
            %   cells with less than some kind of threshold of vertices.
            cellsVertices = sum(obj.voronoiCellVertices, 2);
            toDiscard = find(cellsVertices < 3);
            discarded = length(toDiscard);
            obj.numberOfSites = obj.numberOfSites - discarded;
            if (discarded < 1)
                %   try discarding the smallest site:
                %[~,toDiscard] = min(cellsVertices);
                toDiscard = find(cellsVertices <= 3);
                discarded = length(toDiscard);
                if (discarded > 0)
                    k = ceil(discarded/5);
                    toDiscard = datasample(toDiscard,k);
                    obj.numberOfSites = obj.numberOfSites - k;
                else
                    return;
                end
            end
            for c=toDiscard
                obj.sites(c,:) = [];
                obj.voronoiCellVertices(c,:) = [];
                obj.distancesMatrix(c,:) = [];
                obj.summationMatrix(:,c) = [];
            end
        end
        
        function repairCVT(obj)
            %   try discarding sites in hope to repair:
            %threshold = obj.numberOfSites - 1;
            cur_num_of_sites = obj.numberOfSites;
            tries = 1;
            while (~obj.validateCVT())
                disp(['Repair try #' num2str(tries)]);
                obj.discardInvalidCells();
                obj.getVoronoiCells();
                obj.findVoronoiVertices();
                obj.findVoronoiFaces();
                change = cur_num_of_sites - obj.numberOfSites;
                if (change == 0)
                    break
                end
                tries = tries + 1;
                cur_num_of_sites = obj.numberOfSites;
            end
        end
        
        function isValid = validateCVT(obj)
            %   In some cases the CVT isn't correct due to bad regions.
            %   This method verifies that the adjacency matrix is good and
            %   that means the connectivity is good.
            n = nnz(obj.voronoiAdjMatrix);
            %   number of non zero elements in adj matrix should be equal
            %   to 2*|E| = 3*|V| = n.
            isValid = true;
            [v_size,~] = size(obj.voronoiVertices);
            if (n/3 ~= v_size)
                isValid = false;
            end
            
            isSymmetric = obj.voronoiAdjMatrix - obj.voronoiAdjMatrix';
            if (norm(isSymmetric(:)) ~= 0)
                isValid = false;
            end
            
            if (trace(obj.voronoiAdjMatrix) ~= 0)
                isValid = false;
            end
        end
        
        function showResults(obj)
            obj.mesh.showMesh(obj.cells);
            alpha(0.85);
            hold on;
            
            V_v = obj.voronoiVertices;
            V_e = obj.voronoiEdgeTriangles;
            scatter3(obj.sites(:,1),obj.sites(:,2),obj.sites(:,3),100,'filled','yellow');
            scatter3(V_e(:,1),V_e(:,2),V_e(:,3),25,'filled','cyan');
            scatter3(V_v(:,1),V_v(:,2),V_v(:,3),50,'filled','red');
            
            % show normals:
            %N = obj.getFacesNormals();
            %quiver3(obj.sites(:,1), obj.sites(:,2), obj.sites(:,3) ...
            %    ,N(1,:)', N(2,:)', N(3,:)','-k');
        end
        
        function showPolygon(obj, mixed)
            %   show polygonal approximation using adj matrix:
            %   the mixed argument lets you choose weater to show the
            %   original mesh also or not
            if (nargin < 2)
                mixed = false;
            end
            if (mixed)
                obj.mesh.showMesh(obj.cells);
                alpha(0.5);
                hold on;
            else
                figure;
            end
            [numOfV,~] = size(obj.voronoiVertices);
            for i=1:numOfV
                vertices = find(obj.voronoiAdjMatrix(i,1:i));
                if (isempty(vertices))
                    continue;
                end
                for j=vertices
                    edge = [obj.voronoiVertices(i,:)' obj.voronoiVertices(j,:)'];
                    plot3(edge(1,:),edge(2,:),edge(3,:),'-k','LineWidth',3);
                    hold on;
                end
            end
        end
        
    end %   end of methods
    
end %   end of class

