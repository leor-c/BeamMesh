classdef Mesh < handle
    % A class that represent a mesh triangulation.
    %   Detailed explanation goes here
    
    properties
        vertices
        faces
        dimensions
        AdjMatrix
        WeightedAdjMatrix
        f_areas
        v_areas
        f_normals
        grad
    end
    
    methods
        function obj = Mesh(mesh_off_file)
           %    To initialize a mesh you must pass a .off file.
           if (nargin == 0)
               error('To initialize a mesh you must pass a .off file.');
           end
           obj.load_off_file(mesh_off_file);
           obj.faceAreas();
           obj.vertexAreas();
        end
        
        function load_off_file( obj, filename )
            %load a .off file and make a shared vertex data structure out of it
            %   .off files are structured like so: https://en.wikipedia.org/wiki/OFF_%28file_format%29
            %first, lets read the first line witch is always 'OFF'
            f = fopen(filename);
            fgetl(f);
            A = fscanf(f,'%d %d %d',[1,3]);
            %   numOfVertices, numOfFaces, numOfEdges
            sizes = A;
            numOfV = A(1);
            numOfF = A(2);
            obj.vertices = fscanf(f,'%f', [3, numOfV]);
            obj.vertices = obj.vertices';
            obj.faces = fscanf(f,'%u', [4, numOfF]); % this is for triangles! assume first number is 3!
            obj.faces = obj.faces';
            obj.faces = obj.faces(:,2:end);   % cut the first culumn since it only contains redundant 3
            obj.faces = obj.faces + ones(size(obj.faces));
            row = zeros(1,3*numOfF);
            col = zeros(1,3*numOfF);
            len = zeros(1,3*numOfF);
            for i=1:numOfF
                row(3*(i-1)+1) = obj.faces(i,1);
                col(3*(i-1)+1) = obj.faces(i,2);
                row(3*(i-1)+2) = obj.faces(i,1);
                col(3*(i-1)+2) = obj.faces(i,3);
                row(3*(i-1)+3) = obj.faces(i,2);
                col(3*(i-1)+3) = obj.faces(i,3);
                l1 = obj.vertices(obj.faces(i,2),:) - obj.vertices(obj.faces(i,1),:);
                l2 = obj.vertices(obj.faces(i,3),:) - obj.vertices(obj.faces(i,1),:);
                l3 = obj.vertices(obj.faces(i,3),:) - obj.vertices(obj.faces(i,2),:);
                len(3*(i-1)+1) = norm(l1);
                len(3*(i-1)+2) = norm(l2);
                len(3*(i-1)+3) = norm(l3);
            end
            pos = [row' col'];
            pos = unique(sort(pos,2),'rows');
            g = graph(pos(:,1),pos(:,2));
            if sizes(3) == 0
               sizes(3) = height(g.Edges); 
            end
            obj.dimensions = sizes;
            %   regular adjacency matrix:
            obj.AdjMatrix = adjacency(g);
            %   create adjacency matrix where each cell includes the weight of edge
            %   (i,j) for positive weight function (length in this case):
            obj.WeightedAdjMatrix = sparse(row,col,len,numOfV,numOfV);
            fclose(f);
            len_v = length(obj.vertices); len_f = length(obj.faces);
            if obj.dimensions(1) ~= len_v || obj.dimensions(2) ~= len_f
                disp('Wrong size in .off file!');
                obj.dimensions(1) = len_v;
                obj.dimensions(2) = len_f;
            end 
        end
        
        function save_off_file( obj, filename )
        %saves a mesh in a .off file
            f=fopen(filename,'w');
            fprintf(f,'OFF\n');
            fprintf(f,'%d %d %d\n', obj.dimensions);
            fprintf(f,'%.17f %.18f %.18f\n',obj.vertices');
            fprintf(f,'%d   %d %d %d\n',( obj.faces - ones(size(obj.faces)) )');
            fclose(f);
        end
        
        function [genus, numBoundaryEdges]=topoData( obj )
            %	this function calculates the genus of a mesh and the number of boundary
            %	edges of a mesh.
                %   calculate adjacency matrix where each cell i,j counts
                %   the number of faces adjacent to edge i,j.
                %   for that purpose, we will define that an edge goes from lower vertex
                %   index to higher vertex index
                numberOfFaces = obj.dimensions(2);
                rows = zeros(1,3*numberOfFaces);
                cols = zeros(1,3*numberOfFaces);

                for i=1:numberOfFaces
                    %   for every face add 3 edges of that face.
                    rows(3*(i-1)+1)=obj.faces(i,1);
                    cols(3*(i-1)+1)=obj.faces(i,2);
                    rows(3*(i-1)+2)=obj.faces(i,1);
                    cols(3*(i-1)+2)=obj.faces(i,3);
                    rows(3*(i-1)+3)=obj.faces(i,2);
                    cols(3*(i-1)+3)=obj.faces(i,3);
                end

                %   lets keep the edges (indices of 2 vertices) in matrix:
                pos = [rows' cols'];
                %   sort each row such that each edge is from lower vertex to higher:
                sorted_pos = sort(pos, 2);
                unique_pos = unique(sorted_pos, 'rows');
                numberOfEdges = length(unique_pos);
                numBoundaryEdges = 0;

                %   count edges that have only 1 adjacent face:
                for i=1:length(unique_pos(1,:))
                    relevant_rows = find(sorted_pos == unique_pos(i,1));
                    relevant_cols = find(sorted_pos == unique_pos(i,2));
                    if length(intersect(relevant_rows,relevant_cols) ) == 1
                       numBoundaryEdges = numBoundaryEdges + 1;
                    end
                end

                numberOfVertices = obj.dimensions(1);
                euler_chr = numberOfVertices + numberOfFaces - numberOfEdges;
                genus = (2-euler_chr)/2;
        end
        
        function faceAreas( obj )
        % this function gets shared vertex data structure (vertices +
        % faces), and calculates the faces areas.
            fLen=obj.dimensions(2);
            areas=zeros(fLen,1);
            for i=1:fLen
                l11=obj.vertices(obj.faces(i,2),:)-obj.vertices(obj.faces(i,1),:);
                l22=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,1),:);
                areas(i)=norm(cross(l11,l22))./2;
            end
            obj.f_areas = areas;
        end
            
        function vertexAreas( obj )
        %   Calculate the area of each vertex as the sum of it's 1-ring faces areas
        %   devided by 3
            %   first, calculate the faces areas:
            if isempty(obj.f_areas)
                obj.faceAreas();
            end
            obj.v_areas = zeros(obj.dimensions(1), 1);
            %for each face f_i add the area of this face/3 to every f_i vertex area 
            for f_i=1:obj.dimensions(2)
                for j=1:3
                    obj.v_areas(obj.faces(f_i,j)) = obj.v_areas(obj.faces(f_i,j)) + obj.f_areas(f_i)/3;
                end
            end
        end
        
        function [ I_F_V, vertex_func ] = face_to_vertex_interp(obj, faces_func)
        %interpolates function of faces to function of vertices
        %   create matrix for the interpolation. I_F_V(i,j)=A_F(j)/( 3*A_V(i) )
        %   while I_F_V is the function from faces to vertices we are cumputing.
        %   A_F(j) is the area of face j iff face j belongs to the one-ring of
        %           vertex i and 0 otherwise.
        %   A_V(i) is the area of vertex i (as defined in the vertex area calc.
        %   arguments f_areas, v_areas are optional. if they are missing they will
        %   be computed.
        %   if -1 is passed in faces_func argument, only I_F_V will be computed.
            rows = zeros(3*obj.dimensions(2),1);
            cols = zeros(3*obj.dimensions(2),1);
            vals = zeros(3*obj.dimensions(2),1);
            for i=1:obj.dimensions(2)
                for j=1:3
                    cols(3*(i-1)+j) = i;
                    rows(3*(i-1)+j) = obj.faces(i,j); 
                    vals(3*(i-1)+j) = obj.f_areas(i)/( 3.*obj.v_areas(obj.faces(i,j)) );
                end
                
            end

            I_F_V = sparse(rows, cols, vals, obj.dimensions(1), obj.dimensions(2));
            if (nargin == 2)
                vertex_func = I_F_V * faces_func;
            end
        end
        
        function [ I_V_F, faces_func ] = vertex_to_face_interp( obj, v_func )
        %	Interpolate vertices function to faces funciton.
        %   I_V_F = A_F^-1 * (I_F_V)^T * A_V
        %   where:
        %       I_V_F is the matrix of desiered interpolation
        %       A_F^-1 is the inverse of the diagonal matrix with areas of faces
        %       I_F_V^T is the transpose of the face to vertex interpolation matrix
        %       A_V is the diagonal matrix with areas of vertices
        %   arguments f_areas, v_areas are optional. if they are missing they will
        %   be computed.
        %   if -1 is passed in faces_func argument, only I_V_F will be computed.
            N_V = obj.dimensions(1);
            N_F = obj.dimensions(2);
            A_F = sparse(1:N_F, 1:N_F, obj.f_areas, N_F, N_F);
            A_V = sparse(1:N_V, 1:N_V, obj.v_areas, N_V, N_V);
            I_F_V = obj.face_to_vertex_interp();
            I_V_F = (A_F\(I_F_V'))*A_V;
            if (nargin == 2)
               faces_func = I_V_F * v_func; 
            end
        end
        
        function showMesh( obj, f, createNewFig )
        %	Load the mesh from the .off file, then display it.
            %   default will compute faces area and display it, if other
            %   function is desired it should be passed to 'f' parameter.
            %   f should be the size of vertices of function of vertices or
            %   number of faces if f is a function of faces.
            %   Also, f should be a column vector.
            %   display:
            if nargin == 1
               f = obj.f_areas; 
            end
            
            if nargin < 3 || (nargin >= 3 && createNewFig == true)
                figure;
            end
            [m,n] = size(f);
            if m == obj.dimensions(1)
                %   display mesh with interpolation for vertices:
                patch('Vertices', obj.vertices, 'Faces', obj.faces,'CData',f,'FaceColor','interp');
                colorbar;
            elseif m == obj.dimensions(2)
                %display mesh with f on faces:
                patch('Vertices', obj.vertices, 'Faces', obj.faces, 'CData', f, 'FaceColor', 'flat');
                colorbar;
            else
                disp('Error, f in of wrong size. should be |V| or |F|.');
            end

        end
        
        function faces_centers = getFacesCenters(obj)
            %   Here we asume we work with TRIANGLES.
            %   This function computes the barycentric center of each face,
            %   and return a vector of length |F|x3 with all centers.
            faces_centers = zeros(obj.dimensions(2),3);
            for i=1:obj.dimensions(2)
                faces_centers(i,:) = obj.vertices(obj.faces(i,1),:)./3 + obj.vertices(obj.faces(i,2),:)./3 + obj.vertices(obj.faces(i,3),:)./3;
            end
        end
        
        function valence = getValence( obj )
            %   Calculate the valence of each vertex and return a vector of
            %   size |V| with all valence values.
            valence = obj.AdjMatrix*ones( obj.dimensions(1), 1 );
        end
        
        function min_angles = getMinAngles( obj )
        % This function calculates for each face the minimum angle it
        % contains. Thus, it returns a vector of size |F|.
            fLen=obj.dimensions(2);
            min_angles=zeros(fLen,1);
            for i=1:fLen
                l1=obj.vertices(obj.faces(i,2),:)-obj.vertices(obj.faces(i,1),:);
                l2=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,1),:);
                l3=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,2),:);
                n1 = norm(l1);
                n2 = norm(l2);
                n3 = norm(l3);
                    
                if(n1 <= n2 && n1 <= n3)
                    %   min angle is at vertex 3
                    min_angles(i) = atan2(norm(cross(l2,l3)),dot(l2,l3));
                elseif (n2 <= n1 && n2 <= n3)
                    %   min angle is at vertex 2
                    min_angles(i) = atan2(norm(cross(-1.*l1,l3)),dot(-1.*l1,l3));
                else
                    %   min angle is at vertex 1
                    min_angles(i) = atan2(norm(cross(l1,l2)),dot(l1,l2));
                end
            end
        end
        
        function N = getFacesNormals(obj)
            %   note that faces in shared vertex data structure determine
            %   the direction of the normal to the face. we will use this
            %   assumption.
            if ~isempty(obj.f_normals)
                N = obj.f_normals;
                return;
            end
            N = zeros(obj.dimensions(2),3);
            for i=1:obj.dimensions(2)
                %   for each face calculate the norm:
                l1=obj.vertices(obj.faces(i,2),:)-obj.vertices(obj.faces(i,1),:);
                l2=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,1),:);
                normal = cross(l1,l2);
                normal = normal./(norm(normal));
                N(i,:) = normal;
            end
            obj.f_normals = N;
        end
        
        function N_V = getVertexNormals(obj)
            %   calculate the vertices normals as defined in the
            %   assignment. n_v(i) = sum_j(A_F(j)*n_F(j))/norm(). so sum
            %   face area times face normal on all 1-ring faces of vertex i.
            N_V = zeros(obj.dimensions(1),3);
            if isempty(obj.f_normals)
                obj.getFacesNormals();
            end
            for i=1:obj.dimensions(2)
               %    for each face add to each vertex the A_F*n_F 
               curFaceData = obj.f_normals(i,:).*obj.f_areas(i);
               for j=1:3
                   N_V(obj.faces(i,j),:) = N_V(obj.faces(i,j),:) + curFaceData;
               end
            end
            %   normalize each vector:
            for i=1:obj.dimensions(1)
                N_V(i,:) = N_V(i,:)./norm(N_V(i,:));
            end
        end
        
        function showVectorField(obj, vf, f, createNewFig) 
            %   This function displays a vector field on top of a mesh (it
            %   draws the mesh first and then draw the vector field on
            %   top). Note that if no vector field is passed, default will
            %   be faces normals.
            if nargin == 1
               vf = obj.getFacesNormals();
            end
            if nargin < 3
               f = ones(obj.dimensions(2),1); 
            end
            if nargin < 4 || (nargin >= 4 && createNewFig == true)
                figure;
            end
            %   show mesh:
            obj.showMesh( f, false);
            %   show vector field on top of it:
            hold on;
            [m,n] = size(vf);
            if m == obj.dimensions(2)
                faces_centers = obj.getFacesCenters();
                quiver3(faces_centers(:,1), faces_centers(:,2), faces_centers(:,3), vf(:,1), vf(:,2), vf(:,3),'-k');
            elseif m == obj.dimensions(1)
                quiver3(obj.vertices(:,1), obj.vertices(:,2), obj.vertices(:,3), vf(:,1), vf(:,2), vf(:,3),'-k');
            end
        end
        
        function grad = getGradient(obj)
            %   This function creates a linear operator (matrix) that
            %   represents the Gradient operator for the mesh.
            %   Given a function f of vertices, grad*f result the gradient
            %   of f on the mesh. note that size of grad is 3*|F| X |V|.
            %   and that the vector v = grad*f will be 3*|F| X 1. 
            %   Note: the above v will be as follows: 
            %   first |F| values will be X coordinates, then next |F|
            %   values will be Y coordinates, and finally, last |F| values
            %   will be Z coordinates. Also, grad is sparse matrix.
            if ~isempty(obj.grad)
                %   if the grad already calculated, return it.
               grad = obj.grad;
               return;
            end
            N_F = obj.dimensions(2);
            row = zeros(9*N_F,1);
            column = zeros(9*N_F,1);
            value = zeros(9*N_F,1);
            if isempty(obj.f_normals)
                obj.getFacesNormals();
            end
            for i=1:N_F
                %   calculate B_a and save it's x,y and z:
                e1=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,2),:);
                grad_B_a = cross(obj.f_normals(i,:), e1);
                %grad_B_a = grad_B_a./(2*obj.f_areas(i));
                for j=1:3
                    row(9*(i-1)+j) = i + (j-1)*N_F;
                    column(9*(i-1)+j) = obj.faces(i,1);
                    value(9*(i-1)+j) = grad_B_a(j);
                end
                
                %   calculate B_b and save it's x,y and z:
                e2=obj.vertices(obj.faces(i,1),:)-obj.vertices(obj.faces(i,3),:);
                grad_B_b = cross(obj.f_normals(i,:), e2);
                %grad_B_b = grad_B_b./(2*obj.f_areas(i));
                for j=1:3
                    row(9*(i-1)+j+3) = i + (j-1)*N_F;
                    column(9*(i-1)+j+3) = obj.faces(i,2);
                    value(9*(i-1)+j+3) = grad_B_b(j);
                end
                
                %   calculate B_c and save it's x,y and z:
                e3=obj.vertices(obj.faces(i,2),:)-obj.vertices(obj.faces(i,1),:);
                grad_B_c = cross(obj.f_normals(i,:), e3);
                %grad_B_c = grad_B_c./(2*obj.f_areas(i));
                for j=1:3
                    row(9*(i-1)+j+6) = i + (j-1)*N_F;
                    column(9*(i-1)+j+6) = obj.faces(i,3);
                    value(9*(i-1)+j+6) = grad_B_c(j);
                end
            end
            
            %   now create the sparse matrix grad:
            E = sparse(row, column, value, 3*N_F, obj.dimensions(1));
            
            G_F_inv = sparse(1:3*N_F, 1:3*N_F, [1./obj.f_areas; 1./obj.f_areas; 1./obj.f_areas], 3*N_F, 3*N_F);
            
            grad = 0.5 .* G_F_inv * E;
            if isempty(obj.grad)
               obj.grad = grad; 
            end
        end
        
        function div = getDivergenceOp(obj)
            %   This function calculates the divergence operator matrix.
            N_V = obj.dimensions(1);
            N_F = obj.dimensions(2);
            G_V = sparse(1:N_V, 1:N_V, obj.v_areas, N_V, N_V);
            G_F = sparse(1:3*N_F, 1:3*N_F, [obj.f_areas; obj.f_areas; obj.f_areas], 3*N_F, 3*N_F);
            if isempty(obj.grad)
                obj.getGradient();
            end
            div = -inv(G_V)*obj.grad'*G_F;
        end
        
        function L = getLaplacian(obj)
            %   calculates the Laplace-Beltmori Operator for the mesh using
            %   L = -div * grad
            div = obj.getDivergenceOp();
            grad_op = obj.getGradient();
            L = -div*grad_op;
        end
        
        function W = getCotWMatrix(obj)
            %   return the cotangent weights matrix using the following
            %   facts: L = 0.25 * G_V^-1 * E^T * G_F^-1 * E, 
            %   W = 0.25 * E^T * G_F^-1 * E => W = G_V * L
            L = obj.getLaplacian();
            N_V = obj.dimensions(1);
            G_V = sparse(1:N_V,1:N_V,obj.v_areas,N_V,N_V);
            W = G_V*L;
        end
        
        function H = getMeanCurvature(obj)
            %   calculate mean curvature using laplacian. H(vi)=0.5norm(L*vi)
            L = obj.getLaplacian();
            H = zeros(obj.dimensions(1),1);
            LV = L*obj.vertices;
            for i=1:obj.dimensions(1)
                H(i) = 0.5 .* norm(LV(i,:));
            end
        end
        
        function K = getGaussianCurvature(obj)
            %   Calculate the gauusian curvature of a mesh
            N_F=obj.dimensions(2);
            N_V=obj.dimensions(1);
            K=zeros(N_V,1);
            for i=1:N_F
                l1=obj.vertices(obj.faces(i,2),:)-obj.vertices(obj.faces(i,1),:);
                l2=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,1),:);
                l3=obj.vertices(obj.faces(i,3),:)-obj.vertices(obj.faces(i,2),:);
                %   angle at verteices:
                angle1 = acos( (l1*l2')/(norm(l1)*norm(l2)) );
                angle2 = acos( ((-1.*l1)*l3')/(norm(l1)*norm(l3)) );
                angle3 = acos( (l3*l2')/(norm(l3)*norm(l2)) );
                %   adding to the sum of angles:
                K(obj.faces(i,1)) = K(obj.faces(i,1)) + angle1;
                K(obj.faces(i,2)) = K(obj.faces(i,2)) + angle2;
                K(obj.faces(i,3)) = K(obj.faces(i,3)) + angle3;
            end
            K = 2*pi*ones(N_V,1) - K;
            K = K ./ (obj.v_areas);
        end
        
        function p = getFaceBarycentricPoint(obj, face_idx, weights)
            %   get a point on the face with index 'face_idx' with
            %   weights = [a, b, c]
            %   barycentric coordinate a*v_1 + b*v_2 + c*v_3
            %   note that a,b,c are >= 0 and sums to 1 (a+b+c=1).
            p = zeros(1,3);
            for i=1:3
                v_i = obj.faces(face_idx,i);
                p = p + weights(i)*obj.vertices(v_i,:);
            end
        end
    end
end

