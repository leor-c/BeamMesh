classdef Transformations3D
    %   Do 3d transformations 
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function translatedShape = translate3DShape(shape, target)
            %   move the given shape (matrix of 3d point coordinates) to
            %   the given target location. the initial location is assumed
            %   to be the mean of all shape's points.
            %   matrix is supposed to contain each point in a column, and
            %   if it is not in homogenus coordinates, this function will
            %   handle it.
            %   * target should be column vector also.
            n = size(shape,1);
            if (n == 3)
                %   transform to homogenus coordinates:
                shape = Transformations3D.getHomogenusCoordinates(shape);
            end
            translationMat = Transformations3D.getTranslationMat(shape, target);
            
            translatedShape = translationMat * shape;
            
            if (n == 3)
                translatedShape = translatedShape(1:3,:);
            end
        end
        
        function translationMat = getTranslationMat(shape, target)
            %   Same as translate3DShape above, just returns the matrix
            %   without performing operations on the shape.
            center = mean(shape,2);
            
            translationMat = speye(4);
            translationMat(1:3,4) = target(1:3) - center(1:3);
        end
        
        function homogenusA = getHomogenusCoordinates(A)
            %    transforms euclidean coordinates of set of points to
            %    homogenus ones. assumes all points in A are columns of
            %    A.
            homogenusA = [A; ones(1,size(A,2))];
        end
        
        function newShape = operateOnShape(opMatrix, shape)
            %   this funciton is doing the operations on the shape, and if
            %   necessarry converts to homogenus coords and back
            n = size(shape,1);
            if (n == 3)
                %   transform to homogenus coordinates:
                shape = Transformations3D.getHomogenusCoordinates(shape);
            end
            newShape = opMatrix * shape;
            if(n == 3)
                newShape = newShape(1:3,:);
            end
        end
        
        function rotatedShape = rotateAroundXAxis(shape, theta)
            %   this is the raw operation of rotation around X axis in 3d.
            rotationMat = speye(4);
            cosTheta = cos(theta);
            sinTheta = sin(theta);
            
            rotationMat(2:3,2:3) = [cosTheta -sinTheta;...
                                    sinTheta cosTheta];
            rotatedShape = Transformations3D.operateOnShape(rotationMat, shape);
        end
        
        function rotatedShape = rotateAroundYAxis(shape, theta)
            %   this is the raw operation of rotation around Y axis in 3d.
            rotationMat = speye(4);
            cosTheta = cos(theta);
            sinTheta = sin(theta);
            
            rotationMat([1,3],[1,3]) = [cosTheta -sinTheta;...
                                    sinTheta cosTheta];
            rotatedShape = Transformations3D.operateOnShape(rotationMat, shape);
        end
        
        function rotatedShape = rotateAroundZAxis(shape, theta)
            %   this is the raw operation of rotation around Z axis in 3d.
            rotationMat = Transformations3D.rotateAroundZAxisMat(theta);
            
            rotatedShape = Transformations3D.operateOnShape(rotationMat, shape);
        end
        
        function rotationMat = rotateAroundZAxisMat(theta)
            %   this is the raw operation of rotation around Z axis in 3d.
            rotationMat = speye(4);
            cosTheta = cos(theta);
            sinTheta = sin(theta);
            
            rotationMat(1:2,1:2) = [cosTheta -sinTheta;...
                                    sinTheta cosTheta];
        end
        
        function rotatedShape = rotateAroundArbitraryAxis(shape, axis, theta)
            %   this is the raw operation of rotation around arbitrary axis
            %   (given as input) in 3d.
            %   according to page 203 in book.
            %   normalize axis vector:
            axis = axis ./ norm(axis);
            
            %   align rotation axis with z axis:
            rotateToZ = Transformations3D.alignWithZAxisMat(axis);
            rotateFromZ = Transformations3D.alignWithZAxisMat(-axis);
            rotateAroundZ = Transformations3D.rotateAroundZAxisMat(theta);
            
            translationMat = Transformations3D.getTranslationMat(shape, [0;0;0]);
            translationBackMat = speye(4);
            translationBackMat(1:3,4) = -translationMat(1:3,4);
            
            rotationMat = translationBackMat * rotateFromZ * rotateAroundZ * ...
                rotateToZ * translationMat;
            
            rotatedShape = Transformations3D.operateOnShape(rotationMat, shape);
        end
        
        function rotatedShape = alignWithZAxis(shape, axis)
            %   align shape axis with z axis:
            rotationMat = Transformations3D.alignWithZAxisMat(axis);
            rotatedShape = Transformations3D.operateOnShape(rotationMat, shape);
        end
        
        function rotationMat = alignWithZAxisMat(axis)
            %   assuming axis is a given column vector we want to align
            %   with Z axis:
            %   projection of shape axis on YZ plane:
            yz_projection = norm(axis(2:3));
            if(yz_projection == 0)
                %   axis is X axis. rotate to Z axis..
                cosTheta_X = 0;
                sinTheta_X = 1; 
                R_x = speye(4);
                R_x(2:3,2:3) = [cosTheta_X -sinTheta_X;...
                                sinTheta_X cosTheta_X];
                rotationMat = R_x;
            else
                %rotate around X:
                cosTheta_X = axis(3) / yz_projection;
                sinTheta_X = axis(2) / yz_projection; 
                R_x = speye(4);
                R_x(2:3,2:3) = [cosTheta_X -sinTheta_X;...
                                sinTheta_X cosTheta_X];
                %   rotate around Y:            
                cosTheta_Y = yz_projection;
                sinTheta_Y = axis(1); 
                R_y = speye(4);
                R_y([1,3],[1,3]) = [cosTheta_Y -sinTheta_Y;...
                                    sinTheta_Y cosTheta_Y];
                rotationMat = R_y * R_x;
            end
        end
    end
    
end

