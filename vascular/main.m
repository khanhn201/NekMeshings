addpath('../utils');


R_graft = 3;
R_main = 4;
angle = pi/4;
Nr = 7;
Ntheta = 3;

Nquad = (Nr-1)*(Ntheta-1)*2 + (Ntheta-1)^2;

[elements, boundaries] = meshCirc(Nr, Ntheta, R_graft);
top = elements;
top(:, :, 3) = top(:, :, 3) + 30;

% plotElements2D(elements, boundaries);

% z_levels = linspace(0, 8, 20);
% [elements, boundaries] = extrudeZ(elements, z_levels);
% plotElements3D(elements);

[intersect, boundaries] = meshIntersect(Nr, Ntheta, R_graft, R_main, angle);

p = linspace(0, 1, 10);
[elements, boundaries] = meshFaceToFace(intersect, top, p);
% [elements, boundaries] = meshFaceToFace(top,bottom, p);



[left, boundaries] = meshCirc(Nr, Ntheta, R_main);
left = rotateY(left, pi);
left(:, :, 3) = left(:, :, 3) + 30;
left = rotateY(left, angle);
[elements1, boundaries] = meshFaceToFace(left(2*Nquad+1:3*Nquad,:,:),intersect(2*Nquad+1:3*Nquad,:,:), p);
elements = [elements; elements1];

[right, boundaries] = meshCirc(Nr, Ntheta, R_main);
right = rotateY(right, pi);
right(:, :, 3) = right(:, :, 3) + 40;
right = rotateY(right, -pi+angle);
[elements1, boundaries] = meshFaceToFace(right(1:Nquad,:,:),intersect(1:Nquad,:,:), p);
elements = [elements; elements1];

[bypass, boundaries] = meshCirc(Nr, Ntheta, R_main);
bypass = rotateY(bypass, angle);
[elements1, boundaries] = meshFaceToFace(right(2*Nquad+1:3*Nquad,:,:),bypass(2*Nquad+1:3*Nquad,:,:), p);
elements = [elements; elements1];

bypass = rotateY(bypass, pi);
[elements1, boundaries] = meshFaceToFace(left(1:Nquad,:,:),bypass(1:Nquad,:,:), p);
elements = [elements; elements1];


checkLeftHanded(elements)


N = size(elements,1);
X = permute(elements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk('outer2.vtk', X,Hexes, boundaries,'');

