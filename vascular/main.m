addpath('../utils');


R_graft = 3;
R_main = 4;
D_um = 3;
D_depth = 0.5;
angle = 1*pi/4;
Nr = 4;
Ntheta = 4;

Nquad = (Nr-1)*(Ntheta-1)*2 + (Ntheta-1)^2;
Nz_left = 20;
Nz_right = 40;
Nz_graft = 15;

[elements, boundariesCirc] = meshCirc(Nr, Ntheta, R_graft);
top = elements;
top(:, :, 3) = top(:, :, 3) + 20;

elements = [];
boundaries = [];

[intersect, _, intersectCirc] = meshIntersect(Nr, Ntheta, R_graft, R_main, angle);
[elements1, boundaries1] = meshFaceToFace(intersect, top, linspace(0, 1, Nz_graft), boundariesCirc, [1, 1]);
elements1 = rotateY(elements1, angle);
intersect = rotateY(intersect, angle);
intersectCirc = rotateY(intersectCirc, angle);
boundaries1 = boundaries1(boundaries1(:,3) != 5, :);
boundaries1(boundaries1(:,3) == 6, 3) = 2;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];



% Left pipe
[lefto, boundariesCirc] = meshCirc(Nr, Ntheta, R_main);
lefto = rotateY(lefto, pi);
lefto(:, :, 3) = lefto(:, :, 3) + 30;
[left, boundariesCirc] = meshCirc(Nr, Ntheta, R_main);
left = rotateY(left, pi);
for k=1:size(left,1)
  left(:, :, 3) = -left(:, :, 1);
end
left(:, :, 3) = left(:, :, 3) + 15;

[elements1, boundaries1] = meshFaceToFace(lefto,left, linspace(0, 1, Nz_graft), boundariesCirc, [1, 1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
boundaries1(boundaries1(:,3) == 5, 3) = 3;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];




[leftSideM, boundariesLM] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(3,1,:))', -pi+pi/2-pi/8, 0);
[leftSideP, boundariesLP] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(4,1,:))', pi-pi/2+pi/8, 1);
[bypass, boundariesQuad] = meshQuad(Nr, Ntheta, R_main, pi/2+pi/8, 3*pi/2-pi/8);
bypass = rotateX(bypass, pi);
[elements1, boundaries1] = meshFaceToFace(left(1:4*Nquad,:,:),
                                        [bypass;
                                          leftSideM;
                                          intersect(2*Nquad+1:3*Nquad,:,:);
                                          leftSideP], 
                                        linspace(0, 1, Nz_left),
                                        boundariesCirc, [1,1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
boundaries1 = boundaries1(boundaries1(:,3) != 5, :);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];


[right, boundariesCirc] = meshCirc(Nr, Ntheta, R_main);
right(:, :, 3) = right(:, :, 3) - 40;


[rightSideP, boundariesRP] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(1,1,:))', pi-pi/2+pi/8, 0);
[rightSideM, boundariesRM] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(2,1,:))', -pi+pi/2-pi/8, 1);
[bypass, boundariesQuad] = meshQuad(Nr, Ntheta, R_main, pi/2+pi/8, 3*pi/2-pi/8);

[elements1, boundaries1] = meshFaceToFace(
                                right,
                                [intersect(1:Nquad,:,:);
                                 rightSideM;
                                 bypass;
                                 rightSideP;
                                ], linspace(0, 1, Nz_right), boundariesCirc, [1,1]);

boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
boundaries1(boundaries1(:,3) == 5, 3) = 4;
boundaries1(:, 1) += size(elements, 1);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];



[elements1, boundaries1] = meshCap(Nr, Ntheta, R_main, intersect(Nquad+1:2*Nquad,:,:), leftSideM, rightSideM, intersectCirc(2,:,:), -1, 1.0, 0);
boundaries1(:, 1) += size(elements, 1);
boundaries1(boundaries1(:,3) == 6, 3) = 1;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];

[elements1, boundaries1] = meshCap(Nr, Ntheta, R_main, intersect(3*Nquad+1:4*Nquad,:,:), rightSideP, leftSideP, intersectCirc(4,:,:), 1, 1.0, 0);
boundaries1(:, 1) += size(elements, 1);
boundaries1(boundaries1(:,3) == 6, 3) = 1;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];


% Project
for k=size(intersect,1)*Nz_graft+1:size(elements, 1)
  for m=1:8
    p = squeeze(elements(k,m,:))';
    R_shrink = sqrt(p(1)^2+p(2)^2);
    if (R_shrink == 0)
      continue
    end
    circ_shrink = intersectCirc*R_shrink/R_main;
    [d, pos1, pos2] = distToCirc(p, circ_shrink, R_shrink, 4);
    D_um_loc = D_um + 4*norm(squeeze(intersectCirc(pos1, pos2,:)-circ_shrink(pos1,pos2,:)));
    % D_um_loc = D_um*(R_main/R_shrink)^2;
    R_local = 1 - D_depth/R_main*sin(min(1,max(d/D_um_loc,0))*pi);

    % R_local = R_local*(p(1)^2+p(2)^2)/R_main^2;
    elements(k,m,1:2) = R_local*elements(k,m,1:2);
  end
end;

checkLeftHanded(elements);


N = size(elements,1);
X = permute(elements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk('outer2.vtk', X,Hexes, boundaries,'');

