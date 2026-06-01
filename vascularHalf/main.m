addpath('../utils');


R_graft = 3;
R_main = 4;
D_um = 3;
D_depth = 1.5;
angle = 1*pi/4;

Nr = 4;
Ntheta = 4;

R_rat = 0.6;
mult_r = 1.5;

Nz_leftfaro = 6;
Nz_leftfar = 30;
Nz_leftclose = 10;
Nz_rightfar = 35;
Nz_rightclose = 10;
Nz_graftclose = 30;
Nz_graftfar = 6;


center = [0.0,0,0.0];


Nquad = (Nr-1)*(Ntheta-1)*2 + (Ntheta-1)^2;
Nfan = (Nr-1)*(Ntheta-1)*2;

[grafto, boundariesCirc] = meshCirc(Nr, Ntheta, R_graft, R_rat, mult_r);
grafto(:, :, 3) = grafto(:, :, 3) + 80;
[graft, boundariesCirc] = meshCirc(Nr, Ntheta, R_graft, R_rat, mult_r);
graft(:, :, 3) = graft(:, :, 3) + 70;

elements = [];
boundaries = [];

[intersect, _, intersectCirc] = meshIntersect(Nr, Ntheta, R_graft, R_main, angle, R_rat, mult_r, center);
[elements1, boundaries1] = meshFaceToFace(intersect, graft, linspace(0, 1, Nz_graftclose), boundariesCirc, [1, 1]);
boundaries1 = boundaries1(boundaries1(:,3) != 5, :);
boundaries1(boundaries1(:,3) == 6, 3) = 5;
boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];

[elements1, boundaries1] = meshFaceToFace(graft, grafto, linspace(0, 1, Nz_graftfar), boundariesCirc, [1, 1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) != 5, :);
boundaries1(boundaries1(:,3) == 6, 3) = 2;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];

elements = rotateY(elements, angle);
intersect = rotateY(intersect, angle);
intersectCirc = rotateY(intersectCirc, angle);


Ry = [ cos(angle), 0, sin(angle);
       0,          1, 0;
      -sin(angle), 0, cos(angle)];
center = (Ry * center')';

% Left pipe
[lefto, boundariesCirc] = meshCirc(Nr, Ntheta, R_main, R_rat, mult_r);
lefto = rotateY(lefto, pi);
lefto(:, :, 3) = lefto(:, :, 3) + 80;
[leftc, boundariesCirc] = meshCirc(Nr, Ntheta, R_main, R_rat, mult_r);
leftc = rotateY(leftc, pi);
leftc(:, :, 3) = leftc(:, :, 3) + 70;

[left, boundariesCirc] = meshCirc(Nr, Ntheta, R_main, R_rat, mult_r);
left = rotateY(left, pi);
left(:, :, 3) = -1.5*max(0,left(:, :, 1));
left(:, :, 3) = left(:, :, 3) + 18;

[elements1, boundaries1] = meshFaceToFace(lefto,leftc, linspace(0, 1, Nz_leftfaro), boundariesCirc, [1, 1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
boundaries1(boundaries1(:,3) == 5, 3) = 3;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];

[elements1, boundaries1] = meshFaceToFace(leftc,left, linspace(0, 1, Nz_leftfar), boundariesCirc, [1, 1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
boundaries1(boundaries1(:,3) == 5, 3) = 6;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];



angle1 = atan2(intersect(3*Nquad+Nfan,3,2), intersect(3*Nquad+Nfan,3,1));

[bypass1, boundariesLM] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(1,1,:))', pi, 1, R_rat, mult_r, center);
[bypass2, boundariesLM] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(3,1,:))', -pi, 0, R_rat, mult_r, center);
[elements1, boundaries1] = meshFaceToFace(left,
                                        [bypass1; bypass2;
                                          intersect(2*Nquad+1:4*Nquad,:,:)], 
                                        linspace(0, 1, Nz_leftclose),
                                        boundariesCirc, [1,1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) == 1, :);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];


% Right pipe
[righto, boundariesCirc] = meshCirc(Nr, Ntheta, R_main, R_rat, mult_r);
righto(:, :, 3) = righto(:, :, 3) - 80;
[right, boundariesCirc] = meshCirc(Nr, Ntheta, R_main, R_rat, mult_r);
% right(:, :, 3) = 1.0*right(:, :, 1);
right(:, :, 3) = right(:, :, 3) - 6;

[elements1, boundaries1] = meshFaceToFace(righto,right, linspace(0, 1, Nz_rightfar), boundariesCirc, [1, 1]);
boundaries1(:, 1) += size(elements, 1);
boundaries1 = boundaries1(boundaries1(:,3) != 6, :);
boundaries1(boundaries1(:,3) == 5, 3) = 4;
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];

[bypass1, boundariesLM] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(3,1,:))', -pi, 1, R_rat, mult_r, center);
[bypass2, boundariesLM] = meshSide(Nr, Ntheta, R_main, squeeze(intersectCirc(1,1,:))', pi, 0, R_rat, mult_r, center);
[elements1, boundaries1] = meshFaceToFace(
                                right,
                                [intersect(0*Nquad+1:2*Nquad,:,:);
                                  bypass1; bypass2 
                                ], linspace(0, 1, Nz_rightclose), boundariesCirc, [1,1]);

boundaries1 = boundaries1(boundaries1(:,3) == 1, :);
boundaries1(:, 1) += size(elements, 1);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];



% Project
center_cast = reshape(center,1,1,3);
for k=size(intersect,1)*(Nz_graftclose+Nz_graftfar)+1:size(elements, 1)
  for m=1:8
    p = squeeze(elements(k,m,:))';
    R_shrink = sqrt(p(1)^2+p(2)^2);
    if (R_shrink == 0)
      continue
    end
    circ_shrink = (intersectCirc-center_cast)*R_shrink/R_main + center_cast;
    [d, pos1, pos2] = distToCirc(p, circ_shrink, R_shrink, 4);
    v1 = squeeze(intersectCirc(pos1, pos2,:)-circ_shrink(pos1,pos2,:));
    v2 = p'-squeeze(circ_shrink(pos1,pos2,:));
    v1n = v1(:);
    v2n = v2(:);
    proj = (v1n' * v2n) / (v2n' * v2n) * v2n;

    D_um_loc = D_um + 2*norm(proj);
    % D_um_loc = D_um*(R_main/R_shrink)^2;
    R_local = 1 - D_depth/R_main*sin(min(1,max(d/D_um_loc,0))*pi)*...
                  max(cos(atan2(p(2), p(1))),0);

    % R_local = R_local*(p(1)^2+p(2)^2)/R_main^2;
    elements(k,m,1:2) = R_local*elements(k,m,1:2);
  end
end;

checkLeftHanded(elements);


size(elements)
size(boundaries)
N = size(elements,1);
X = permute(elements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk('outer2.vtk', X,Hexes, boundaries,'');

exportRE2("ab", elements, boundaries)
