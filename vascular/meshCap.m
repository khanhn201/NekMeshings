function [elements, boundaries] = meshCap(Nr, Ntheta, R_main, top_face, left_face, right_face, top_circ, dir, D_um, D_depth)
% plotElements2D(left_face, [])
Nring = (Nr-1)*(Ntheta-1);


pface = fermatPoint([
  squeeze(top_face(Nring, 3, :))';
  squeeze(right_face(Nring, 3, :))';
  squeeze(left_face(Nring, 3, :))';
]);
x = pface(1);
y = dir*sqrt(R_main^2-x^2);
z = pface(3);
d = distToCirc(x,y,z, top_circ, R_main);
R_local = 1 - D_depth/R_main*sin(min(1,max(d/D_um,0))*pi);
x = R_local*x;
y = R_local*y;

pface(1) = x;
pface(2) = y;
pface(3) = z;
pcenter = fermatPoint([
  pface;
  squeeze(right_face(Nring+1, 1, :))';
  squeeze(left_face(Nring+1, 1, :))';
  squeeze(top_face(Nring+1, 3, :))';
]);

theta_p = linspace(0,1,Ntheta);
r_p = linspace(0,1,Nr);

bottom = (1-theta_p)'*right_face(1, 1, :) + theta_p'*top_face(Nring+1, 1, :);
top = (1-theta_p)'*right_face(Ntheta-1, 2,:) + theta_p'*pcenter;
[center_face, boundariesf] = meshSideToSide(bottom, top, theta_p);
[elements, boundaries] = meshFaceToFace(left_face(2*Nring+1:end,:,:), center_face, theta_p, boundariesf, []);




theta1 = atan2(right_face(Nring, 3,2),right_face(Nring, 3,1));
theta2 = atan2(pface(2),pface(1));
theta = linspace(theta1, theta2, Ntheta);
right_cir = zeros(Ntheta, 3);
right_cir(:, 1) = R_main*cos(theta);
right_cir(:, 2) = R_main*sin(theta);
right_cir(:, 3) = (1-theta_p)'*right_face(Nring, 3,3) + theta_p'*pface(3);

surface = zeros(Ntheta, Ntheta, 3);
for i = 1:Ntheta
    theta1 = atan2(right_cir(i,2),right_cir(i,1));
    theta2 = atan2(top_circ(1,i,2),top_circ(1,i,1));
    theta = linspace(theta1, theta2, Ntheta);
    for j = 1:Ntheta
        rat = (j-1)/(length(theta)-1);
        z = (1-rat)*right_cir(i,3) + rat*top_circ(1,i,3);
        x = R_main*cos(theta(j));
        y = R_main*sin(theta(j));
        d = distToCirc(x,y,z, top_circ, R_main);
        R_local = R_main - D_depth*sin(min(1,max(d/D_um,0))*pi);
        x = R_local*cos(theta(j));
        y = R_local*sin(theta(j));
        surface(i, j, :) = [x, y, z];
    end
end
surf_elements = [];
for i = 1:Ntheta-1
    for j = 1:Ntheta-1
        surf_elements(end+1,1,:) = squeeze(surface(i,   j, :))';
        surf_elements(end,2,:) = squeeze(surface(i,   j+1,   :))';
        surf_elements(end,3,:) = squeeze(surface(i+1, j+1,   :))';
        surf_elements(end,4,:) = squeeze(surface(i+1, j, :))';
    end
end
bottom = (1-theta_p)'*right_face(Ntheta-1, 2,:) + theta_p'*right_face(1, 1, :);
top = (1-theta_p)'*pcenter + theta_p'*top_face(Nring+1, 1, :);
[center_face, boundariesf] = meshSideToSide(bottom, top, theta_p);
[elements1, boundaries1] = meshFaceToFace(center_face, surf_elements, r_p, boundariesf, []);
boundaries1(:, 1) += size(elements, 1);
elements = [elements; elements1];
boundaries = [boundaries1];


theta1 = atan2(pface(2),pface(1));
theta2 = atan2(left_face(Nring, 3,2),left_face(Nring, 3,1));
theta = linspace(theta1, theta2, Ntheta);
left_cir = zeros(Ntheta, 3);
left_cir(:, 1) = R_main*cos(theta);
left_cir(:, 2) = R_main*sin(theta);
left_cir(:, 3) = (1-theta_p)'*pface(3) + theta_p'*left_face(Nring, 3,3);

surface = zeros(Ntheta, Ntheta, 3);
for i = 1:Ntheta
    theta1 = atan2(left_cir(i,2),left_cir(i,1));
    theta2 = atan2(top_circ(1,i+Ntheta,2),top_circ(1,i+Ntheta,1));
    theta = linspace(theta1, theta2, Ntheta);
    for j = 1:Ntheta
        rat = (j-1)/(length(theta)-1);
        z = (1-rat)*left_cir(i,3) + rat*top_circ(1,i+Ntheta,3);
        x = R_main*cos(theta(j));
        y = R_main*sin(theta(j));
        d = distToCirc(x,y,z, top_circ, R_main);
        R_local = R_main - D_depth*sin(min(1,max(d/D_um,0))*pi);
        x = R_local*cos(theta(j));
        y = R_local*sin(theta(j));

        surface(i, j, :) = [x, y, z];
    end
end
surf_elements = [];
for i = 1:Ntheta-1
    for j = 1:Ntheta-1
        surf_elements(end+1,1,:) = squeeze(surface(i,   j, :))';
        surf_elements(end,2,:) = squeeze(surface(i,   j+1,   :))';
        surf_elements(end,3,:) = squeeze(surface(i+1, j+1,   :))';
        surf_elements(end,4,:) = squeeze(surface(i+1, j, :))';
    end
end
bottom = (1-theta_p)'*pcenter + theta_p'*top_face(Nring+1, 1, :);
top = (1-theta_p)'*left_face(Nring+1, 1, :) + theta_p'*left_face(Nring+Ntheta-1, 2,:);
[center_face, boundariesf] = meshSideToSide(bottom, top, theta_p);
[elements1, boundaries1] = meshFaceToFace(center_face, surf_elements, r_p, boundariesf, []);
boundaries1(:, 1) += size(elements, 1);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];


theta1 = atan2(left_face(1, 1,2),left_face(1, 1,1));
theta2 = atan2(left_face(Nring, 3,2),left_face(Nring, 3,1));
theta = linspace(theta1, theta2, Ntheta);
botleft_cir = zeros(Ntheta, 3);
botleft_cir(:, 1) = R_main*cos(theta);
botleft_cir(:, 2) = R_main*sin(theta);
botleft_cir(:, 3) = theta_p'*left_face(Nring, 3,3);


surface = zeros(Ntheta, Ntheta, 3);
for i = 1:Ntheta
    theta1 = atan2(botleft_cir(i,2),botleft_cir(i,1));
    theta2 = atan2(right_cir(i,2),right_cir(i,1));
    theta = linspace(theta1, theta2, Ntheta);
    for j = 1:Ntheta
        rat = (j-1)/(length(theta)-1);
        z = (1-rat)*botleft_cir(i,3) + rat*right_cir(i,3);
        x = R_main*cos(theta(j));
        y = R_main*sin(theta(j));
        d = distToCirc(x,y,z, top_circ, R_main);
        R_local = R_main - D_depth*sin(min(1,max(d/D_um,0))*pi);
        x = R_local*cos(theta(j));
        y = R_local*sin(theta(j));
        surface(i, j, :) = [x, y, z];
    end
end
surf_elements = [];
for i = 1:Ntheta-1
    for j = 1:Ntheta-1
        surf_elements(end+1,1,:) = squeeze(surface(i,   j, :))';
        surf_elements(end,2,:) = squeeze(surface(i,   j+1,   :))';
        surf_elements(end,3,:) = squeeze(surface(i+1, j+1,   :))';
        surf_elements(end,4,:) = squeeze(surface(i+1, j, :))';
    end
end
bottom = (1-theta_p)'*right_face(Nring+Ntheta-1, 2, :) + theta_p'*right_face(Nring+1, 1, :);
top = (1-theta_p)'*left_face(Nring+1, 1, :) + theta_p'*pcenter;
[center_face, boundariesf] = meshSideToSide(bottom, top, theta_p);
[elements1, boundaries1] = meshFaceToFace(center_face, surf_elements, r_p, boundariesf, []);
boundaries1(:, 1) += size(elements, 1);
elements = [elements; elements1];
boundaries = [boundaries; boundaries1];



boundaries = boundaries(boundaries(:,3) == 6, :);


end

function P = fermatPoint(points)
f = @(P) norm(P - points);
P0 = mean(points);

P = fminsearch(f, P0);
end

function d = distToCirc(x,y,z, circ, R)
theta = atan2(y,x)+2*pi;
s = R * theta;
Nc = size(circ,2);
d = inf;

for k = 1:Nc-1
    a3 = squeeze(circ(1,k,:))';
    b3 = squeeze(circ(1,k+1,:))';

    ta = atan2(a3(2), a3(1))+2*pi;
    tb = atan2(b3(2), b3(1))+2*pi;

    sa = R * ta;
    sb = R * tb;

    a = [sa, a3(3)];
    b = [sb, b3(3)];
    p = [s, z];
    ab = b - a;
    t = dot(p-a, ab) / dot(ab, ab);
    t = max(0, min(1, t));
    proj = a + t*ab;

    d = min(d, norm(p - proj));
end
end
