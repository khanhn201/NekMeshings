function [elements, boundaries] = meshSide(Nr, Ntheta, R_main, point, angle, rev, R_rat, mult_r)
theta_p = linspace(0,1,Ntheta);

pleft = point;
pbottom = [R_main*cos(angle), R_main*sin(angle), 0];
theta_center = (atan2(pleft(2),pleft(1))+angle)/2;
pleftc = [R_main*cos(theta_center), R_main*sin(theta_center), pleft(3)/2];

points = zeros((Ntheta-1)*2+1, 3);
theta = linspace(angle, atan2(pleft(2),pleft(1)), Ntheta*2-1);
z_arr = linspace(pbottom(3), pleft(3), Ntheta*2-1);

if rev == 1 
  ptemp = pleft;
  pleft = pbottom;
  pbottom = ptemp;
  theta = flip(theta);
  z_arr = flip(z_arr);
end

for k=1:length(theta)
  x = R_main*cos(theta(k));
  y = R_main*sin(theta(k));
  points(k, :) = [x, y, z_arr(k)];
end



bottom2 = pleft*R_rat;
bottom0 = pbottom*R_rat;
bottom1 = fermatPoint([bottom0; bottom2; pleftc]);

r_p = 1:Nr;
r_p = 1-(mult_r.^(Nr-(r_p-1))-mult_r)/(mult_r^(Nr)-mult_r);

bottom = (1-theta_p)'*bottom0 + theta_p'*bottom1;
top    = points(1:Ntheta, :);
[elements1, boundaries1] = meshSideToSide(bottom, top, r_p);
boundaries1 = boundaries1(boundaries1(:,3) == 3, :);
boundaries1(:, 3) = 1;

bottom = (1-theta_p)'*bottom1 + theta_p'*bottom2;
top    = points(Ntheta:end, :);
[elements2, boundaries2] = meshSideToSide(bottom, top, r_p);
boundaries2 = boundaries2(boundaries2(:,3) == 3, :);
boundaries2(:, 3) = 1;
boundaries2(:, 1) += size(elements1, 1);


bottom = (1-theta_p)'*[0,0,0] + theta_p'*bottom2;
top = (1-theta_p)'*bottom0 + theta_p'*bottom1;
[elements3, boundaries3] = meshSideToSide(bottom, top, theta_p);

elements = [elements1;elements2; elements3];
boundaries = [boundaries1; boundaries2];


end





function P = fermatPoint(points)
f = @(P) norm(P - points);
P0 = mean(points);

P = fminsearch(f, P0);
end
