function [elements, boundaries, circ] = meshIntersect(Nr, Ntheta, R1, R2, angle)
a = [-1, 0, tan(pi/2-angle)];
a = a / norm(a);
ax = a(1); az = a(3);


theta_p = -linspace(0,1,(Ntheta-1)*8+1)*2*pi + pi/4;
points = zeros((Ntheta-1)*8+1, 3);
for k=1:length(theta_p)
  x = R1*cos(theta_p(k));
  y = R1*sin(theta_p(k));
  
  A = 1 - az^2;
  B = -2*ax*az*x;
  C = x^2 + y^2 - ax^2*x^2 - R2^2;

  disc = B^2 - 4*A*C;
  z = (-B + sqrt(disc)) / (2*A);
  points(k, :) = [x, y, z];
end
elements = [];
boundaries = [];
circ = [];


r_p = linspace(0,1,Nr);
theta_p_loc = linspace(0,1,Ntheta);
for i=1:4
  j = (i-1)*2 + 1;
  k0 = (j-1)*(Ntheta-1) + 1;
  k1 = j*(Ntheta-1) + 1;
  k2 = (j+1)*(Ntheta-1) + 1;

  bottom0 = points(k0, :) / 2;
  bottom2 = points(k2, :) / 2;
  bottom1 = fermatPoint([bottom0; bottom2; points(k1, :)]);

  bottom = (1-theta_p_loc)'*bottom0 + theta_p_loc'*bottom1;
  top    = points(k0:k1, :);
  [elements1, boundaries] = meshSideToSide(bottom, top, r_p);
  elements = [elements;elements1];

  bottom = (1-theta_p_loc)'*bottom1 + theta_p_loc'*bottom2;
  top    = points(k1:k2, :);
  [elements1, boundaries] = meshSideToSide(bottom, top, r_p);
  elements = [elements;elements1];

  bottom = (1-theta_p_loc)'*[0,0,0] + theta_p_loc'*bottom2;
  top = (1-theta_p_loc)'*bottom0 + theta_p_loc'*bottom1;
  [elements1, boundaries] = meshSideToSide(bottom, top, theta_p_loc);
  elements = [elements;elements1];

  circ(i,:,:) = [points(k0:k1, :); points(k1:k2, :)];
end

% bottom = points(Ntheta:2*Ntheta)/2;
% top = points(Ntheta:2*Ntheta)
% [elements2, boundaries2] = meshSideToSide(bottom, top, r_p);
%
% bottom = (1-theta_p')*template(4,:) + theta_p'*template(1,:);
% top = (1-theta_p')*template(3,:) + theta_p'*template(2,:);
% [elements3, boundaries3] = meshSideToSide(bottom, top, theta_p);
%
% n1 = size(elements1,1);
% n2 = size(elements2,1);
% boundaries2(:,1) = boundaries2(:,1) + n1;
% boundaries3(:,1) = boundaries3(:,1) + n1 + n2;
%
% elements = [elements1;elements2;elements3];
% boundaries = [boundaries1; boundaries2];
% boundaries = boundaries(boundaries(:,3) == 3, :);




end





function P = fermatPoint(points)
f = @(P) norm(P - points);
P0 = mean(points);

P = fminsearch(f, P0);
end
