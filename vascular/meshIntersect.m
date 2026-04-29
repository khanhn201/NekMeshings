function [elements, boundaries, circ] = meshIntersect(Nr, Ntheta, R1, R2, angle, R_rat, mult_r)
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


r_p = 1:Nr;
r_p = 1-(mult_r.^(Nr-(r_p-1))-mult_r)/(mult_r^(Nr)-mult_r);
theta_p_loc = linspace(0,1,Ntheta);
for i=1:4
  j = (i-1)*2 + 1;
  k0 = (j-1)*(Ntheta-1) + 1;
  k1 = j*(Ntheta-1) + 1;
  k2 = (j+1)*(Ntheta-1) + 1;

  bottom0 = points(k0, :)*R_rat;
  bottom2 = points(k2, :)*R_rat;
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

end





function P = fermatPoint(points)
f = @(P) norm(P - points);
P0 = mean(points);

P = fminsearch(f, P0);
end
