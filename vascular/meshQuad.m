function [elements, boundaries] = meshQuad(Nr, Ntheta, Ro, theta1, theta2)
Ri = Ro/2;
template = [
  Ri*cos(theta1), Ri*sin(theta1), 0;
  fermatPoint([Ri*cos(theta1), Ri*sin(theta1), 0],[Ri*cos(theta2), Ri*sin(theta2), 0],[Ro*cos(theta2/2+theta1/2), Ro*sin(theta2/2+theta1/2),0]);
  Ri*cos(theta2), Ri*sin(theta2), 0;
  0, 0, 0;
];


boundaries = [];
theta_p = linspace(0,1,Ntheta);
r_p = linspace(0,1,Nr);

bottom = (1-theta_p')*template(3,:) + theta_p'*template(2,:);
top = [
  Ro*cos(-theta_p*(theta2-theta1)/2 + theta2); 
  Ro*sin(-theta_p*(theta2-theta1)/2 + theta2); 
  0*theta_p
]';
[elements1, boundaries1] = meshSideToSide(bottom, top, r_p);
bottom = (1-theta_p')*template(2,:) + theta_p'*template(1,:);
top = [
  Ro*cos(-theta_p*(theta2-theta1)/2 + theta2 - (theta2-theta1)/2); 
  Ro*sin(-theta_p*(theta2-theta1)/2 + theta2 - (theta2-theta1)/2); 
  0*theta_p
]';
[elements2, boundaries2] = meshSideToSide(bottom, top, r_p);


bottom = (1-theta_p')*template(4,:) + theta_p'*template(1,:);
top = (1-theta_p')*template(3,:) + theta_p'*template(2,:);
[elements3, boundaries3] = meshSideToSide(bottom, top, theta_p);

n1 = size(elements1,1);
n2 = size(elements2,1);
boundaries2(:,1) = boundaries2(:,1) + n1;
boundaries3(:,1) = boundaries3(:,1) + n1 + n2;

elements = [elements1;elements2;elements3];
% elements = [elements1;elements2];
boundaries = [boundaries1; boundaries2];
boundaries = boundaries(boundaries(:,3) == 3, :);

boundaries(:, 3) = 1;

end





function P = fermatPoint(A, B, C)
a = norm(B - C);
b = norm(A - C);
c = norm(A - B);

angleA = acos((b^2 + c^2 - a^2)/(2*b*c));
angleB = acos((a^2 + c^2 - b^2)/(2*a*c));
angleC = acos((a^2 + b^2 - c^2)/(2*a*b));

f = @(P) norm(P - A) + norm(P - B) + norm(P - C);
P0 = (A + B + C)/3;

P = fminsearch(f, P0);
end
