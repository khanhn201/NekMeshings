function [elements, boundaries] = meshCirc(Nr, Ntheta, Ro)
Ri = Ro/2;
template = [
  Ri/sqrt(2), -Ri/sqrt(2), 0;
  fermatPoint([Ri/sqrt(2), -Ri/sqrt(2), 0],[Ri/sqrt(2), Ri/sqrt(2), 0],[Ro, 0,0]);
  Ri/sqrt(2), Ri/sqrt(2), 0;
  0, 0, 0;
];


boundaries = [];
theta_p = linspace(0,1,Ntheta);
r_p = linspace(0,1,Nr);

bottom = (1-theta_p')*template(3,:) + theta_p'*template(2,:);
top = [Ro*cos(-theta_p*pi/4+pi/4); Ro*sin(-theta_p*pi/4+pi/4); 0*theta_p]';
[elements1, boundaries1] = meshSideToSide(bottom, top, r_p);
bottom = (1-theta_p')*template(2,:) + theta_p'*template(1,:);
top = [Ro*cos(-theta_p*pi/4); Ro*sin(-theta_p*pi/4); 0*theta_p]';
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


elements_all = [];
boundaries_all = [];

Ne = size(elements,1);
for k = 0:3
    theta = -k * pi/2;
    R = [cos(theta), -sin(theta), 0;
         sin(theta),  cos(theta), 0;
         0,           0,          1];

    elements_rot = elements;

    for e = 1:size(elements,1)
        for n = 1:4
            p = squeeze(elements(e,n,:))';
            elements_rot(e,n,:) = (R * p')';
        end
    end

    elements_all = [elements_all; elements_rot];
    boundaries_rot = boundaries;
    boundaries_rot(:,1) = boundaries_rot(:,1) + k*Ne;
    boundaries_all = [boundaries_all; boundaries_rot];
end

elements = elements_all;
boundaries = boundaries_all;


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
