% function d = distToCirc(p, circ, R, m=1)
% theta = atan2(p(2),p(1));
% s = R * theta;
% Nc = size(circ,2);
% d = inf;
%
% for j = 1:m
% for k = 1:Nc-1
%     a3 = squeeze(circ(j,k,:))';
%     b3 = squeeze(circ(j,k+1,:))';
%
%     ta = atan2(a3(2), a3(1));
%     tb = atan2(b3(2), b3(1));
%
%     sa = R * ta;
%     sb = R * tb;
%
%     a = [sa, a3(3)];
%     b = [sb, b3(3)];
%     pp = [s, p(3)];
%     ab = b - a;
%     t = dot(pp-a, ab) / dot(ab, ab);
%     t = max(0, min(1, t));
%     proj = a + t*ab;
%
%     d = min(d, norm(pp - proj));
% end
% end
% end

function [d, pos1, pos2] = distToCirc(p, circ, R, m=1)
theta = atan2(p(2),p(1));
s = R * theta;
Nc = size(circ,1)*size(circ,2);
d = inf;
s_circ = zeros(Nc,1);
z_circ = zeros(Nc,1);

for k = 1:size(circ,1)
for j = 1:size(circ,2)
    pt = squeeze(circ(k,j,:))';
    t = atan2(pt(2), pt(1));
    s_circ(k*size(circ,2) + j) = R * t;
    z_circ(k*size(circ,2) + j) = pt(3);
end
end

s_circ = unwrap(s_circ / R) * R;
pp = [R*atan2(p(2),p(1)), p(3)];

[in, on] = inpolygon(pp(1), pp(2), s_circ, z_circ);

pos1 = 1;
pos2 = 1;
if in || on
    d = 0;
    return;
end

for j = 1:m
for k = 1:size(circ,2)-1
    a3 = squeeze(circ(j,k,:))';
    b3 = squeeze(circ(j,k+1,:))';

    % unwrap endpoints
    ta = atan2(a3(2), a3(1));
    tb = atan2(b3(2), b3(1));

    sa = R * ta;
    sb = R * tb;

    a = [sa, a3(3)];
    b = [sb, b3(3)];
    pp = [s, p(3)];

    % --- handle periodic wrap ---
    for shift = [-2*pi*R, 0, 2*pi*R]
        a_shift = a; a_shift(1) = a_shift(1) + shift;
        b_shift = b; b_shift(1) = b_shift(1) + shift;

        ab = b_shift - a_shift;
        t = dot(pp - a_shift, ab) / dot(ab, ab);
        t = max(0, min(1, t));

        proj = a_shift + t*ab;

        d = min(d, norm(pp - proj));
        pos1 = j;
        pos2 = k;
    end
end
end
end
