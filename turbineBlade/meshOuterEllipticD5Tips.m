function [elements,boundaries, pp_coarse] = meshOuterEllipticD5Tips(pp, arc_length, arc_length_at_max_y, flipped)
    config
    nx=2*n_top + 4*k_inner + 4;
    ny=k_outer;
    maxit=100000;
    Ermax=10^-9;

    point_at_min = ppval(pp, 0);
    x = point_at_min(1);
    i = 0:nx/2;
    % i = 0:(nx/2);
    cheb_nodes = cos(i*pi / (nx/2));
    cheb_nodes_mapped = (cheb_nodes + 1) / 2;
    s_fine = cheb_nodes_mapped * arc_length_at_max_y;
    s_fine1 = flip(s_fine)(1:end-1);
    s_fine = cheb_nodes_mapped * (arc_length-arc_length_at_max_y) + arc_length_at_max_y;
    s_fine2 = flip(s_fine)(1:end-1);
    s_fine = [s_fine1(1), s_fine1(3:end-1),...
              s_fine2(1), s_fine2(3:end-1)];
    % s_fine = [s_fine1, s_fine2];

    s_fine = [linspace(0, arc_length_at_max_y, n_top+2*k_inner + 1)(1:end-1), ...
          linspace(arc_length_at_max_y, arc_length, n_bottom+2*k_inner + 1)(1:end-1)];
    if flipped == true 
        s_fine(s_fine >= arc_length_at_max_y) = s_fine(s_fine >= arc_length_at_max_y) - arc_length;
    end
    all_points = ppval(pp, s_fine)';
    if flipped == true
        s_fine1 = circshift(s_fine, -n_top - 2*k_inner);
        all_points1 = circshift(all_points, -n_top - 2*k_inner);
        s_fine1 = [s_fine1'; arc_length_at_max_y];
        all_points1 = [all_points1; all_points1(1, :)];
    else
        s_fine1 = [s_fine'; arc_length];
        all_points1 = [all_points; all_points(1,:)];
    end
    all_points1 = [
        (all_points1(2,:)-all_points1(1,:))/(s_fine1(2)-s_fine1(1)); 
        all_points1; 
        (all_points1(end,:)-all_points1(end-1,:))/(s_fine1(end)-s_fine1(end-1));];
    pp_coarse = spline(s_fine1, all_points1');


    layer_next = [];
    for t = 1:length(all_points)/2
        p1 = all_points(t, :);
        if t == 1
            p2 = all_points(end, :);
        else
            p2 = all_points(t-1, :);
        end
        p3 = all_points(t+1, :);

        [q1, q2, q3] = findQuarterBisectors(p1, p2, p3);
        [q4, q5] = findTrisects(p1, p3, p2);
        if t == 1
            U = p1 + first_layer_thickness*q4;
            layer_next = [layer_next; U];
            U = p1 + first_layer_thickness*q2;
            layer_next = [layer_next; U];
            U = p1 + first_layer_thickness*q5;
            layer_next = [layer_next; U];
        else
            U = p1 + first_layer_thickness*q2;
            layer_next = [layer_next; U];
        end
    end
    for t = length(all_points)/2+1:length(all_points)
        p1 = all_points(t, :);
        p2 = all_points(t-1, :);
        if t == length(all_points)
            p3 = all_points(1, :);
        else
            p3 = all_points(t+1, :);
        end

        [q1, q2, q3] = findQuarterBisectors(p1, p2, p3);
        [q4, q5] = findTrisects(p1, p3, p2);
        if t == length(all_points)/2+1
            U = p1 + first_layer_thickness*q4;
            layer_next = [layer_next; U];
            U = p1 + first_layer_thickness*q2;
            layer_next = [layer_next; U];
            U = p1 + first_layer_thickness*q5;
            layer_next = [layer_next; U];
        else
            U = p1 + first_layer_thickness*q2;
            layer_next = [layer_next; U];
        end
    end
    layer_next = circshift(layer_next, -1); 
    
    elements = [];
    boundaries = [];

    element = [];
    element(1,:) = layer_next(1, :);
    element(2,:) = layer_next(end, :);
    element(3,:) = all_points(1, :);
    element(4,:) = layer_next(2, :);
    elements(end+1, :, :) = element;
    for t=2:length(all_points)/2+1
        element = [];
        element(1,:) = layer_next(t+1, :);
        element(2,:) = layer_next(t, :);
        element(3,:) = all_points(t-1, :);
        element(4,:) = all_points(t, :);
        elements(end+1, :, :) = element;
        boundaries(end+1, :) = [size(elements, 1); 1;];
    end
    element = [];
    element(1,:) = layer_next(length(all_points)/2+3, :);
    element(2,:) = layer_next(length(all_points)/2+2, :);
    element(3,:) = all_points(length(all_points)/2+1, :);
    element(4,:) = layer_next(length(all_points)/2+4, :);
    elements(end+1, :, :) = element;
    for t=length(all_points)/2+4:length(all_points)+3
        tm2 = t-2;
        if t == length(all_points)+3
            tm2 = 1;
        end
        element = [];
        element(1,:) = layer_next(t+1, :);
        element(2,:) = layer_next(t, :);
        element(3,:) = all_points(t-3, :);
        element(4,:) = all_points(tm2, :);
        elements(end+1, :, :) = element;
        boundaries(end+1, :) = [size(elements, 1); 1;];
    end

    all_points = layer_next;





    % c1=[all_points(1,2)*(1-t)+t*(-R_a);all_points(1,3)*(1-t)+t*(0)];        %left
    % c1=zeros(2,ny);
    c2=[all_points(:,2)';all_points(:,3)'];          %bottom

    c4=[[-R_a+zeros(1,k_inner)';linspace(-R_a,R_a,n_top+2+1)(1:end-1)'; R_a+zeros(1,k_inner)'; ...
          R_a+zeros(1,k_inner)';linspace(R_a,-R_a,n_top+2+1)(1:end-1)';-R_a+zeros(1,k_inner)';...
        ]'; ...
        [linspace(0,R_t,k_inner+1)(1:end-1)';R_t+zeros(1,n_top+2)'; ...
         linspace(R_t,R_b,2*k_inner+1)(1:end-1)';R_b+zeros(1,n_top+2)';linspace(R_b,0,k_inner+1)(1:end-1)']'];        %top
    % c3=[all_points(end,2)*(1-t)+t*(c4(1,end));all_points(end,3)*(1-t)+t*(c4(2,end))]; %right
    % c3=zeros(2,ny);

    alpha=zeros(nx,ny);
    beta=zeros(nx,ny);
    gamma=zeros(nx,ny);
    X=zeros(nx,ny);
    i = 1:nx;
    j = 1:ny;
    [I, J] = meshgrid(i, j);  % J is ny×nx, I is ny×nx

    % Interpolate in the vertical direction
    X = (1 - (J' - 1) / (ny - 1)) .* reshape(c2(1,I'), [nx, ny]) + ((J' - 1) / (ny - 1)) .* reshape(c4(1,I'), [nx, ny]);
    Y = (1 - (J' - 1) / (ny - 1)) .* reshape(c2(2,I'), [nx, ny]) + ((J' - 1) / (ny - 1)) .* reshape(c4(2,I'), [nx, ny]);
    newX=X; newY=Y;
    Er1=zeros(1,maxit);
    Er2=zeros(1,maxit);

    weight(i,j) = repmat(1 + 9 * exp(-j/ny*10), nx, 1);

    for t=1:maxit
        i = 1:nx;
        j = 2:ny-1;

        im = i - 1;
        ip = i + 1;

        im(im < 1) = nx;
        ip(ip > nx) = 1;

        x_eta = (X(i,j+1)-X(i,j-1))/2;
        x_xi = (X(ip,j)-X(im,j))/2;
        x_xixi = X(ip,j)+X(im,j)-2*X(i,j);
        x_etaeta = X(i,j+1)+X(i,j-1)-2*X(i,j);
        y_eta = (Y(i,j+1)-Y(i,j-1))/2;
        y_xi = (Y(ip,j)-Y(im,j))/2;
        y_xixi = Y(ip,j)+Y(im,j)-2*Y(i,j);
        y_etaeta = Y(i,j+1)+Y(i,j-1)-2*Y(i,j);

        alpha(i,j)=x_eta.^2 + y_eta.^2;
        beta(i,j)=x_xi.*x_eta + y_xi.*y_eta;
        gamma(i,j)=x_xi.^2 + y_xi.^2;
        
        newX(i,j)=(-0.5)./(alpha(i,j)+gamma(i,j)+10^-9).*(...
            2*beta(i,j).*(X(ip,j+1)-X(im,j+1)-X(ip,j-1) + X(im,j-1))/4 ...
            -alpha(i,j).*(X(ip,j)+X(im,j))-gamma(i,j).*(X(i,j+1)+X(i,j-1))...
            + weight(i,j).*x_xi 
            + weight(i,j).*x_eta
        );
        newY(i,j)=(-0.5)./(alpha(i,j)+gamma(i,j)+10^-9).*(...
            2*beta(i,j).*(Y(ip,j+1)-Y(im,j+1)-Y(ip,j-1) + Y(im,j-1))/4 ...
            -alpha(i,j).*(Y(ip,j)+Y(im,j))-gamma(i,j).*(Y(i,j+1)+Y(i,j-1))...
            + weight(i,j).*y_xi 
            + weight(i,j).*y_eta
        );
        Er1(1,t)=max(max(abs(newX-X)));
        Er2(1,t)=max(max(abs(newY-Y)));
        X=newX;
        Y=newY;
        % min(Er1(t),Er2(t))
        if Er1(t)<Ermax &&Er2(t)<Ermax
          break
        end
    end
    if t==maxit
        warning('convergence not reached')
    end

    % clf 
    % hold on
    % axis equal
    % for m=1:nx
    % plot(X(m,:),Y(m,:),'b');
    % end
    % for m=1:ny
    % plot(X(:,m),Y(:,m),'Color',[0 0 0]);
    % end

    for i=1:nx
    for j=2:ny
        inext = i+1;
        if i == nx
            inext = 1;
        end
        element = zeros(4,3);
        element(1,:) = [x, X(inext,j), Y(inext,j)];
        element(2,:) = [x, X(i,j), Y(i,j)];
        element(3,:) = [x, X(i,j-1), Y(i,j-1)];
        element(4,:) = [x, X(inext,j-1), Y(inext,j-1)];
        elements(end+1, :, :) = element;
        if j == ny
            if i > k_inner && i < k_inner + n_top + 3
                boundaries(end+1, :) = [size(elements, 1); 2;];
            end
            if i > 3*k_inner + n_top + 2 && i < 3*k_inner + 2*n_top + 5
                boundaries(end+1, :) = [size(elements, 1); 3;];
            end
        end
    end
    end
    checkCounterClockwise(elements)
end

function checkCounterClockwise(elements)
    for k = 1:size(elements, 1)
        element = squeeze(elements(k, :, :));
        for i = 1:4
            v1 = element(i, :);
            v2 = element(mod(i, 4) + 1, :);
            v3 = element(mod(i + 1, 4) + 1, :);
            edge1 = v2 - v1;
            edge2 = v3 - v2;
            
            cross_prod = cross(edge1, edge2);
            
            if cross_prod(1) <= 0
                fprintf('Element %d is not counterclockwise at corner %d.\n', k, i);
            end
        end
    end
end
function bisector = findBisect(p1, p2, p3)
    vec1 = p2 - p1;
    vec2 = p3 - p1;
    vec1 = vec1 / norm(vec1);
    vec2 = vec2/ norm(vec2);
    
    bisector = vec1 + vec2;
    bisector = bisector' / norm(bisector);
    
    cross_prod = cross([vec1], [vec2]);
    
    if cross_prod(1) > 0
        bisector = -bisector;
    end
    bisector = bisector';
end
function [q1, q2, q3] = findQuarterBisectors(p1, p2, p3)
    vec1 = p2 - p1;
    vec2 = p3 - p1;
    vec1 = vec1 / norm(vec1);
    vec2 = vec2 / norm(vec2);

    % First bisector (between vec1 and vec2)
    bisector = vec1 + vec2;
    bisector = bisector / norm(bisector);
    cross_prod = cross(vec1, vec2);
    if cross_prod(1) > 0
        bisector = -bisector;
    end

    % First quarter: between vec1 and bisector
    q1 = vec1 + bisector;
    q1 = q1 / norm(q1);

    % Third quarter: between bisector and vec2
    q3 = bisector + vec2;
    q3 = q3 / norm(q3);

    % Ensure correct orientation (optional, similar to your original code)

    % Return as row vectors
    q1 = q1(:)';
    q2 = bisector(:)';
    q3 = q3(:)';
end

function [trivec1, trivec2]=findTrisects(p1, p2, p3)
    vec1 = p2 - p1;
    vec2 = p3 - p1;
    vec1 = vec1 / norm(vec1);
    vec2 = vec2 / norm(vec2);
    
    angle = atan2d(vec2(3), vec2(2)) - atan2d(vec1(3), vec1(2));
    if angle < 0
        angle = angle + 360; 
    end
    
    angle2 = angle / 3;
    angle1 = 2 * angle / 3;
    
    trivec1 = rotateVector(vec1, angle1)';
    trivec2 = rotateVector(vec1, angle2)';
end

function rotated_vec = rotateVector(vec, angle)
    rotation_matrix = [
        0, 0, 0;
        0, cosd(angle), -sind(angle); 
        0, sind(angle), cosd(angle)];
    rotated_vec = rotation_matrix * vec(:);
end
