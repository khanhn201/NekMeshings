% Mesh slices
% slicesCoord = readSlices('iea15.mat');
slicesCoord = readSlices('nrel5mw2.mat');
config
zss = slicesCoord(:, 1, 3);
zs = linspace(z_shift, zss(end), n_slices);
slicesCoord = interp1(zss, slicesCoord, zs);

zs = [];
sliceElements = [];
sliceBoundaries = [];
sliceSplines = {};

[elements, boundaries, pp_coarse] = meshHub();
z = 0;
sliceElements(end+1, :, :, :) = elements;
sliceBoundaries = boundaries;
zs(end+1) = 0;
sliceSplines{size(sliceElements,1)} = pp_coarse;
for i = 1:n_slices
    slice = squeeze(slicesCoord(i, :, :));
    slice(:, 3) = slice(1, 3);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elements, boundaries, pp_coarse, X, Y, first_layer] = meshOuterElliptic(pp, arc_length, arc_length_at_max_y);

    z = slice(1,3);
    % [elements] = smoothMesh(elements, boundaries);
    sliceElements(end+1, :, :, :) = elements;
    sliceBoundaries = boundaries;
    zs(end+1) = z;
    sliceSplines{size(sliceElements,1)} = pp_coarse;

    endcap_slice = slice;
end
% End cap outer
% for j = 1:length(R_end_caps)
%     slice = endcap_slice;
%     R_x = R_end_caps(j);
%     slice(:,3) = R_x;
%     [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
%     [elements, boundaries, pp_coarse] = meshOuterElliptic(pp, arc_length, arc_length_at_max_y);
%     sliceElements(end+1, :, :, :) = elements;
%     sliceBoundaries = boundaries;
%     zs(end+1) = R_x;
%     % sliceSplines{size(sliceElements,1)} = pp_coarse;
% end

% Add connecting splines
connectingSplines = {};
for i = 1:size(sliceSplines{1}.breaks, 2)
    splinePoints = [];
    z = [];
    for j=1:length(sliceSplines)
        splinePoints(end+1, :) = ppval(sliceSplines{j}, sliceSplines{j}.breaks(i));
        z(end+1) = splinePoints(end, 3);
    end
    splinePoints = [[0;0;1]';splinePoints;[0;0;1]'];
    pp_coarse = spline(z, splinePoints');
    connectingSplines{i} = pp_coarse;
end
% plotSplines
% asdfasf(sadfasdf)



% Connect slices
disp("connect slices")
elements = [];
boundaries = [];
surfaces = struct('elem', {}, 'face', {}, 'splineCoeffs', {}, 'splineEnds', {});
[numSlices, numElements, numVertices, dim] = size(sliceElements);
for k = 1:(numSlices - 1)
    k
% for k = 1:1
    count_wall = 0;
    for elem = 1:numElements
        vertices_layer_k = squeeze(sliceElements(k, elem, :, :)); % 4 x 3
        vertices_layer_k1 = squeeze(sliceElements(k+1, elem, :, :)); % 4 x 3
        element = [vertices_layer_k; vertices_layer_k1]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx, 2);
            if tag == 1
                if k < numSlices
                    spline1 = sliceSplines{k};
                    spline3 = sliceSplines{k+1};
                    boundaries(end+1, :) = [size(elements,1); 3; 1];

                    count_wall = count_wall + 1;
                    spline2 = connectingSplines{count_wall+1};
                    spline4 = connectingSplines{count_wall};

                    j = count_wall;
                    spline1piece = spline1.coefs((j-1)*3 + 1: (j-1)*3 + 3,:);
                    spline1start = spline1.breaks(j);
                    spline1end = spline1.breaks(j + 1);

                    spline3piece = spline3.coefs((j-1)*3 + 1: (j-1)*3 + 3,:);
                    spline3start = spline3.breaks(j);
                    spline3end = spline3.breaks(j + 1);

                    spline2piece = spline2.coefs((k-1)*3 + 1: (k-1)*3 + 3,:);
                    spline2start = spline2.breaks(k);
                    spline2end = spline2.breaks(k+1);

                    spline4piece = spline4.coefs((k-1)*3 + 1: (k-1)*3 + 3,:);
                    spline4start = spline4.breaks(k);
                    spline4end = spline4.breaks(k+1);
                    
                    if k > 1
                        surfaceStruct.elem = size(elements,1);
                        surfaceStruct.splineCoeffs = zeros(4, 3, 4);
                        surfaceStruct.splineCoeffs(1,:,:) = spline1piece;
                        surfaceStruct.splineCoeffs(2,:,:) = spline2piece;
                        surfaceStruct.splineCoeffs(3,:,:) = spline3piece;
                        surfaceStruct.splineCoeffs(4,:,:) = spline4piece;
                        surfaceStruct.splineEnds = zeros(4, 2);
                        surfaceStruct.splineEnds(1,:) = [spline1start, spline1end];
                        surfaceStruct.splineEnds(2,:) = [spline2start, spline2end];
                        surfaceStruct.splineEnds(3,:) = [spline3start, spline3end];
                        surfaceStruct.splineEnds(4,:) = [spline4start, spline4end];
                        surfaceStruct.face = 3;
                        surfaces(end+1) = surfaceStruct;
                    end
                end
            end
            if tag == 2 
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
            if tag == 3 && R_downstream == 0
                boundaries(end+1, :) = [size(elements,1); 1; 4];
            end
        end
        % if k == numSlices - 1
        %     boundaries(end+1, :) = [size(elements,1); 6; 3];
        % end
    end
end



layer_count = size(X, 2);
layer_size = size(X, 1);
z = z(end);
Xmod = X(:, :);
Ymod = Y(:, :);
Xmod(1, 1:end-1) = (Xmod(1, 1:end-1) + Xmod(1, 2:end))/2;
Ymod(1, 1:end-1) = (Ymod(1, 1:end-1) + Ymod(1, 2:end))/2;
Xmod(layer_size/2+1, 1:end-1) = (Xmod(layer_size/2+1, 1:end-1) + Xmod(layer_size/2+1, 2:end))/2;
Ymod(layer_size/2+1, 1:end-1) = (Ymod(layer_size/2+1, 1:end-1) + Ymod(layer_size/2+1, 2:end))/2;

z_diag = [z+first_layer_thickness];
for j = 2:layer_count
    m = 1.9^(j-1)*first_layer_thickness;
    mp = 1.9^(j-2)*first_layer_thickness;
    for i = 1:layer_size
        inext = i + 1;
        if i == layer_size
            inext = 1;
        end
        element = zeros(8, 3);
        element(1,:) = [X(inext,j),   Y(inext,j), z];
        element(2,:) = [X(i,j),       Y(i,j), z];
        element(3,:) = [X(i,j-1),     Y(i,j-1), z];
        element(4,:) = [X(inext,j-1), Y(inext,j-1), z];


        element(5,:) = [Xmod(inext,j),   Ymod(inext,j), z+m];
        element(6,:) = [Xmod(i,j),       Ymod(i,j), z+m];
        element(7,:) = [Xmod(i,j-1),     Ymod(i,j-1), z+mp];
        element(8,:) = [Xmod(inext,j-1), Ymod(inext,j-1), z+mp];
        if j==layer_count
            element(5,:) = [Xmod(inext,j),   Ymod(inext,j), z+mp];
            element(6,:) = [Xmod(i,j),       Ymod(i,j), z+mp];
        end
        checkLeftHanded(element);
        elements(end+1, :, :) = element;

        if j==layer_count
            boundaries(end+1, :) = [size(elements,1); 6; 3];
            if i >= k_outer -2 && i < k_outer + n_top
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
            if i >= n_top + 2*k_outer + 2 && i <=  2*n_top + 2*k_outer + 3
                boundaries(end+1, :) = [size(elements,1); 1; 4];
            end
        end
    end
    z_diag(end+1) = z+m;
end


% Add end caps inner
[pp, arc_length, arc_length_at_max_y] = fitSpline(endcap_slice);
[elementsInner, boundariesInner] = meshInnerRec(pp, arc_length);
inner_size = size(elementsInner, 1);
elementsInner = [elementsInner; first_layer];

for j = 1:layer_count-1
% for j = 1:1
    p_coord = [X(:, j), Y(:, j), z*ones(layer_size,1)];
    if j > 1
        p_coord = [Xmod(:, j-1), Ymod(:, j-1), z_diag(j-1)*ones(layer_size,1)];
    end
    ux = solveLaplace(elementsInner, p_coord, Xmod(:,j)-p_coord(:,1));
    uy = solveLaplace(elementsInner, p_coord, Ymod(:,j)-p_coord(:,2));
    elementsInnerNext = elementsInner;
    elementsInnerNext(:, :, 1) += ux;
    elementsInnerNext(:, :, 2) += uy;
    elementsInnerNext(:, :, 3) = z_diag(j);
    p_coord = [Xmod(:, j), Ymod(:, j), z_diag(j)*ones(layer_size,1)];
    elementsInnerNext = relaxQuadMesh(elementsInnerNext, p_coord, 500);

    for elem = 1:size(elementsInner, 1)
        layer_k = squeeze(elementsInner(elem,:, :)); 
        layer_k1 = squeeze(elementsInnerNext(elem,:, :)); 

        element = [layer_k; layer_k1]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        if j == 1 && elem <= inner_size
            boundaries(end+1, :) = [size(elements,1); 5; 1];
        end
        if j == layer_count-1
            boundaries(end+1, :) = [size(elements,1); 6; 3];
        end
    end
    elementsInner = elementsInnerNext;
end
zs(end+1) = z_diag(end-1);

% Wrap cylinder
zs = squeeze(zs);
zs = zs(:);
% xs = [xs; R_end_caps'];
[cylElements, cylBoundaries] = wrapFanDiamond(zs);
config
ys = linspace(R_b, R_t, k_inner*2 + 1)(:);
for k = 2:size(ys,1)
    y_prev = ys(k-1);
    y = ys(k);
    for elem = 1:size(cylElements,1)
        cylElemk = squeeze(cylElements(elem, :, :));
        cylElemk(:, 2) = y_prev;
        cylElemk1 = squeeze(cylElements(elem, :, :));
        cylElemk1(:, 2) = y;
        element = [cylElemk1; cylElemk];

        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, cylBoundaries(:, 1));
        if isBoundary
            boundaries(end+1, :) = [size(elements,1); 1; 3];
        end
        if k == 2 && R_downstream == 0
            boundaries(end+1, :) = [size(elements,1); 6; 4];
        end
        if k == size(ys,1)
            boundaries(end+1, :) = [size(elements,1); 5; 2];
        end
    end
end
% plotElements3D(elements)

% Extends downstream
if R_downstream < 0
    r = (R_b / R_downstream)^(1 / (k_downstream - 1));  % common ratio
    ys = R_downstream * r.^((0:k_downstream-1))(:);
    % ys = linspace(R_downstream, R_b, k_downstream)(:);
    for k = 2:size(ys,1)
        y_prev = ys(k-1);
        y = ys(k);
        for elem = 1:size(cylElements,1)
            cylElemk = squeeze(cylElements(elem, :, :));
            cylElemk(:, 2) = y_prev;
            cylElemk1 = squeeze(cylElements(elem, :, :));
            cylElemk1(:, 2) = y;
            element = [cylElemk1; cylElemk];

            checkLeftHanded(element);

            elements(end+1, :, :) = element;
            [isBoundary, idx] = ismember(elem, cylBoundaries(:, 1));
            if isBoundary
                boundaries(end+1, :) = [size(elements,1); 1; 3];
            end
            if k == 2
                boundaries(end+1, :) = [size(elements,1); 6; 4];
            end
        end
    end
    [gridElements, gridBoundaries] = gridBlade(zs);
    for k = 2:size(ys,1)
        y_prev = ys(k-1);
        y = ys(k);
        for elem = 1:size(gridElements,1)
            cylElemk = squeeze(gridElements(elem, :, :));
            cylElemk(:, 2) = y_prev;
            cylElemk1 = squeeze(gridElements(elem, :, :));
            cylElemk1(:, 2) = y;
            element = [cylElemk1; cylElemk];

            checkLeftHanded(element);

            elements(end+1, :, :) = element;
            [isBoundary, idx] = ismember(elem, gridBoundaries(:, 1));
            if isBoundary
                boundaries(end+1, :) = [size(elements,1); 1; 3];
            end
            if k == 2
                boundaries(end+1, :) = [size(elements,1); 6; 4];
            end
        end
    end
end

% plotBC(elements, boundaries)

% Replicate 2 more rotated copies
N_elem = size(elements,1)
groupElements = zeros(3 * N_elem, 8, 3);
groupElements(1:N_elem,:,:) = elements;
groupBoundaries = zeros(3 * size(boundaries,1), 3);
groupBoundaries(1:size(boundaries,1), :) = boundaries;

groupSurfaces = struct('elem', {}, 'face', {}, 'splineCoeffs', {}, 'splineEnds', {});
for k=1:2
    R = [cos(2*k*pi/3), 0, -sin(2*k*pi/3);
         0,           1,            0;
         sin(2*k*pi/3), 0, cos(2*k*pi/3)];
    for i=1:N_elem
        element = squeeze(elements(i, :, :));
        rotated = (R * element')';
        groupElements(k*N_elem + i, :, :) = rotated;
    end
    offset = k * N_elem;
    b_start = k * size(boundaries, 1) + 1;
    b_end = (k+1) * size(boundaries, 1);
    groupBoundaries(b_start:b_end, :) = boundaries;
    groupBoundaries(b_start:b_end, 1) += offset;
end
for k=1:3
    R = [cos(2*(k-1)*pi/3), 0, -sin(2*(k-1)*pi/3);
         0,           1,            0;
         sin(2*(k-1)*pi/3), 0, cos(2*(k-1)*pi/3)];
    offset = (k - 1) * N_elem;
    for i=1:length(surfaces)
        surfaceStruct.elem = surfaces(i).elem + offset;
        surfaceStruct.face = surfaces(i).face;
        surfaceStruct.splineEnds = surfaces(i).splineEnds;
        coeffs = surfaces(i).splineCoeffs;
        for j = 1:4
           surfaceStruct.splineCoeffs(j,:,:) = R * squeeze(coeffs(j,:,:)); 
        end
        groupSurfaces((k-1)*length(surfaces) + i) = surfaceStruct;
    end
end
size(groupElements)
size(groupBoundaries)
size(groupSurfaces)


% rmdir('./surfaces/', 's')
% mkdir('./surfaces/')
% for i=1:length(groupSurfaces)
%     element = groupSurfaces(i).elem;
%     filename = sprintf('surfaces/surface%08d.txt', element);
%     fid = fopen(filename, "w");
%     for j = 1:4
%         for k = 1:3
%             fprintf(fid, '%15.7g %15.7g %15.7g %15.7g\n', groupSurfaces(i).splineCoeffs(j,k,:));
%         end
%         fprintf(fid, '%15.7g %15.7g\n', groupSurfaces(i).splineEnds(j,1), groupSurfaces(i).splineEnds(j,2));
%     end
%     fclose(fid);
% end

% fid = fopen('all_surfaces.txt', 'w');
% for i = 1:length(groupSurfaces)
%     element = groupSurfaces(i).elem;
%     fprintf(fid, 'ELEMENT %d\n', element);  % Optional header for each element
%     for j = 1:4
%         for k = 1:3
%             fprintf(fid, '%15.7g %15.7g %15.7g %15.7g\n', groupSurfaces(i).splineCoeffs(j,k,:));
%         end
%         fprintf(fid, '%15.7g %15.7g\n', groupSurfaces(i).splineEnds(j,1), groupSurfaces(i).splineEnds(j,2));
%     end
%
%     fprintf(fid, '\n');  % Optional spacing between elements
% end
% fclose(fid);

zs(end)
exportSSURF("inner", groupSurfaces);
% exportREA("turbineInner.rea", groupElements, groupBoundaries)
exportRE2("inner", groupElements, groupBoundaries);
% exportToVTK("inner.vtk", groupElements);

N = size(groupElements,1);
X = permute(groupElements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk(X,Hexes, groupBoundaries,'')
% plotBC(groupElements, groupBoundaries)
