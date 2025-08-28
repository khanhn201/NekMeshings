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

disp("mesh slices")
zs(end+1) = 0;
for i = 1:hub_layers
    i
    da = (i-1)/(hub_layers)*(z_shift-R_a/sqrt(3)) + R_a/sqrt(3);
    db = (i-1)/(hub_layers)*z_shift;
    projAngle = atan2(da - db,R_a);
    [elements, boundaries, pp_coarse] = meshHub(projAngle, db);
    sliceElements(end+1, :, :, :) = elements;
    sliceBoundaries = boundaries;
    if i > 1 % 1st layer is assumed to be 0 for wrapFanDiamond
        zs(end+1) = da;
    end
    sliceSplines{size(sliceElements,1)} = pp_coarse;
end
for i = 1:n_slices
    i
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
totalElements = (numSlices - 1) * numElements;


elements = zeros(totalElements, 8, 3);
boundaries = zeros(totalElements, 3);  % Overestimate
surfaces(totalElements) = struct('elem', [], 'face', [], 'splineCoeffs', [], 'splineEnds', []);
elemTags = zeros(numElements, 1);
elemTags(sliceBoundaries(:,1)) = sliceBoundaries(:,2);


v0 = reshape(sliceElements(1:end-1, :, :, :), [], 4, 3);  % (totalElements, 4, 3)
v1 = reshape(sliceElements(2:end, :, :, :), [], 4, 3);
elements = cat(2, v0, v1);  % (totalElements, 8, 3)

elemID = reshape(1:totalElements, numSlices-1, numElements);  % (numSlices-1 x numElements)
bndCount = 0;
surfCount = 0;
count_wall = 0;

for k = 1:(numSlices - 1)
    k
    
    tag_k = elemTags;
    elem_ids_k = elemID(k, :);

    count_wall = 0;

    wall_idx = find(tag_k == 1);  % tag 1 → wall
    if ~isempty(wall_idx)
        n = numel(wall_idx);
        bndRange = bndCount+1 : bndCount+n;
        boundaries(bndRange, :) = [elem_ids_k(wall_idx)', 3*ones(n,1), ones(n,1)];
        bndCount = bndCount + n;

        if k > 1
            spline1 = sliceSplines{k};
            spline3 = sliceSplines{k+1};
            for j = 1:n
                elem_j = elem_ids_k(wall_idx(j));
                count_wall = count_wall + 1;

                spline2 = connectingSplines{count_wall+1};
                spline4 = connectingSplines{count_wall};

                surfaces(surfCount+1).elem = elem_j;
                surfaces(surfCount+1).face = 3;
                surfaces(surfCount+1).splineCoeffs = zeros(4, 3, 4);
                surfaces(surfCount+1).splineCoeffs(1, :, :) = spline1.coefs((j-1)*3 + (1:3), :);
                surfaces(surfCount+1).splineCoeffs(2, :, :) = spline2.coefs((k-1)*3 + (1:3), :);
                surfaces(surfCount+1).splineCoeffs(3, :, :) = spline3.coefs((j-1)*3 + (1:3), :);
                surfaces(surfCount+1).splineCoeffs(4, :, :) = spline4.coefs((k-1)*3 + (1:3), :);
                surfaces(surfCount+1).splineEnds = [ ...
                    spline1.breaks(j), spline1.breaks(j+1); ...
                    spline2.breaks(k), spline2.breaks(k+1); ...
                    spline3.breaks(j), spline3.breaks(j+1); ...
                    spline4.breaks(k), spline4.breaks(k+1)];
                surfCount = surfCount + 1;
            end
        end
    end

    inflow_idx = find(tag_k == 2);
    if ~isempty(inflow_idx)
        n = numel(inflow_idx);
        bndRange = bndCount+1 : bndCount+n;
        boundaries(bndRange, :) = [elem_ids_k(inflow_idx)', ones(n,1), 2*ones(n,1)];
        bndCount = bndCount + n;
    end

    if R_downstream == 0
        outflow_idx = find(tag_k == 3);
        if ~isempty(outflow_idx)
            n = numel(outflow_idx);
            bndRange = bndCount+1 : bndCount+n;
            boundaries(bndRange, :) = [elem_ids_k(outflow_idx)', ones(n,1), 4*ones(n,1)];
            bndCount = bndCount + n;
        end
    end
end

boundaries = boundaries(1:bndCount, :);
surfaces = surfaces(1:surfCount);



disp("capping ends")

layer_count = size(X, 2);
layer_size = size(X, 1);
z = z(end);
Xmod = X(:, :);
Ymod = Y(:, :);
factors = linspace(0.25, 1, ceil(k_inner/2));
for i=1:ceil(k_inner/2)
    factor = factors(i);
    Xmod(i, 1:end-1) = factor*Xmod(i, 1:end-1) + (1-factor)*Xmod(i, 2:end);
    Ymod(i, 1:end-1) = factor*Ymod(i, 1:end-1) + (1-factor)*Ymod(i, 2:end);
    if i != 1
        Xmod(end-i+2, 1:end-1) = factor*Xmod(end-i+2, 1:end-1) + (1-factor)*Xmod(end-i+2, 2:end);
        Ymod(end-i+2, 1:end-1) = factor*Ymod(end-i+2, 1:end-1) + (1-factor)*Ymod(end-i+2, 2:end);
    end
    Xmod(layer_size/2+i, 1:end-1) = factor*Xmod(layer_size/2+i, 1:end-1) + (1-factor)*Xmod(layer_size/2+i, 2:end);
    Ymod(layer_size/2+i, 1:end-1) = factor*Ymod(layer_size/2+i, 1:end-1) + (1-factor)*Ymod(layer_size/2+i, 2:end);
    if i != 1
        Xmod(layer_size/2-i+2, 1:end-1) = factor*Xmod(layer_size/2-i+2, 1:end-1) + (1-factor)*Xmod(layer_size/2-i+2, 2:end);
        Ymod(layer_size/2-i+2, 1:end-1) = factor*Ymod(layer_size/2-i+2, 1:end-1) + (1-factor)*Ymod(layer_size/2-i+2, 2:end);
    end
end

z_diag = [z+first_layer_thickness];
for j = 2:layer_count
    m = first_layer_thickness + (Xmod(1, 1) - Xmod(1, j));
    mp = first_layer_thickness + (Xmod(1, 1) - Xmod(1, j-1));
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
        elements(end+1, :, :) = element;

        if j==layer_count
            boundaries(end+1, :) = [size(elements,1); 6; 3];
            if i >= k_inner + 1 && i <= k_inner + n_top + 2
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
            if i >= n_top + 3*k_inner + 3 && i <=  2*n_top + 3*k_inner + 4
                boundaries(end+1, :) = [size(elements,1); 1; 4];
            end
        end
    end
    z_diag(end+1) = z+m;
end


disp("capping ends inner")
% Add end caps inner
[pp, arc_length, arc_length_at_max_y] = fitSpline(endcap_slice);
[elementsInner, boundariesInner] = meshInnerRec(pp, arc_length, arc_length_at_max_y);
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
    elementsInnerNext = relaxQuadMesh(elementsInnerNext, p_coord, 10);


    numElems = size(elementsInner, 1);
    newElements = cat(2, elementsInner, elementsInnerNext);
    boundaryRows = [];
    if j == 1
        % Add [index; 5; 1] for the first inner_size elements
        ids = size(elements,1) + (1:inner_size)';
        boundaryRows = [ids, 5*ones(inner_size,1), ones(inner_size,1)];
    end

    if j == layer_count - 1
        % Add [index; 6; 3] for all elements in this layer
        ids = size(elements,1) + (1:numElems)';
        boundaryRows = [boundaryRows; ids, 6*ones(numElems,1), 3*ones(numElems,1)];
    end

    elements = cat(1, elements, newElements);
    boundaries = cat(1, boundaries, boundaryRows);

    elementsInner = elementsInnerNext;
end
zs(end+1) = z_diag(end-1);







disp("wrap cylinder")
% Wrap cylinder
zs = squeeze(zs);
zs = zs(:);
% xs = [xs; R_end_caps'];
[cylElements, cylBoundaries] = wrapFanDiamond(zs);
config
ys = linspace(R_b, R_t, k_inner*2 + 1)(:);
numElems = size(cylElements,1);
numLayers = length(ys) - 1;
layerk = reshape(cylElements, [numElems, 4, 3]);  % (N, 4, 3)
layerk = cat(2, layerk, layerk);                  % (N, 8, 3)
elements2 = zeros(numElems*numLayers, 8, 3);
boundaries2 = [];

for k = 2:length(ys)
    y_prev = ys(k-1);
    y = ys(k);
    idx = (k-2)*numElems + (1:numElems);

    % Assign the Y-coordinate values
    layerk_temp = layerk;
    layerk_temp(:,1:4,2) = y;
    layerk_temp(:,5:8,2) = y_prev;

    % Save to elements
    elements2(idx, :, :) = layerk_temp;

    [isBoundary, loc] = ismember((1:numElems)', cylBoundaries(:,1));
    if any(isBoundary)
        boundaries2 = [boundaries2;
                      [idx(isBoundary)', repmat(1, sum(isBoundary), 1), repmat(3, sum(isBoundary), 1)]];
    end
    if k == 2 && R_downstream == 0
        boundaries2 = [boundaries2;
                      [idx', repmat(6, numElems, 1), repmat(4, numElems, 1)]];
    end
    if k == length(ys)
        boundaries2 = [boundaries2;
                      [idx', repmat(5, numElems, 1), repmat(2, numElems, 1)]];
    end
end

boundaries2(:, 1) += size(elements, 1);
% elements = [elements; elements2];
% boundaries = [boundaries; boundaries2];


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

%
% disp("checking lefthand")
% for i=1:size(elements, 1)
%     checkLeftHanded(squeeze(elements(i, :, :)));
% end


% plotBC(elements, boundaries)

disp("cloning")
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


zs(end)
exportSSURF("inner", groupSurfaces);
% exportREA("turbineInner.rea", groupElements, groupBoundaries)
exportRE2("inner", groupElements, groupBoundaries);

% groupElements(:, :, :) = groupElements(:, :, :)/100;
N = size(groupElements,1);
X = permute(groupElements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk(X,Hexes, groupBoundaries,'')
% plotBC(groupElements, groupBoundaries)
