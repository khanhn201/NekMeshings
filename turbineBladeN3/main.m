rmdir('./surfaces/', 's')
mkdir('./surfaces/')

% Mesh slices
% slicesCoord = readSlices('iea15.mat');
slicesCoord = readSlices('nrel5mw.mat');

config
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

for i = 1:slice_spacing:size(slicesCoord, 1)
    i
    slice = squeeze(slicesCoord(i, :, :));
    slice(:, 3) += z_shift;
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elements, boundaries, pp_coarse] = meshOuterElliptic(pp, arc_length, arc_length_at_max_y);

    z = slice(1,3);
    % [elements] = smoothMesh(elements, boundaries);
    sliceElements(end+1, :, :, :) = elements;
    sliceBoundaries = boundaries;
    zs(end+1) = z;
    sliceSplines{size(sliceElements,1)} = pp_coarse;

    endcap_slice = slice;
end
% End cap outer
for j = 1:length(R_end_caps)
    slice = endcap_slice;
    R_x = R_end_caps(j);
    slice(:,3) = R_x;
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elements, boundaries, pp_coarse] = meshOuterElliptic(pp, arc_length, arc_length_at_max_y);
    sliceElements(end+1, :, :, :) = elements;
    sliceBoundaries = boundaries;
    zs(end+1) = R_x;
    % sliceSplines{size(sliceElements,1)} = pp_coarse;
end

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
% asdasfdas(asdas)

% Connect slices
disp("connect slices")
elements = [];
boundaries = [];
[numSlices, numElements, numVertices, dim] = size(sliceElements);
for k = 1:(numSlices - 1)
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
                if k < (numSlices - length(R_end_caps))
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

                    j = count_wall;
                    spline3piece = spline3.coefs((j-1)*3 + 1: (j-1)*3 + 3,:);
                    spline3start = spline3.breaks(j);
                    spline3end = spline3.breaks(j + 1);

                    spline2piece = spline2.coefs((k-1)*3 + 1: (k-1)*3 + 3,:);
                    spline2start = spline2.breaks(k);
                    spline2end = spline2.breaks(k+1);

                    spline4piece = spline4.coefs((k-1)*3 + 1: (k-1)*3 + 3,:);
                    spline4start = spline4.breaks(k);
                    spline4end = spline4.breaks(k+1);

                    filename = sprintf('surfaces/surface%08d.txt', size(elements,1));
                    fid = fopen(filename, "w");
                    fprintf(fid, '%15.7g %15.7g %15.7g %15.7g\n', spline1piece.');
                    fprintf(fid, '%15.7g %15.7g\n', spline1start, spline1end);
                    fprintf(fid, '%15.7g %15.7g %15.7g %15.7g\n', spline2piece.');
                    fprintf(fid, '%15.7g %15.7g\n', spline2start, spline2end);
                    fprintf(fid, '%15.7g %15.7g %15.7g %15.7g\n', spline3piece.');
                    fprintf(fid, '%15.7g %15.7g\n', spline3start, spline3end);
                    fprintf(fid, '%15.7g %15.7g %15.7g %15.7g\n', spline4piece.');
                    fprintf(fid, '%15.7g %15.7g\n', spline4start, spline4end);
                    fclose(fid);
                end
            end
            if tag == 2 
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
            if tag == 3
                boundaries(end+1, :) = [size(elements,1); 1; 4];
            end
        end
        if k == numSlices - 1
            boundaries(end+1, :) = [size(elements,1); 6; 3];
        end
    end
end



% Add end caps inner
[pp, arc_length, arc_length_at_max_y] = fitSpline(endcap_slice);
[elementsInner, boundariesInner] = meshInnerAsym(pp, arc_length);
for j = 1:length(R_end_caps)
    R_x = R_end_caps(j);
    for elem = 1:size(elementsInner, 1)
        layer_k = squeeze(elementsInner(elem,:, :)); 
        layer_k1 = squeeze(elementsInner(elem,:, :)); 
        if j > 1
            layer_k(:,3) = R_prev;
        end
        layer_k1(:,3) = R_x;
        element = [layer_k; layer_k1]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        if j == 1
            boundaries(end+1, :) = [size(elements,1); 5; 1];
        end
        if j == length(R_end_caps)
            boundaries(end+1, :) = [size(elements,1); 6; 3];
        end
    end
    R_prev = R_x;
end

% Wrap cylinder
zs = squeeze(zs);
zs = zs(:);
% xs = [xs; R_end_caps'];
[cylElements, cylBoundaries] = wrapFan(zs);
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
        if k == 2
            boundaries(end+1, :) = [size(elements,1); 6; 4];
        end
        if k == size(ys,1)
            boundaries(end+1, :) = [size(elements,1); 5; 2];
        end
    end
end
% plotElements3D(elements)


% Replicate 2 more rotated copies
N_elem = size(elements,1);
groupElements = elements;
groupBoundaries = boundaries;
for i=1:N_elem
    R = [cos(2*pi/3), 0, -sin(2*pi/3);
         0,           1,            0;
         sin(2*pi/3), 0, cos(2*pi/3)];
    element = squeeze(elements(i, :, :));
    for j=1:8
        element(j, :) = (R*squeeze(element(j, :))')';
    end
    groupElements(end+1, :, :) = element;
end
for i=1:size(boundaries, 1)
    boundary = boundaries(i, :);
    boundary(1, 1) += N_elem;
    groupBoundaries(end+1, :) = boundary;
end
for i=1:N_elem
    R = [cos(4*pi/3), 0, -sin(4*pi/3);
         0,           1,            0;
         sin(4*pi/3), 0, cos(4*pi/3)];
    element = squeeze(elements(i, :, :));
    for j=1:8
        element(j, :) = (R*squeeze(element(j, :))')';
    end
    groupElements(end+1, :, :) = element;
end
for i=1:size(boundaries, 1)
    boundary = boundaries(i, :);
    boundary(1, 1) += 2*N_elem;
    groupBoundaries(end+1, :) = boundary;
end
size(groupElements)
size(groupBoundaries)





% exportREA("turbineInner.rea", groupElements, groupBoundaries)
plotBC(groupElements, groupBoundaries)
