filename = 'slices/slice0049.txt';
filename = 'slices/slice0000.txt';
slice = readSliceFile(filename);
[pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
[elementsOuter, boundariesOuter] = meshOuter(pp, arc_length, arc_length_at_max_y);
[elementsInner, boundariesInner] = meshInner(pp, arc_length, arc_length_at_max_y);
elements = [elementsOuter; elementsInner;];
% boundaries = boundariesOuter;
plotElements(elements, boundaries);


% Mesh slices
xs = [];
sliceElements = [];
sliceBoundaries = [];
for i = 0:8:208
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elements, boundaries] = meshOuter(pp, arc_length, arc_length_at_max_y);
    sliceElements(end+1, :, :, :) = elements;
    sliceBoundaries = boundaries;
    xs(end+1) = slice(1,1);
end
size(sliceElements)
if ~ismember(0, xs) || norm(xs + flip(xs)) > 1e-6
    disp('Error: slices do not contain x=0 or not symmetric around x=0')
end

% Connect slices
elements = [];
boundaries = [];
[numSlices, numElements, numVertices, dim] = size(sliceElements);
for k = 1:(numSlices - 1)
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
                boundaries(end+1, :) = [size(elements,1); 3; 1];
            end
            if tag == 2 || tag == 3
                boundaries(end+1, :) = [size(elements,1); 1; 2];
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
        end
    end
end



% Add end caps
config
for i = 0:0
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elementsOuter, boundariesOuter] = meshOuter(pp, arc_length, arc_length_at_max_y);
    [elementsInner, boundariesInner] = meshInner(pp, arc_length, arc_length_at_max_y);
    for elem = 1:size(elementsInner, 1)
        layer_k = squeeze(elementsInner(elem,:, :)); 
        layer_k1 = squeeze(elementsInner(elem,:, :)); 
        layer_k1(:,1) = -R_x;
        element = [layer_k1; layer_k]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        boundaries(end+1, :) = [size(elements,1); 6; 1];
        boundaries(end+1, :) = [size(elements,1); 5; 2];
    end
    for elem = 1:size(elementsOuter, 1)
        layer_k = squeeze(elementsOuter(elem,:, :)); 
        layer_k1 = squeeze(elementsOuter(elem,:, :)); 
        layer_k1(:,1) = -R_x;
        element = [layer_k1; layer_k]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        boundaries(end+1, :) = [size(elements,1); 5; 2];

        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx, 2);
            if tag == 2 || tag == 3
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
        end
    end
end

for i = 208:208
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elementsOuter, boundariesOuter] = meshOuter(pp, arc_length, arc_length_at_max_y);
    [elementsInner, boundariesInner] = meshInner(pp, arc_length, arc_length_at_max_y);
    for elem = 1:size(elementsInner, 1)
        layer_k = squeeze(elementsInner(elem,:, :)); 
        layer_k1 = squeeze(elementsInner(elem,:, :)); 
        layer_k1(:,1) = R_x;
        element = [layer_k; layer_k1]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        boundaries(end+1, :) = [size(elements,1); 5; 1];
        boundaries(end+1, :) = [size(elements,1); 6; 2];
    end
    for elem = 1:size(elementsOuter, 1)
        layer_k = squeeze(elementsOuter(elem,:, :)); 
        layer_k1 = squeeze(elementsOuter(elem,:, :)); 
        layer_k1(:,1) = R_x;
        element = [layer_k; layer_k1]; % 8 x 3
        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        boundaries(end+1, :) = [size(elements,1); 6; 2];

        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx, 2);
            if tag == 2 || tag == 3
                boundaries(end+1, :) = [size(elements,1); 1; 2];
            end
        end
    end
end

% Wrap cylinder
xs = squeeze(xs);
xs = xs(:);
half_pos = floor(size(xs,1)/2)+1;
xs = xs(half_pos:end);
xs(end+1) = R_x;
[cylElements, cylBoundaries] = wrapCylinder(xs);
config
zs = linspace(R_b, R_t, k_inner*2)(:);
for k = 2:size(zs,1)
    z_prev = zs(k-1);
    z = zs(k);
    for elem = 1:size(cylElements,1)
        cylElemk = squeeze(cylElements(elem, :, :));
        cylElemk(:, 3) = z_prev;
        cylElemk1 = squeeze(cylElements(elem, :, :));
        cylElemk1(:, 3) = z;
        element = [cylElemk; cylElemk1];

        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, cylBoundaries(:, 1));
        if isBoundary
            boundaries(end+1, :) = [size(elements,1); 3; 3];
        end
        if k == 2
            boundaries(end+1, :) = [size(elements,1); 5; 3];
        end
        if k == size(zs,1)
            boundaries(end+1, :) = [size(elements,1); 6; 3];
        end
    end
end

[cylElements, cylBoundaries] = wrapCylinder(xs);
cylElements(:, :, 2) = -cylElements(:, :, 2);
zs = linspace(R_b, R_t, k_inner*2)(:);
for k = 2:size(zs,1)
    z_prev = zs(k-1);
    z = zs(k);
    for elem = 1:size(cylElements,1)
        cylElemk = squeeze(cylElements(elem, :, :));
        cylElemk(:, 3) = z_prev;
        cylElemk1 = squeeze(cylElements(elem, :, :));
        cylElemk1(:, 3) = z;
        element = [cylElemk1; cylElemk];

        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, cylBoundaries(:, 1));
        if isBoundary
            boundaries(end+1, :) = [size(elements,1); 3; 3];
        end
        if k == 2
            boundaries(end+1, :) = [size(elements,1); 6; 3];
        end
        if k == size(zs,1)
            boundaries(end+1, :) = [size(elements,1); 5; 3];
        end
    end
end

size(elements)
size(boundaries)
% plotElements3D(elements)

exportREA("output.rea", elements, boundaries)
plotBC(elements, boundaries)
