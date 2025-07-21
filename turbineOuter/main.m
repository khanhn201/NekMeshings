[sliceElements,sliceBoundaries] = meshTurbineOuter();
% plotElements(sliceElements, sliceBoundaries);

config

k_vec = 1:k_upstream;
r = (mult_stream.^k_vec - mult_stream) / (mult_stream^k_upstream - mult_stream);
xs = (R_upstream - R_t) * r + R_t;

numElems = size(sliceElements, 1);
numLayers = numel(xs);
elements = zeros((numLayers - 1) * numElems, 8, 3);

layer0 = repmat(sliceElements, [numLayers - 1, 1, 1]); % [(k-1)*N x 4 x 3]
layer1 = repmat(sliceElements, [numLayers - 1, 1, 1]);

x_prev_vals = xs(1:end-1);
x_vals = xs(2:end);
x_prev_full = repelem(x_prev_vals(:), numElems, 4);
x_full = repelem(x_vals(:), numElems, 4);
layer0(:, :, 1) = x_prev_full;
layer1(:, :, 1) = x_full;
elements = cat(2, layer0, layer1);  % [numHexes x 8 x 3]

boundaries = [];

[isBoundary, idxs] = ismember((1:numElems)', sliceBoundaries(:, 1));
tags = zeros(numElems, 1);
tags(isBoundary) = sliceBoundaries(idxs(isBoundary), 2);

for k = 2:numLayers
    hexStart = (k-2)*numElems + 1;
    hexEnd = (k-1)*numElems;
    curIdx = hexStart:hexEnd;

    if k == 2 && any(tags == 1)
        mask = (tags == 1);
        boundaries = [boundaries; [curIdx(mask)', repmat([5, 5], nnz(mask), 1)]];
    end

    mask_outflow = ismember(tags, [3, 5, 9]);
    if any(mask_outflow)
        boundaries = [boundaries; [curIdx(mask_outflow)', repmat([2, 2], nnz(mask_outflow), 1)]];
    end

    if any(tags == 7)
        mask = (tags == 7);
        boundaries = [boundaries; [curIdx(mask)', repmat([3, 1], nnz(mask), 1)]];
    end

    if k == numLayers
        boundaries = [boundaries; [curIdx', repmat([6, 2], numElems, 1)]];
    end
end

xs = linspace(R_b, R_t, k_blade);
x_prev_vals = xs(1:end-1);
x_vals = xs(2:end);
numNewLayers = numel(x_prev_vals);
numElems = size(sliceElements, 1);

% Replicate elements
layer0 = repmat(sliceElements, [numNewLayers, 1, 1]);
layer1 = repmat(sliceElements, [numNewLayers, 1, 1]);

% Assign x coordinates
x_prev_full = repelem(x_prev_vals(:), numElems, 4);
x_full = repelem(x_vals(:), numElems, 4);
layer0(:, :, 1) = x_prev_full;
layer1(:, :, 1) = x_full;

% Form blade elements
bladeElements = cat(2, layer0, layer1);  % [M x 8 x 3]

% Identify tags for all elements
elemIDs = repmat((1:numElems)', numNewLayers, 1);
[isBoundary, idxs] = ismember(elemIDs, sliceBoundaries(:, 1));
tags = zeros(size(elemIDs));
tags(isBoundary) = sliceBoundaries(idxs(isBoundary), 2);

% Only keep if tag ≠ 1
validMask = tags ~= 1;
bladeElements = bladeElements(validMask, :, :);
tags = tags(validMask);
elemIDs = elemIDs(validMask);

curIdx = (size(elements, 1) + 1):(size(elements, 1) + size(bladeElements, 1));
curIdx = curIdx(:);
elements = cat(1, elements, bladeElements);

% Boundary tagging (same as before)
b_outflow = ismember(tags, [3, 5, 9]);
b_inlet   = tags == 7;
b_wall    = ismember(tags, [2, 4, 8]);
b_axis    = tags == 6;

boundaries = [
    boundaries;
    [curIdx(b_outflow), repmat([2, 2], nnz(b_outflow), 1)];
    [curIdx(b_inlet),   repmat([3, 1], nnz(b_inlet), 1)];
    [curIdx(b_wall),    repmat([4, 3], nnz(b_wall), 1)];
    [curIdx(b_axis),    repmat([1, 3], nnz(b_axis), 1)];
];

k = 1:k_downstream;
r = (mult_stream.^k-mult_stream)/(mult_stream^k_downstream-mult_stream);
xs = (R_downstream-R_b) * r + R_b;

numElems = size(sliceElements, 1);
numLayers = numel(xs);
elements2 = zeros((numLayers - 1) * numElems, 8, 3);

layer0 = repmat(sliceElements, [numLayers - 1, 1, 1]); % [(k-1)*N x 4 x 3]
layer1 = repmat(sliceElements, [numLayers - 1, 1, 1]);

x_prev_vals = xs(1:end-1);
x_vals = xs(2:end);
x_prev_full = repelem(x_prev_vals(:), numElems, 4);
x_full = repelem(x_vals(:), numElems, 4);
layer0(:, :, 1) = x_full;
layer1(:, :, 1) = x_prev_full;
elements2 = cat(2, layer0, layer1);  % [numHexes x 8 x 3]

boundaries2 = [];

[isBoundary, idxs] = ismember((1:numElems)', sliceBoundaries(:, 1));
tags = zeros(numElems, 1);
tags(isBoundary) = sliceBoundaries(idxs(isBoundary), 2);

for k = 2:numLayers
    hexStart = (k-2)*numElems + 1;
    hexEnd = (k-1)*numElems;
    curIdx = (size(elements, 1) + hexStart):(size(elements, 1) + hexEnd);

    if k == 2 && any(tags == 1)
        mask = (tags == 1);
        boundaries2 = [boundaries2; [curIdx(mask)', repmat([6, 5], nnz(mask), 1)]];
    end

    mask_outflow = ismember(tags, [3, 5, 9]);
    if any(mask_outflow)
        boundaries2 = [boundaries2; [curIdx(mask_outflow)', repmat([2, 2], nnz(mask_outflow), 1)]];
    end

    if any(tags == 7)
        mask = (tags == 7);
        boundaries2 = [boundaries2; [curIdx(mask)', repmat([3, 1], nnz(mask), 1)]];
    end

    if k == numLayers
        boundaries2 = [boundaries2; [curIdx', repmat([5, 4], numElems, 1)]];
    end
end

elements = cat(1, elements, elements2);
boundaries = cat(1, boundaries, boundaries2);

size(elements)
size(boundaries)
% plotElements3D(elements)

tmp = elements(:, :, 1);
elements(:, :, 1) = -elements(:, :, 2);
elements(:, :, 2) = tmp;


% exportREA("turbineOuter.rea", elements, boundaries)
% plotBC(elements, boundaries)
exportRE2("outer", elements, boundaries);
% exportToVTK("inner.vtk", groupElements);

N = size(elements,1);
X = permute(elements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk(X,Hexes, boundaries,'')
