addpath(fullfile('..','utils'));
[elements, boundaries] = importRE2('input.re2');

% Swap y and z
tmp = elements(:,:,2);
elements(:,:,2) = elements(:,:,3);
elements(:,:,3) = tmp;
elements(:,:,1) = -elements(:,:,1);

% Change the side condition on y to W
idx = (boundaries(:,3) == 7);
boundaries(idx,3) = 1;

% Extrude the top layer upward
ftonMap = getFaceToNodeMap(); 
idx = find(boundaries(:,3) == 6);

elemToExtrude = boundaries(idx, 1);
faceToExtrude = boundaries(idx, 2);
numFaces = length(idx);
faceElements = zeros(numFaces, 4, 3);
boundaryFaces = cell(numFaces, 1);
for i = 1:numFaces
    face = faceToExtrude(i);
    nodes = ftonMap(face, :);
    faceElements(i, :, :) = elements(elemToExtrude(i), nodes, :);
    
    % Find all boundary faces for this element (excluding the extruded face)
    elementBoundaries = boundaries(boundaries(:, 1) == elemToExtrude(i) & boundaries(:, 2) ~= faceToExtrude(i), 2);
    boundaryFaces{i} = elementBoundaries;
end

dy = 0.25;
numLayers = 12;
totalNewElements = numFaces * numLayers;

newElements = zeros(totalNewElements, 8, 3);
totalBoundaryConditions = 0;
for i = 1:numFaces
    totalBoundaryConditions = totalBoundaryConditions + length(boundaryFaces{i});
end
totalBoundaryConditions = totalBoundaryConditions * numLayers + numFaces;

newBoundaries = zeros(totalBoundaryConditions, 3);

% Create all layers at once
elementIdx = 1;
boundaryIdx = 1;

faceMap = [1, 2, 2, 4, 4, 4];
for k = 1:numLayers
    k
    face1 = faceElements;
    face2 = faceElements;
    face1(:, :, 2) = face1(:, :, 2) + (k-1) * dy;
    face2(:, :, 2) = face2(:, :, 2) + k * dy;
    
    newElements(elementIdx:elementIdx+numFaces-1, 1:4, :) = face1;
    newElements(elementIdx:elementIdx+numFaces-1, 5:8, :) = face2;

    currentElementNumbers = size(elements, 1) + (elementIdx:elementIdx+numFaces-1);
    
    for i = 1:numFaces
        for face=boundaryFaces{i}'
            newBoundaries(boundaryIdx, 1) = currentElementNumbers(i);
            newBoundaries(boundaryIdx, 2) = faceMap(face);
            newBoundaries(boundaryIdx, 3) = 1;
            boundaryIdx = boundaryIdx + 1;
        end
    end
    if k == numLayers
        newBoundaries(boundaryIdx:boundaryIdx+numFaces-1, 1) = currentElementNumbers';
        newBoundaries(boundaryIdx:boundaryIdx+numFaces-1, 2) = 6;
        newBoundaries(boundaryIdx:boundaryIdx+numFaces-1, 3) = 1;
        boundaryIdx = boundaryIdx + numFaces;
    end

    
    elementIdx = elementIdx + numFaces;
end

elements = [elements; newElements];
boundaries = [boundaries; newBoundaries];

boundaries(idx,:) = [];




exportVTK('test.vtk', elements, boundaries)
exportRE2('cbarr', elements, boundaries)
