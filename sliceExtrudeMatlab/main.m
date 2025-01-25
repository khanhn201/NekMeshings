filename = 'slices/slice0049.txt';
slice = readSliceFile(filename);
% plotSlice(slice)
[pp, arc_length, arc_length_at_max_y] = fitSpline(slice);

[elements, boundaries] = extrudeSliceQuad(pp, arc_length, arc_length_at_max_y);
plotElements(elements, boundaries);

sliceElements = [];
sliceBoundaries = [];
for i = 0:10:100
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elements, boundaries] = extrudeSliceQuad(pp, arc_length, arc_length_at_max_y);
    sliceElements(end+1, :, :, :) = elements;
    sliceBoundaries = boundaries;
end
size(sliceElements)

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
                boundaries(end+1, :, :) = [size(elements,1); 3; 1];
            elseif tag == 6 
                boundaries(end+1, :, :) = [size(elements,1); 4; 2];
            else
                boundaries(end+1, :, :) = [size(elements,1); 1; 2];
            end
        end
    end
end

% plotElements3D(elements)
% plotBC(elements, boundaries)


% figure;
% hold on;
% axis equal;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('3D Objects from Consecutive Layers');
%
% for k = 1:(numLayers - 1)
%     for elem = 1:numElements
%         % Get the 4 vertices of the element in layer k
%         vertices_layer_k = squeeze(layers(k, elem, :, :)); % 4 x 3
%
%         % Get the 4 vertices of the element in layer k+1
%         vertices_layer_k1 = squeeze(layers(k+1, elem, :, :)); % 4 x 3
%
%         % Combine them to form the 3D object (8 vertices, 3 coordinates)
%         vertices = [vertices_layer_k; vertices_layer_k1]; % 8 x 3
%
%         % Define the faces of the cuboid (using vertex indices)
%         faces = [
%             1, 2, 6, 5; % Bottom face
%             2, 3, 7, 6; % Side face
%             3, 4, 8, 7; % Top face
%             4, 1, 5, 8; % Opposite side face
%             1, 2, 3, 4; % Front face
%             5, 6, 7, 8  % Back face
%         ];
%
%         % Plot the 3D object using patch
%         patch('Vertices', vertices, 'Faces', faces, ...
%               'FaceColor', 'cyan', 'EdgeColor', 'black', 'FaceAlpha', 0.5);
%     end
% end
% hold off;
