filename = 'slices/slice0049.txt';
slice = readSliceFile(filename);
% plotSlice(slice)
[pp, arc_length, arc_length_at_max_y] = fitSpline(slice);

[elements, boundaries] = extrudeSliceQuad(pp, arc_length, arc_length_at_max_y);
% size(elements)
size(boundaries)
plotElements(elements, boundaries);

layers = [];
for i = 0:5:100
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elements, boundaries] = extrudeSliceQuad(pp, arc_length, arc_length_at_max_y);
    layers(end+1, :, :, :) = elements;
end
size(layers)
elements = [];
% for i = 1:size(layers, 1)
% Input: layers (size: 21 x 236 x 4 x 3)
[numLayers, numElements, numVertices, dim] = size(layers);

% Iterate over consecutive layers
for k = 1:(numLayers - 1)
    for elem = 1:numElements
        vertices_layer_k = squeeze(layers(k, elem, :, :)); % 4 x 3
        vertices_layer_k1 = squeeze(layers(k+1, elem, :, :)); % 4 x 3
        object3D = [vertices_layer_k; vertices_layer_k1]; % 8 x 3
        
        elements(end+1, :, :) = object3D;
    end
end

% plotElements3D(elements)


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
