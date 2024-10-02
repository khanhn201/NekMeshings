% Internal boundary: 2
% faceData = [];
% for i = 1:N
%   for face = 1:6
%     coordinates = reshape(elements(i, face_node_map(face, :), :), [4, 3]);
%     coordinates = sortrows(coordinates, [1 2 3]);
%     faceData = [faceData; [reshape(coordinates, [1, 12]), i, face]];
%   end
% end
% faceData = sortrows(faceData, [1 2 3 4 5 6 7 8 9 10 11 12]);
% max_width = 10;
% for i = 1:(size(faceData,1)-max_width)
%     for k = 1:max_width
%       element_i = faceData(i, 13);
%       face_i = faceData(i, 14);
%       coordinates_i = faceData(i, 1:12);
%       element_j = faceData(i+k, 13);
%       face_j = faceData(i+k, 14);
%       coordinates_j = faceData(i+k, 1:12);
%       if isequal(coordinates_i, coordinates_j)
%           boundaries(element_i, face_i, 1) = 2;
%           boundaries(element_i, face_i, 2) = element_j;
%           boundaries(element_i, face_i, 3) = face_j;
%           boundaries(element_j, face_j, 1) = 2;
%           boundaries(element_j, face_j, 2) = element_i;
%           boundaries(element_j, face_j, 3) = face_i;
%           break
%       endif
%   end
% end
