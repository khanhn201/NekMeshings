% function plotElements3D(elements)
%   figure;
%   hold on;
%   axis equal;
%   view(3); % Set 3D view
%   xlabel('X');
%   ylabel('Y');
%   zlabel('Z');
%
%   for i = 1:size(elements, 1)
%       lower_face = squeeze(elements(i, 1:4, :));
%       upper_face = squeeze(elements(i, 5:8, :));
%       fill3(lower_face(:, 1), lower_face(:, 2), lower_face(:, 3), 'r', 'FaceAlpha', 0.5);
%       fill3(upper_face(:, 1), upper_face(:, 2), upper_face(:, 3), 'b', 'FaceAlpha', 0.5);
%       for j = 1:4
%           next_idx = mod(j, 4) + 1; % Wrap around index for the next vertex
%           side_vertices = [lower_face(j, :); lower_face(next_idx, :); upper_face(next_idx, :); upper_face(j, :)];
%           fill3(side_vertices(:, 1), side_vertices(:, 2), side_vertices(:, 3), 'g', 'FaceAlpha', 0.5);
%       end
%   end
%
%   hold off;
% end

% function plotElements3D(elements)
%   figure;
%   hold on;
%   axis equal;
%   view(3); % Set 3D view
%   xlabel('X');
%   ylabel('Y');
%   zlabel('Z');
%
%   for i = 1:size(elements, 1)
%       lower_face = squeeze(elements(i, 1:4, :));
%       upper_face = squeeze(elements(i, 5:8, :));
%
%       % Plot edges of lower face
%       for j = 1:4
%           next_idx = mod(j, 4) + 1; % Wrap around index for the next vertex
%           plot3([lower_face(j, 1), lower_face(next_idx, 1)], ...
%                 [lower_face(j, 2), lower_face(next_idx, 2)], ...
%                 [lower_face(j, 3), lower_face(next_idx, 3)], 'k');
%       end
%
%       % Plot edges of upper face
%       for j = 1:4
%           next_idx = mod(j, 4) + 1; % Wrap around index for the next vertex
%           plot3([upper_face(j, 1), upper_face(next_idx, 1)], ...
%                 [upper_face(j, 2), upper_face(next_idx, 2)], ...
%                 [upper_face(j, 3), upper_face(next_idx, 3)], 'k');
%       end
%
%       % Plot edges connecting lower and upper faces
%       for j = 1:4
%           plot3([lower_face(j, 1), upper_face(j, 1)], ...
%                 [lower_face(j, 2), upper_face(j, 2)], ...
%                 [lower_face(j, 3), upper_face(j, 3)], 'k');
%       end
%   end
%
%   hold off;
% end
function plotElements3D(elements)
  figure;
  hold on;
  axis equal;
  view(3); % Set 3D view
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  
  % Preallocate arrays for storing line coordinates
  x_lines = [];
  y_lines = [];
  z_lines = [];
  
  for i = 1:size(elements, 1)
      lower_face = squeeze(elements(i, 1:4, :));
      upper_face = squeeze(elements(i, 5:8, :));
      
      % Store edges of lower face
      for j = 1:4
          next_idx = mod(j, 4) + 1; % Wrap around index for the next vertex
          x_lines = [x_lines; lower_face(j, 1), lower_face(next_idx, 1)];
          y_lines = [y_lines; lower_face(j, 2), lower_face(next_idx, 2)];
          z_lines = [z_lines; lower_face(j, 3), lower_face(next_idx, 3)];
      end
      
      % Store edges of upper face
      for j = 1:4
          next_idx = mod(j, 4) + 1; % Wrap around index for the next vertex
          x_lines = [x_lines; upper_face(j, 1), upper_face(next_idx, 1)];
          y_lines = [y_lines; upper_face(j, 2), upper_face(next_idx, 2)];
          z_lines = [z_lines; upper_face(j, 3), upper_face(next_idx, 3)];
      end
      
      % Store edges connecting lower and upper faces
      for j = 1:4
          x_lines = [x_lines; lower_face(j, 1), upper_face(j, 1)];
          y_lines = [y_lines; lower_face(j, 2), upper_face(j, 2)];
          z_lines = [z_lines; lower_face(j, 3), upper_face(j, 3)];
      end
  end
  
  % Plot all lines at once using the line function
  line(x_lines', y_lines', z_lines', 'Color', 'k', 'LineWidth', 1);
  
  hold off;
end




