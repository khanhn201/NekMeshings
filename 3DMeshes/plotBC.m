

function plotBC(elements, bc_wall, bc_inflow, bc_outflow)
  face_node_map = [
    1 2 5 6; % Face 1
    2 3 7 6; % Face 2
    3 4 7 8; % Face 3
    1 4 5 8; % Face 4
    1 2 3 4; % Face 5
    5 6 7 8 % Face 6
  ];
  figure;
  hold on;
  color_wall = [0.8, 0.2, 0.2];    % Red for walls
  color_inflow = [0.2, 0.8, 0.2];  % Green for inflow
  color_outflow = [0.2, 0.2, 0.8]; % Blue for outflow

  % Plot bc_wall
  for i = 1:size(bc_wall, 1)
      elem = bc_wall(i, 1);
      face = bc_wall(i, 2);
      coordinates = elements(elem, face_node_map(face, :), :);
      x = coordinates(:,:, 1);
      y = coordinates(:,:, 2);
      z = coordinates(:,:, 3);
      patch(x, y, z, color_wall, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
  end

  % Plot bc_inflow
  for i = 1:size(bc_inflow, 1)
      elem = bc_inflow(i, 1);
      face = bc_inflow(i, 2);
      coordinates = elements(elem, face_node_map(face, :), :);
      x = coordinates(:,:, 1);
      y = coordinates(:,:, 2);
      z = coordinates(:,:, 3);
      patch(x, y, z, color_inflow, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
  end

  % Plot bc_outflow
  for i = 1:size(bc_outflow, 1)
      elem = bc_outflow(i, 1);
      face = bc_outflow(i, 2);
      coordinates = elements(elem, face_node_map(face, :), :);
      x = coordinates(:,:, 1);
      y = coordinates(:,:, 2);
      z = coordinates(:,:, 3);
      patch(x, y, z, color_outflow, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
  end

  % Add labels and legend
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title('Boundary Conditions Visualization');
  view(3); % 3D view
  axis equal;
  grid on;

  hold off;
end
