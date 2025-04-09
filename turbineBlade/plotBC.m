

function plotBC(elements, boundaries)
  face_node_map = [
    1 2 6 5; % Face 1
    2 3 7 6; % Face 2
    3 4 8 7; % Face 3
    1 4 8 5; % Face 4
    1 2 3 4; % Face 5
    5 6 7 8 % Face 6
  ];
  figure;
  view(3); % 3D view
  hold on;
  color_matrix = [
      0.8, 0.2, 0.2;  % Tag 1: Wall (Red)
      0.2, 0.8, 0.2;  % Tag 2: Inflow (Green)
      0.2, 0.2, 0.8;  % Tag 3: Outflow (Blue)
      0.8, 0.2, 0.8   % Tag 4: Internal (Purple)
      0.5, 0.5, 0.5   % Tag 4: Internal (Purple)
  ];
  for elem = [45561]
      for face = 1:6
          coordinates = elements(elem, face_node_map(face, :), :);
          x = coordinates(:,:, 1);
          y = coordinates(:,:, 2);
          z = coordinates(:,:, 3);
          color = color_matrix(5, :);
          fill3(x, y, z, color, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
      end
  end
  % Plot bc_wall
  for i = 1:size(boundaries, 1)
      elem = boundaries(i, 1);
      face = boundaries(i, 2);
      tag = boundaries(i, 3);
      coordinates = elements(elem, face_node_map(face, :), :);
      if tag == 4
          x = coordinates(:,:, 1);
          y = coordinates(:,:, 2);
          z = coordinates(:,:, 3);
          color = color_matrix(tag, :);
          fill3(x, y, z, color, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
      end
  end
  for i = 1:size(boundaries, 1)
      elem = boundaries(i, 1);
      face = boundaries(i, 2);
      tag = boundaries(i, 3);
      if tag == 1
          coordinates = elements(elem, face_node_map(face, :), :);
          x = coordinates(:,:, 1);
          y = coordinates(:,:, 2);
          z = coordinates(:,:, 3);
          color = color_matrix(tag, :);
          fill3(x, y, z, color, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
      end
  end
  for i = 1:size(boundaries, 1)
      elem = boundaries(i, 1);
      face = boundaries(i, 2);
      tag = boundaries(i, 3);
      coordinates = elements(elem, face_node_map(face, :), :);
      if tag == 2
          x = coordinates(:,:, 1);
          y = coordinates(:,:, 2);
          z = coordinates(:,:, 3);
          color = color_matrix(tag, :);
          fill3(x, y, z, color, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
      end
  end
  for i = 1:size(boundaries, 1)
      elem = boundaries(i, 1);
      face = boundaries(i, 2);
      tag = boundaries(i, 3);
      coordinates = elements(elem, face_node_map(face, :), :);
      if tag == 3
          x = coordinates(:,:, 1);
          y = coordinates(:,:, 2);
          z = coordinates(:,:, 3);
          color = color_matrix(tag, :);
          fill3(x, y, z, color, 'FaceAlpha', 0.5, 'EdgeColor', 'black'); % Plot faces with transparency
      end
  end
  % Add labels and legend
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  axis equal;
  grid on;

  hold off;
end
