function plotElements3D(elements)
  figure;
  hold on;
  axis equal;
  view(3); % Set 3D view
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  title('QuadSphere Mesh Plot');

  for i = 1:size(elements, 1)
      lower_face = squeeze(elements(i, 1:4, :));
      upper_face = squeeze(elements(i, 5:8, :));
      fill3(lower_face(:, 1), lower_face(:, 2), lower_face(:, 3), 'r', 'FaceAlpha', 0.5);
      fill3(upper_face(:, 1), upper_face(:, 2), upper_face(:, 3), 'b', 'FaceAlpha', 0.5);
      for j = 1:4
          next_idx = mod(j, 4) + 1; % Wrap around index for the next vertex
          side_vertices = [lower_face(j, :); lower_face(next_idx, :); upper_face(next_idx, :); upper_face(j, :)];
          fill3(side_vertices(:, 1), side_vertices(:, 2), side_vertices(:, 3), 'g', 'FaceAlpha', 0.5);
      end
  end

  hold off;
end


