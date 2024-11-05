function plotMesh(elements, bc, plotElements, showNumber)
  if nargin < 3
    plotElements = false;
  end
  if nargin < 4
    showNumber = false;
  end

  face_node_map = [
    1 2; % Face 1
    2 3; % Face 2
    3 4; % Face 3
    1 4; % Face 4
  ];
  figure;
  hold on;

  % Plot elements
  if plotElements
    for i = 1:size(elements, 1)
      element = squeeze(elements(i, :, :));
      x_points = [element(:, 1); element(1, 1)];
      y_points = [element(:, 2); element(1, 2)];
      
      plot(x_points, y_points, '-black');
      if showNumber
        centroid_x = mean(element(:, 1));
        centroid_y = mean(element(:, 2));
        text(centroid_x, centroid_y, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
      end
    end
  end



  % Plot bc
  colors = jet(max(bc(:, 3)));
  for i = 1:size(bc, 1)
    elem = bc(i, 1);
    face = bc(i, 2);
    color = colors(bc(i, 3), :);

    coordinates = elements(elem, face_node_map(face, :), :);
    x = coordinates(:,:, 1);
    y = coordinates(:,:, 2);
    plot(x, y, '-', 'Color', color, 'LineWidth', 1.5);
  end

  axis equal;
  grid on;

  hold off;
end
