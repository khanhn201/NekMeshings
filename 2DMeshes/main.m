% main.m

addpath('shapes');


function isCCW = isCounterClockwise(element)
    x = element(:, 1);
    y = element(:, 2);
    area = 0.5 * (x(1)*y(2) + x(2)*y(3) + x(3)*y(4) + x(4)*y(1) - ...
                  (y(1)*x(2) + y(2)*x(3) + y(3)*x(4) + y(4)*x(1)));
    isCCW = area > 0;
end

function elements = addLayer(elements, r, thickness, refine)
  inner_lines = circle(r, refine);
  outer_lines = circle(r + thickness, refine);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    if ~isCounterClockwise(squeeze(element))
        fprintf("not counter clockwise layer\n")
    end
    elements = [elements; element];
  end
end

function elements = addTransitionalLayer(elements, r, thickness, refine, alpha)
  circle_lines = circle(r, refine);
  square_lines = square(r + thickness, refine);
  for i = 1:size(circle_lines, 1)
      intermediate_lines = alpha * square_lines(i, :, :) + (1 - alpha) * circle_lines(i, :, :);

      element = zeros(1, 4, 2);
      element(1, 1, :) = circle_lines(i, 2, :);
      element(1, 2, :) = circle_lines(i, 1, :);
      element(1, 3, :) = intermediate_lines(1, 1, :);
      element(1, 4, :) = intermediate_lines(1, 2, :);
      if ~isCounterClockwise(squeeze(element))
        fprintf("not counter clockwise transit\n")
      end

      elements = [elements; element];

      element = zeros(1, 4, 2);
      element(1, 1, :) = intermediate_lines(1, 2, :);
      element(1, 2, :) = intermediate_lines(1, 1, :);
      element(1, 3, :) = square_lines(i, 1, :);
      element(1, 4, :) = square_lines(i, 2, :);
      if ~isCounterClockwise(squeeze(element))
        fprintf("not counter clockwise transit\n")
      end

      elements = [elements; element];
  end
end

function elements = addElementSide(elements, N, R, R_max, height, thickness, multiplier)
  inner_lines = geometricLine(N, R, R_max, multiplier, height);
  outer_lines = geometricLine(N, R, R_max, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    if R > 0
        element(1, 1, :) = outer_lines(i, 1, :);
        element(1, 2, :) = outer_lines(i, 2, :);
        element(1, 3, :) = inner_lines(i, 2, :);
        element(1, 4, :) = inner_lines(i, 1, :);
    else
        element(1, 1, :) = inner_lines(i, 1, :);
        element(1, 2, :) = inner_lines(i, 2, :);
        element(1, 3, :) = outer_lines(i, 2, :);
        element(1, 4, :) = outer_lines(i, 1, :);
    end
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise side\n")
    end
    elements = [elements; element];
  end
end

elements = [];
bc = []; % Element Face Type; 1: Wall, 2: Inflow, 3: Outflow, 4: Sym, 5: Moving Wall

R = 0.5;
refine = 8;


% Circle wall
elements = addLayer(elements, R, 0.025, refine);
R += 0.025;
for m = 1:size(elements, 1)
  bc = [bc; m, 1, 1];
end
% Expand circle
% for i = 1:1
%   elements = addLayer(elements, R, 0.025, refine);
%   R += 0.025;
% end
% Circle -> Square
elements = addTransitionalLayer(elements, R, 0.025, refine, 0.5);
R += 0.025
top_layer = square(R, refine);
top_layer = top_layer(floor(length(top_layer)/4)+1:2*floor(length(top_layer)/4), :, :);

% Expand upward
R_top = 20;
count_top = 15;
multiplier = 1.2;
for i = 1:refine
  elements = addElementSide(elements, count_top, R, R_top, -R + (i-1)*R/refine*2, R/refine*2, multiplier);
  bc = [bc; size(elements, 1), 2, 4];
end

R_bottom = -0.6;
count_bottom = 2;
% Expand downward
for i = 1:refine
  elements = addElementSide(elements, count_bottom, -R, R_bottom, -R + (i-1)*R/refine*2, R/refine*2, multiplier);
  bc = [bc; size(elements, 1), 2, 5];
end


% Extrude right
height = R;
thickness = 0.25;
extrude_mult = 1.2;
extrude_count = 25;
for k = 1:extrude_count
  prev_count = size(elements, 1);
  inner_lines = top_layer;
  inner_lines(:, :, 1) = height;
  outer_lines = top_layer;
  outer_lines(:, :, 1) = height + thickness;
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 1, :);
    element(1, 2, :) = inner_lines(i, 2, :);
    element(1, 3, :) = outer_lines(i, 2, :);
    element(1, 4, :) = outer_lines(i, 1, :);
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise right mid\n")
    end
    elements = [elements; element];
  end
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc = [bc; m, 3, 3];
    end
  end

  inner_lines = geometricLine(count_top, R, R_top, multiplier, height);
  outer_lines = geometricLine(count_top, R, R_top, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise right top\n")
    end
    elements = [elements; element];
  end
  bc = [bc; size(elements, 1), 4, 4];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc = [bc; m, 3, 3];
    end
  end


  inner_lines = geometricLine(count_bottom, -R, R_bottom, multiplier, height);
  outer_lines = geometricLine(count_bottom, -R, R_bottom, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 1, :);
    element(1, 2, :) = inner_lines(i, 2, :);
    element(1, 3, :) = outer_lines(i, 2, :);
    element(1, 4, :) = outer_lines(i, 1, :);
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise right bottom\n")
    end
    elements = [elements; element];
  end
  bc = [bc; size(elements, 1), 2, 5];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc = [bc; m, 3, 3];
    end
  end

  height += thickness;
  thickness = thickness*extrude_mult;
end


% Extrude left
height = -R;
thickness = -0.25;
extrude_mult = 1.2;
extrude_count = 20;
for k = 1:extrude_count
  prev_count = size(elements, 1);
  inner_lines = top_layer;
  inner_lines(:, :, 1) = height;
  outer_lines = top_layer;
  outer_lines(:, :, 1) = height + thickness;
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise left mid\n")
    end
    elements = [elements; element];
  end
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc = [bc; m, 3, 2];
    end
  end

  inner_lines = geometricLine(count_top, R, R_top, multiplier, height);
  outer_lines = geometricLine(count_top, R, R_top, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 1, :);
    element(1, 2, :) = inner_lines(i, 2, :);
    element(1, 3, :) = outer_lines(i, 2, :);
    element(1, 4, :) = outer_lines(i, 1, :);
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise left top\n")
    end
    elements = [elements; element];
  end
  bc = [bc; size(elements, 1), 2, 4];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc = [bc; m, 3, 2];
    end
  end


  inner_lines = geometricLine(count_bottom, -R, R_bottom, multiplier, height);
  outer_lines = geometricLine(count_bottom, -R, R_bottom, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    if ~isCounterClockwise(squeeze(element))
      fprintf("not counter clockwise left bottom\n")
    end
    elements = [elements; element];
  end
  bc = [bc; size(elements, 1), 4, 5];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc = [bc; m, 3, 2];
    end
  end

  height += thickness;
  thickness = thickness*extrude_mult;
end


N = size(elements, 1)
boundaries = zeros(N, 6, 5);
AR = 0.4;
elements(:, :, 1) = elements(:, :, 1)* AR;
% BC
for i = 1:size(bc, 1)
  boundaries(bc(i,1), bc(i,2), 1) = bc(i,3);
end

histID = fopen('output.his', 'w');
count = 30;
fprintf(histID, '%3d\n', count);
theta = linspace(0, 2*pi, count+1);
for i = 1:count
  x = 0.5 * cos(theta(i))*AR;
  y = 0.5 * sin(theta(i));
  z = 0.0;
  fprintf(histID, '%13.5e %13.5e %13.5e\n', x, y, z);
end
fclose(histID);
exportREA("output.rea", elements, boundaries)

plotMesh(elements, bc, true);
