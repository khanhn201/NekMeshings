% main.m

addpath('shapes');

function elements = addLayer(elements, r, thickness, refine)
  inner_lines = circle(r, refine);
  outer_lines = circle(r + thickness, refine);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
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

      elements = [elements; element];

      element = zeros(1, 4, 2);
      element(1, 1, :) = intermediate_lines(1, 2, :);
      element(1, 2, :) = intermediate_lines(1, 1, :);
      element(1, 3, :) = square_lines(i, 1, :);
      element(1, 4, :) = square_lines(i, 2, :);

      elements = [elements; element];
  end
end

function elements = addElementSide(elements, N, R, R_max, height, thickness, multiplier)
  inner_lines = geometricLine(N, R, R_max, multiplier, height);
  outer_lines = geometricLine(N, R, R_max, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    if R < 0
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
    elements = [elements; element];
  end
end

elements = [];
bc_wall = [];
bc_outflow = [];
bc_inflow = [];

R = 1;
refine = 8;

N = size(elements, 1);
boundaries = zeros(N, 6, 5);

% Circle wall
elements = addLayer(elements, R, 0.25, refine);
R += 0.25;
for m = 1:size(elements, 1)
  bc_wall = [bc_wall; m, 1];
end
% Expand circle
for i = 1:1
  elements = addLayer(elements, R, 0.25, refine);
  R += 0.25;
end
% Circle -> Square
elements = addTransitionalLayer(elements, R, 0.5, refine, 0.5);
R += 0.5;
top_layer = square(R, refine);
top_layer = top_layer(1:floor(length(top_layer)/4), :, :);

% Expand sideway
R_max = 5;
sideway_count = 12;
multiplier = 0.9;
for i = 1:refine
  elements = addElementSide(elements, sideway_count, R, R_max, -R + (i-1)*R/refine*2, R/refine*2, multiplier);
  bc_wall = [bc_wall; size(elements, 1), 2];
  elements = addElementSide(elements, sideway_count, -R, -R_max, -R + (i-1)*R/refine*2, R/refine*2, multiplier);
  bc_wall = [bc_wall; size(elements, 1), 2];
end

% Extrude upward
height = R;
thickness = 1;
extrude_mult = 1.3;
extrude_count = 10;
for k = 1:extrude_count
  prev_count = size(elements, 1);
  inner_lines = top_layer;
  inner_lines(:, :, 2) = height;
  outer_lines = top_layer;
  outer_lines(:, :, 2) = height + thickness;
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    elements = [elements; element];
  end
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc_outflow = [bc_outflow; m, 3];
    end
  end

  inner_lines = geometricLine(sideway_count, R, R_max, multiplier, height);
  outer_lines = geometricLine(sideway_count, R, R_max, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 1, :);
    element(1, 2, :) = inner_lines(i, 2, :);
    element(1, 3, :) = outer_lines(i, 2, :);
    element(1, 4, :) = outer_lines(i, 1, :);
    elements = [elements; element];
  end
  bc_wall = [bc_wall; size(elements, 1), 2];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc_outflow = [bc_outflow; m, 3];
    end
  end


  inner_lines = geometricLine(sideway_count, -R, -R_max, multiplier, height);
  outer_lines = geometricLine(sideway_count, -R, -R_max, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    elements = [elements; element];
  end
  bc_wall = [bc_wall; size(elements, 1), 4];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc_outflow = [bc_outflow; m, 3];
    end
  end

  height += thickness;
  thickness = thickness*extrude_mult;
end


% Extrude downward
height = -R;
thickness = -1;
extrude_mult = 1.2;
extrude_count = 4;
for k = 1:extrude_count
  prev_count = size(elements, 1);
  inner_lines = top_layer;
  inner_lines(:, :, 2) = height;
  outer_lines = top_layer;
  outer_lines(:, :, 2) = height + thickness;
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 1, :);
    element(1, 2, :) = inner_lines(i, 2, :);
    element(1, 3, :) = outer_lines(i, 2, :);
    element(1, 4, :) = outer_lines(i, 1, :);
    elements = [elements; element];
  end
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc_inflow = [bc_inflow; m, 3];
    end
  end

  inner_lines = geometricLine(sideway_count, R, R_max, multiplier, height);
  outer_lines = geometricLine(sideway_count, R, R_max, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 2, :);
    element(1, 2, :) = inner_lines(i, 1, :);
    element(1, 3, :) = outer_lines(i, 1, :);
    element(1, 4, :) = outer_lines(i, 2, :);
    elements = [elements; element];
  end
  bc_wall = [bc_wall; size(elements, 1), 4];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc_inflow = [bc_inflow; m, 3];
    end
  end


  inner_lines = geometricLine(sideway_count, -R, -R_max, multiplier, height);
  outer_lines = geometricLine(sideway_count, -R, -R_max, multiplier, height + thickness);
  for i = 1:size(inner_lines, 1)
    element = zeros(1, 4, 2);
    element(1, 1, :) = inner_lines(i, 1, :);
    element(1, 2, :) = inner_lines(i, 2, :);
    element(1, 3, :) = outer_lines(i, 2, :);
    element(1, 4, :) = outer_lines(i, 1, :);
    elements = [elements; element];
  end
  bc_wall = [bc_wall; size(elements, 1), 2];
  if k == extrude_count
    for m = prev_count+1:size(elements, 1)
      bc_inflow = [bc_inflow; m, 3];
    end
  end

  height += thickness;
  thickness = thickness*extrude_mult;
end


N = size(elements, 1)

% BC
for i = 1:size(bc_wall, 1)
  boundaries(bc_wall(i,1), bc_wall(i,2), 1) = 1;
end
for i = 1:size(bc_inflow, 1)
  boundaries(bc_inflow(i,1), bc_inflow(i,2), 1) = 2;
end
for i = 1:size(bc_outflow, 1)
  boundaries(bc_outflow(i,1), bc_outflow(i,2), 1) = 3;
end

exportREA("output.rea", elements, boundaries)

% plotMesh(elements, bc_wall, bc_inflow, bc_outflow, true);
