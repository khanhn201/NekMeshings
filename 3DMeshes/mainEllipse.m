disp("ellipse")

function isRightHanded = checkIfRightHanded(element)
  v1 = element(1, 1, :);  % Vertex 1
  v2 = element(1, 2, :);  % Vertex 2
  v4 = element(1, 4, :);  % Vertex 4
  v5 = element(1, 5, :);  % Vertex 5

  vec12 = v2 - v1;
  vec14 = v4 - v1;
  vec15 = v5 - v1;
  cross_vec = cross(vec12, vec14);
  dot_prod = dot(cross_vec, vec15);

  isRightHanded = dot_prod > 0;
end

function elements = addQuadLayer(elements, r, thickness, refine, type1, type2)
  lower_faces = quadSphere(r, refine, type1);
  upper_faces = quadSphere(r + thickness, refine, type2);
  for i = 1:size(lower_faces, 1)
    element = zeros(1, 8, 3);
    element(1, 1:4, :) = lower_faces(i, :, :);
    element(1, 5:8, :) = upper_faces(i, :, :);
    if checkIfRightHanded(element) == false
      disp('Left-handed element detected')
    endif
    elements = [elements; element];
  end
end

function thicknesses = thicknessesGeometric(R_start, R_end, N, multiplier)
    thicknesses = zeros(N, 1);
    total_length = R_end - R_start;
    ratio = (multiplier^N - 1) / (multiplier - 1);
    thicknesses(1) = total_length / ratio;
    for i = 2:N
        thicknesses(i) = thicknesses(i-1) * multiplier;
    end
end

scale = 0.004;
refine = 4;
elements = [];
bc_wall = [];
bc_outflow = [];
bc_inflow = [];

R_start = 1*scale;
R_end = 50*scale;
layer_count = 10;
multiplier = 1.2;
thicknesses = thicknessesGeometric(R_start, R_end, layer_count, multiplier);

R = R_start;
for k = 1:layer_count
  thickness = thicknesses(k);
  prev_count = size(elements, 1);
  elements = addQuadLayer(elements, R, thickness, refine, 'sphere', 'sphere');
  R = R + thickness;
  if k == 1
    for m = (prev_count + 1):size(elements, 1)
      bc_wall = [bc_wall; m, 5];
    end
  endif
  if k == layer_count
      for m = (prev_count + 1):size(elements, 1)
        bc_outflow = [bc_outflow; m, 6];
      end
  endif
end



N = size(elements, 1)
boundaries = zeros(N, 6, 5);

% Wall boundary: 1
for i = 1:size(bc_wall, 1)
  boundaries(bc_wall(i,1), bc_wall(i,2), 1) = 1;
end
% Inflow boundary: 2
for i = 1:size(bc_inflow, 1)
  boundaries(bc_inflow(i,1), bc_inflow(i,2), 1) = 2;
end
% Outflow boundary: 3
for i = 1:size(bc_outflow, 1)
  boundaries(bc_outflow(i,1), bc_outflow(i,2), 1) = 3;
end



exportREA("output.rea", elements, boundaries)
% plotElements(elements)
plotBC(elements, bc_wall, bc_inflow, bc_outflow);
