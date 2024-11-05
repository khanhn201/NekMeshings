disp("cyl")

function elements = extrudeLayer(elements, faces, height, thickness)
  lower_faces = faces;
  lower_faces(:, :, 2) = lower_faces(:, :, 2) + height;
  upper_faces = lower_faces;
  upper_faces(:, :, 2) = lower_faces(:, :, 2) + thickness;
  for i = 1:size(lower_faces, 1)
    element = zeros(1, 8, 3);
    if (thickness < 0)
      element(1, 5:8, :) = lower_faces(i, :, :);
      element(1, 1:4, :) = upper_faces(i, :, :);
    else
      element(1, 1:4, :) = lower_faces(i, :, :);
      element(1, 5:8, :) = upper_faces(i, :, :);
    endif
    if checkIfRightHanded(element) == false
      disp('Left-handed element detected')
    endif
    elements = [elements; element];
  end
end

function filtered_faces = filterYFacingFaces(faces, pos)
    filtered_faces = [];
    for i = 1:size(faces, 1)
        face = faces(i, :, :);
        face_y_coords = face(:, :, 2);
        if all(face_y_coords == pos)
            face(:, :, 2) = 0;
            filtered_faces = [filtered_faces; face];
        end
    end
end
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

function graphFaces(faces)
  figure;
  hold on;
  view(3); % Set 3D view
  for i = 1:size(faces, 1)
    face = squeeze(faces(i, :, :));
    patch(face(:, 1), face(:, 2), face(:, 3), 'g');
  endfor
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  axis equal;
  grid on;
  hold off;
end

R = 0.5;
thickness = 0.15;
thickness_y = 0.5;
mult = 1.5;
mult_y = 1.05;
layer_count = 40;
surrounding_count = 13;
refine = 2;

elements = [];
bc_wall = [];
bc_outflow = [];
bc_inflow = [];

center_faces = quadSphere(R, refine, 'cyl_y');
center_faces = filterYFacingFaces(center_faces, R);

external_faces = [];
surrounding_faces = [];
first_surrounding_faces = [];

for k = 1:surrounding_count
  cyl_faces = disc(R, R + thickness, 2^(refine+2), 'y');
  R = R + thickness;
  thickness = thickness*mult;
  if k == surrounding_count
      external_faces = cyl_faces;
  elseif k == 1
      first_surrounding_faces = cyl_faces
  else
      surrounding_faces = cat(1, surrounding_faces, cyl_faces);
  end
end

% graphFaces(center_faces);
% graphFaces(surrounding_faces);
% graphFaces(external_faces);

cur_height = 0.0;
elements = extrudeLayer(elements, center_faces, cur_height, thickness_y);
for m = 1:size(elements, 1)
    bc_wall = [bc_wall; m, 5];
end
prev_count = size(elements, 1);
elements = extrudeLayer(elements, surrounding_faces, cur_height, thickness_y);
elements = extrudeLayer(elements, first_surrounding_faces, cur_height, thickness_y);
prev_count = size(elements, 1);
elements = extrudeLayer(elements, external_faces, cur_height, thickness_y);
for m = (prev_count + 1):size(elements, 1)
    bc_inflow = [bc_inflow; m, 2];
end
cur_height = cur_height + thickness_y;

% plotElements(elements)

% Extrude cylinder up
for k = 1:layer_count
    prev_count = size(elements, 1);
    elements = extrudeLayer(elements, center_faces, cur_height, thickness_y);
    elements = extrudeLayer(elements, surrounding_faces, cur_height, thickness_y);
    elements = extrudeLayer(elements, first_surrounding_faces, cur_height, thickness_y);

    prev_count_2 = size(elements, 1);
    elements = extrudeLayer(elements, external_faces, cur_height, thickness_y);
    for m = (prev_count_2 + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 2];
    end
    cur_height = cur_height + thickness_y;
    thickness_y = thickness_y*mult_y;
    if k == layer_count
        for m = (prev_count + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 6];
        end
    endif
end

cur_height
% Elevate
cur_height = 0.0;
thickness_y = -2.0;
layer_count = 3;
for k = 1:layer_count
    prev_count = size(elements, 1);
    elements = extrudeLayer(elements, first_surrounding_faces, cur_height, thickness_y);
    for m = (prev_count + 1):size(elements, 1)
        bc_wall = [bc_wall; m, 4];
    end
    elements = extrudeLayer(elements, surrounding_faces, cur_height, thickness_y);

    prev_count_2 = size(elements, 1);
    elements = extrudeLayer(elements, external_faces, cur_height, thickness_y);
    for m = (prev_count_2 + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 2];
    end
    cur_height = cur_height + thickness_y;
    if k == layer_count
        for m = (prev_count + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 5];
        end
    endif
end




% Postprocess
R
cur_height

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
% plotMesh(elements)
plotBC(elements, bc_wall, bc_inflow, bc_outflow);




