disp("plume")

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

function elements = addSurroundingCylinder(elements, height, r, thickness, refine)
  bottom_faces = disc(r, r + thickness, 2^(refine+2), 'y');
  bottom_faces(:, :, 2) = bottom_faces(:, :, 2) - height/2;
  top_faces = bottom_faces;
  for i = 1:(2^(refine))
    top_faces(:, :, 2) = bottom_faces(:, :, 2) + height/(2^(refine));
    for j = 1:size(bottom_faces, 1)
      element = zeros(1, 8, 3);
      element(1, 1:4, :) = bottom_faces(j, :, :);
      element(1, 5:8, :) = top_faces(j, :, :);
      if checkIfRightHanded(element) == false
        disp('Left-handed element detected SurroundingCylinder')
      endif
      elements = [elements; element];
    end
    bottom_faces = top_faces;
  end
end

function elements = extrudeCylinder(elements, faces, height, thickness)
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

scale = 0.004;
R = 1*scale;
elements = [];
refine = 2;
bc_wall = [];
bc_outflow = [];
bc_inflow = [];

thickness = 0.4*scale;
% Spherical boundary layers
for k = 1: 2 % 1 layers of 0.25 thickness
  elements = addQuadLayer(elements, R, thickness, refine, 'cyl_y', 'cyl_y');
  R = R + thickness;
  if k == 1
    for m = 1:size(elements, 1)
      bc_wall = [bc_wall; m, 5];
    end
  endif
end

% Transition to cylinder
elements = addQuadLayer(elements, R, thickness, refine, 'cyl_y', 'cyl_y');
R = R + thickness;
cyl_height = R;

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

faces = quadSphere(cyl_height, refine, 'cyl_y');
y_extrude_faces = filterYFacingFaces(faces, cyl_height);
external_cyl_faces = [];

% Expand cylinder
thickness = 0.7*scale
for k = 1: 9
  prev_count = size(elements, 1);
  elements = addSurroundingCylinder(elements, 2*cyl_height, R, thickness, refine);
  cyl_faces = disc(R, R + thickness, 2^(refine+2), 'y');
  R = R + thickness;
  thickness = thickness*1.4;
   if k == 9
      for m = (prev_count + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 2];
      end
      external_cyl_faces = cyl_faces;
   else
      y_extrude_faces = cat(1, y_extrude_faces, cyl_faces);
   endif
end

% Extrude cylinder down
cur_height = -cyl_height;
thickness = 1*scale;
for k = 1:3
  prev_count = size(elements, 1);
  elements = extrudeCylinder(elements, y_extrude_faces, cur_height, -thickness);

  prev_count_2 = size(elements, 1);
  elements = extrudeCylinder(elements, external_cyl_faces, cur_height, -thickness);
  for m = (prev_count_2 + 1):size(elements, 1)
     bc_inflow = [bc_inflow; m, 2];
  end

  cur_height = cur_height - thickness;
  thickness = thickness*1.3;
   if k == 3
      for m = (prev_count + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 5];
      end
   endif
end

% Extrude cylinder up
cur_height = cyl_height;
thickness = 0.7*scale;
for k = 1:30
  prev_count = size(elements, 1);
  elements = extrudeCylinder(elements, y_extrude_faces, cur_height, thickness);

  prev_count_2 = size(elements, 1);
  elements = extrudeCylinder(elements, external_cyl_faces, cur_height, thickness);
  for m = (prev_count_2 + 1):size(elements, 1)
     bc_inflow = [bc_inflow; m, 2];
  end

  cur_height = cur_height + thickness;
  thickness = thickness*1.1;
   if k == 30
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
# plotMesh(elements)
plotBC(elements, bc_wall, bc_inflow, bc_outflow);




