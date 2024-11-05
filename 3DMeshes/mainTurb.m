disp("cyl")

function elements = extrudeLayer(elements, faces, height, thickness)
  if size(faces, 1) == 0
      return
  endif
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

R = 0.5;
thickness = 1;
thickness_y = 1;
mult = 1;
layer_count = 2;
blade_count = 2;
blade_y_count = 2;
surrounding_count = 2;
refine = 3;
cur_height = 0;

elements = [];
bc_wall = [];
bc_outflow = [];
bc_inflow = [];
bc_int = [];

center_faces = circle(R, 0, 2^(refine+1));

external_faces = [];

surrounding_faces = cell(1, 8);
surrounding_faces_inner = cell(1, 8);

cyl_faces_inner = cell(1, 8);
cyl_faces_arr = cell(1, 8);
cyl_faces_left = cell(1, 8);
cyl_faces_right = cell(1, 8);
cyl_faces_outer = cell(1, 8);
for kk = 1:8
    surrounding_faces{kk} = [];
    surrounding_faces_inner{kk} = [];
    cyl_faces_inner{kk} = [];
    cyl_faces_arr{kk} = [];
    cyl_faces_left{kk} = [];
    cyl_faces_right{kk} = [];
    cyl_faces_outer{kk} = [];
end

for k = 1:blade_count
  cyl_faces = disc(R, R + thickness, 2^(refine+2), 'y');
    for kk = 1:8
      start_idx = (kk - 1) * 2^(refine-1) + 1;
      end_idx = kk * 2^(refine-1);
      if k == 1
        cyl_faces_inner{kk} = cat(1, cyl_faces_inner{kk}, cyl_faces(start_idx:end_idx, :, :));
      elseif k == blade_count
        cyl_faces_outer{kk} = cat(1, cyl_faces_outer{kk}, cyl_faces(start_idx:end_idx, :, :));
      else
        cyl_faces_arr{kk} = cat(1, cyl_faces_arr{kk}, cyl_faces(start_idx+1:end_idx-1, :, :));
        cyl_faces_left{kk} = cat(1, cyl_faces_left{kk}, cyl_faces(start_idx, :, :));
        cyl_faces_right{kk} = cat(1, cyl_faces_right{kk}, cyl_faces(end_idx, :, :));
      end
    end
  R = R + thickness;
  thickness = thickness*mult;
end

mult = 1.5;
for k = 1:surrounding_count
  cyl_faces = disc(R, R + thickness, 2^(refine+2), 'y');

  R = R + thickness;
  thickness = thickness*mult;
  for kk = 1:8
      start_idx = (kk - 1) * 2^(refine-1) + 1;
      end_idx = kk * 2^(refine-1);

      if k == surrounding_count
          external_faces = cyl_faces;
      elseif k == 1
          surrounding_faces_inner{kk} = cat(1, surrounding_faces_inner{kk}, cyl_faces(start_idx:end_idx, :, :));
      else
          surrounding_faces{kk} = cat(1, surrounding_faces{kk}, cyl_faces(start_idx:end_idx, :, :));
      end
  end
end


% Begin extruding
for k = 1:blade_y_count
    for kk = 1:8
        if mod(kk, 2) == 0
            continue
        end
        elements = extrudeLayer(elements, cyl_faces_arr{kk}, cur_height, thickness_y);
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, cyl_faces_inner{kk}, cur_height, thickness_y);
        for m = (prev_count + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 4];
        end
        bc_wall = [bc_wall; (prev_count + 1), 1];
        bc_wall = [bc_wall; size(elements,1), 3];
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, cyl_faces_right{kk}, cur_height, thickness_y);
        for m = (prev_count + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 3];
        end
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, cyl_faces_left{kk}, cur_height, thickness_y);
        for m = (prev_count + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 1];
        end
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, cyl_faces_outer{kk}, cur_height, thickness_y);
        % for m = (prev_count + 1):size(elements, 1)
        %     bc_wall = [bc_wall; m, 2];
        % end
        bc_wall = [bc_wall; (prev_count + 1), 1];
        bc_wall = [bc_wall; size(elements,1), 3];
    end
    for kk = 1:8
        elements = extrudeLayer(elements, surrounding_faces{kk}, cur_height, thickness_y);
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, surrounding_faces_inner{kk}, cur_height, thickness_y);
        if mod(kk, 2) == 0
            for m = (prev_count + 1):size(elements, 1)
                bc_wall = [bc_wall; m, 4];
            end
        end
    end
    prev_count = size(elements, 1);
    elements = extrudeLayer(elements, external_faces, cur_height, thickness_y);
    for m = (prev_count + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 2];
    end
    cur_height = cur_height + thickness_y;
end

thickness_y = 2;
for k = 1:layer_count
    prev_count_t = size(elements, 1);
    elements = extrudeLayer(elements, center_faces, cur_height, thickness_y);
    if k == 1
        for m = (prev_count_t + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 5];
        end
    end
    for kk = 1:8
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, cyl_faces_arr{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_inner{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_outer{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_left{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_right{kk}, cur_height, thickness_y);
        if k == 1 && mod(kk,2)==0
            for m = (prev_count + 1):size(elements, 1)
                bc_wall = [bc_wall; m, 5];
            end
        end
        elements = extrudeLayer(elements, surrounding_faces{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, surrounding_faces_inner{kk}, cur_height, thickness_y);
    end

    prev_count_2 = size(elements, 1);
    elements = extrudeLayer(elements, external_faces, cur_height, thickness_y);
    for m = (prev_count_2 + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 2];
    end
    cur_height = cur_height + thickness_y;
    if k == layer_count
        for m = (prev_count_t + 1):size(elements, 1)
            bc_inflow = [bc_inflow; m, 6];
        end
    endif
end

cur_height = 0
thickness_y = -thickness_y
for k = 1:layer_count
    prev_count_t = size(elements, 1);
    elements = extrudeLayer(elements, center_faces, cur_height, thickness_y);
    if k == 1
        for m = (prev_count_t + 1):size(elements, 1)
            bc_wall = [bc_wall; m, 6];
        end
    end
    for kk = 1:8
        prev_count = size(elements, 1);
        elements = extrudeLayer(elements, cyl_faces_arr{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_inner{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_outer{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_left{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, cyl_faces_right{kk}, cur_height, thickness_y);
        if k == 1 && mod(kk,2)==0
            for m = (prev_count + 1):size(elements, 1)
                bc_wall = [bc_wall; m, 6];
            end
        end
        elements = extrudeLayer(elements, surrounding_faces{kk}, cur_height, thickness_y);
        elements = extrudeLayer(elements, surrounding_faces_inner{kk}, cur_height, thickness_y);
    end

    prev_count_2 = size(elements, 1);
    elements = extrudeLayer(elements, external_faces, cur_height, thickness_y);
    for m = (prev_count_2 + 1):size(elements, 1)
        bc_inflow = [bc_inflow; m, 2];
    end
    cur_height = cur_height + thickness_y;
    if k == layer_count
        for m = (prev_count_t + 1):size(elements, 1)
            bc_inflow = [bc_inflow; m, 5];
        end
    endif
end



bc_int = bc_inflow;
bc_inflow = [];

% Postprocess

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
% Int boundary: 4
for i = 1:size(bc_int, 1)
  boundaries(bc_int(i,1), bc_int(i,2), 1) = 4;
end



exportREA("output.rea", elements, boundaries)
% plotElements(elements)
plotBC(elements, bc_wall, bc_inflow, bc_outflow, bc_int);




