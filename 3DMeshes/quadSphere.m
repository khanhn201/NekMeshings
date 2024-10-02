function faces = quadSphere(r, n, shape)
    faces = [];
    cube_faces = cube(r);
    if isscalar(n)
        n = repmat(n, 1, 6);
    else
        error('n must be either a scalar or a 6-element array.');
    end
    for f = 1:6
        current_cube_faces = cube_faces(f, :, :);
        for i = 1:n(f)
          new_faces = [];
          for j = 1:size(current_cube_faces,1)
            face = current_cube_faces(j, :, :);
            face = reshape(face, 4, 3);

            subfaces = subdivideAndProject(r, face, shape);
            new_faces = cat(1, new_faces, subfaces);
          end
          current_cube_faces = new_faces;
        end
        faces = cat(1, faces, current_cube_faces);
    end
end

function faces = subdivideAndProject(r, face, shape)
    midpoints = zeros(4, 3);
    for i = 1:4
        face(i, :) = projectInto(r, shape, face(i, :));
    end
    midpoints(1, :) = (face(1, :) + face(2, :)) / 2;
    midpoints(2, :) = (face(2, :) + face(3, :)) / 2;
    midpoints(3, :) = (face(3, :) + face(4, :)) / 2;
    midpoints(4, :) = (face(4, :) + face(1, :)) / 2;

    center = mean(face, 1);

    for i = 1:4
        midpoints(i, :) = adjustMidpointIfOnRadius(r, shape, face(i, :), face(mod(i, 4) + 1, :), midpoints(i, :));
    end
    for i = 1:4
        midpoints(i, :) = projectInto(r, shape, midpoints(i, :));
    end
    center = projectInto(r, shape, center);

    faces = zeros(4, 4, 3);
    faces(1, :, :) = [face(1, :); midpoints(1, :); center; midpoints(4, :)];
    faces(2, :, :) = [midpoints(1, :); face(2, :); midpoints(2, :); center];
    faces(3, :, :) = [center; midpoints(2, :); face(3, :); midpoints(3, :)];
    faces(4, :, :) = [midpoints(4, :); center; midpoints(3, :); face(4, :)];
end

function midpoint = adjustMidpointIfOnRadius(r, shape, point1, point2, midpoint)
    if isPointOnRadius(r, shape, point1) && isPointOnRadius(r, shape, point2)
        switch shape
            case 'cyl_x'
                normVal = norm(midpoint(2:3));
                if normVal ~= 0
                    midpoint(2:3) = midpoint(2:3) * r / normVal;
                end
            case 'cyl_y'
                normVal = norm([midpoint(1), midpoint(3)]);
                if normVal ~= 0
                    midpoint([1, 3]) = midpoint([1, 3]) * r / normVal;
                end
            case 'cyl_z'
                normVal = norm(midpoint(1:2));
                if normVal ~= 0
                    midpoint(1:2) = midpoint(1:2) * r / normVal;
                end
        end
    end
end

function onRadius = isPointOnRadius(r, shape, point)
    switch shape
        case 'cyl_x'
            onRadius = abs(norm(point(2:3)) - r) < 1e-6;
        case 'cyl_y'
            onRadius = abs(norm([point(1), point(3)]) - r) < 1e-6;
        case 'cyl_z'
            onRadius = abs(norm(point(1:2)) - r) < 1e-6;
        otherwise
            onRadius = false;
    end
end

function point = projectInto(r, shape, point)
        switch shape
            case 'sphere'
                point = point.*r / norm(point);
            case 'cube'
                point = point;
            case 'cyl_x'
                normVal = norm(point(2:3));
                if (abs(point(1)) != r) || (abs(point(2)) == r) || (abs(point(3)) == r)
                    if normVal
                        point(2:3) = point(2:3) * r / normVal;
                    endif
                endif
            case 'cyl_y'
                normVal = norm([point(1), point(3)]);
                if (abs(point(2)) != r) || (abs(point(1)) == r) || (abs(point(3)) == r)
                    if normVal
                        point([1, 3]) = point([1, 3]) * r / normVal;
                    endif
                endif
            case 'cyl_z'
                normVal = norm(point(1:2));
                if (abs(point(3)) != r) || (abs(point(1)) == r) || (abs(point(2)) == r)
                    if normVal
                        point(1:2) = point(1:2) * r / normVal;
                    endif
                endif
            otherwise
                error('Invalid type');
        end
end

function faces = cube(r)
  faces = zeros(6, 4, 3);
  vertices = [
    1,  1,  1;
    -1,  1,  1;
    -1, -1,  1;
    1, -1,  1;
    1,  1, -1;
    -1,  1, -1;
    -1, -1, -1;
    1, -1, -1
  ].*r;

  faceIndices = [
    1, 2, 3, 4;
    1, 4, 8, 5;
    1, 5, 6, 2;
    2, 6, 7, 3;
    4, 3, 7, 8;
    5, 8, 7, 6
  ];

  for i = 1:6
      faces(i,:,:) = vertices(faceIndices(i,:), :);
  end
end
