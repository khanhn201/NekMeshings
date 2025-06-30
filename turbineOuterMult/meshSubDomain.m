function [elements,boundaries] = meshSubDomain(front, back, left, right, shiftStream, shiftPer, 
    upstream_layers, downstream_layers,
    boundariesMap % Front: v, Back: O, Left: SYM, Right: SYM 
)
[sliceElements,sliceBoundaries] = meshTurbineOuterVar(left, right, shiftPer);
% plotElements(sliceElements, sliceBoundaries);

front = front - shiftStream;
back = back - shiftStream;

elements = []
boundaries = []
config
k = 1:upstream_layers;
r = (mult_stream.^k-mult_stream)/(mult_stream^upstream_layers-mult_stream);
xs = (front-R_t) * r + R_t;
for k = 2:size(xs,2)
    x_prev = xs(k-1);
    x = xs(k)
    slice = sliceElements(:, :, :);
    slicek1 = sliceElements(:, :, :);
    slice(:, :, 1) = x_prev;
    slicek1(:, :, 1) = x;
    for elem = 1:size(sliceElements,1)
        cylElemk = squeeze(slice(elem, :, :));
        cylElemk1 = squeeze(slicek1(elem, :, :));
        element = [cylElemk; cylElemk1];

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx,2);
            if tag ==  1 && k == 2
                boundaries(end+1, :) = [size(elements,1); 5; 5];
            end
            if tag == 3 % Top
                boundaries(end+1, :) = [size(elements,1); 2; 6];
            end
            if tag == 9 && boundariesMap(3) != 0 % Left
                boundaries(end+1, :) = [size(elements,1); 2; boundariesMap(3)];
            end
            if tag == 5 && boundariesMap(4) != 0 % Right
                boundaries(end+1, :) = [size(elements,1); 2; boundariesMap(4)];
            end
            if tag == 7
                boundaries(end+1, :) = [size(elements,1); 3; 1];
            end
        end
        if k == size(xs,2)
            boundaries(end+1, :) = [size(elements,1); 6; boundariesMap(1)];
        end
    end
end

xs = linspace(R_b, R_t, k_blade);
for k = 2:size(xs,2)
    x_prev = xs(k-1);
    x = xs(k)
    slice = sliceElements(:, :, :);
    slicek1 = sliceElements(:, :, :);
    slice(:, :, 1) = x_prev;
    slicek1(:, :, 1) = x;
    for elem = 1:size(sliceElements,1)
        cylElemk = squeeze(slice(elem, :, :));
        cylElemk1 = squeeze(slicek1(elem, :, :));
        element = [cylElemk; cylElemk1];

        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx,2);
            if tag !=  1
                elements(end+1, :, :) = element;
                if tag == 3 % Top
                    boundaries(end+1, :) = [size(elements,1); 2; 6];
                end
                if tag == 9 && boundariesMap(3) != 0 % Left
                    boundaries(end+1, :) = [size(elements,1); 2; boundariesMap(3)];
                end
                if tag == 5 && boundariesMap(4) != 0 % Right
                    boundaries(end+1, :) = [size(elements,1); 2; boundariesMap(4)];
                end
                if tag == 7
                    boundaries(end+1, :) = [size(elements,1); 3; 1];
                end
                if tag ==  2 || tag == 4|| tag == 8
                    boundaries(end+1, :) = [size(elements,1); 4; 3];
                end
                if tag == 6
                    boundaries(end+1, :) = [size(elements,1); 1; 3];
                end
            end
        else
            elements(end+1, :, :) = element;
        end
    end
end
% r = -(R_b - R_downstream) / R_t;
% xs = R_b * (r .^ linspace(0, 1, k_downstream))
k = 1:downstream_layers;
r = (mult_stream.^k-mult_stream)/(mult_stream^downstream_layers-mult_stream);
xs = (back-R_b) * r + R_b;
for k = 2:size(xs,2)
    x_prev = xs(k-1);
    x = xs(k)
    slice = sliceElements(:, :, :);
    slicek1 = sliceElements(:, :, :);
    slice(:, :, 1) = x_prev;
    slicek1(:, :, 1) = x;
    for elem = 1:size(sliceElements,1)
        cylElemk = squeeze(slice(elem, :, :));
        cylElemk1 = squeeze(slicek1(elem, :, :));
        element = [cylElemk1; cylElemk];


        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx,2);
            if tag ==  1 && k == 2
                boundaries(end+1, :) = [size(elements,1); 6; 5];
            end
            if tag == 3 % Top
                boundaries(end+1, :) = [size(elements,1); 2; 6];
            end
            if tag == 9 && boundariesMap(3) != 0 % Left
                boundaries(end+1, :) = [size(elements,1); 2; boundariesMap(3)];
            end
            if tag == 5 && boundariesMap(4) != 0 % Right
                boundaries(end+1, :) = [size(elements,1); 2; boundariesMap(4)];
            end
            if tag == 7
                boundaries(end+1, :) = [size(elements,1); 3; 1];
            end
        end
        if k == size(xs,2)
            boundaries(end+1, :) = [size(elements,1); 5; boundariesMap(2)];
        end
    end
end
for k = 1:size(elements,1)
    checkLeftHanded(squeeze(elements(k, :, :)));
end
% plotElements3D(elements)

elements(:, :, 1) += shiftStream;
tmp = elements(:, :, 1);
elements(:, :, 1) = -elements(:, :, 2);
elements(:, :, 2) = tmp;


% exportREA("turbineOuter.rea", elements, boundaries)
% plotBC(elements, boundaries)
