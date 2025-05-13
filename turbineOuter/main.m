[sliceElements,sliceBoundaries] = meshTurbineOuter();
% plotElements(sliceElements, sliceBoundaries);

elements = []
boundaries = []
config
k = 1:k_upstream;
r = (mult_stream.^k-mult_stream)/(mult_stream^k_upstream-mult_stream);
xs = (R_upstream-R_t) * r + R_t;
for k = 2:size(xs,2)
    x_prev = xs(k-1);
    x = xs(k);
    for elem = 1:size(sliceElements,1)
        cylElemk = squeeze(sliceElements(elem, :, :));
        cylElemk(:, 1) = x_prev;
        cylElemk1 = squeeze(sliceElements(elem, :, :));
        cylElemk1(:, 1) = x;
        element = [cylElemk; cylElemk1];

        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx,2);
            if tag ==  1 && k == 2
                boundaries(end+1, :) = [size(elements,1); 5; 5];
            end
            if tag ==  3 || tag == 5|| tag == 9
                boundaries(end+1, :) = [size(elements,1); 2; 2];
            end
            if tag == 7
                boundaries(end+1, :) = [size(elements,1); 3; 1];
            end
        end
        if k == size(xs,2)
            boundaries(end+1, :) = [size(elements,1); 6; 2];
        end
    end
end

xs = linspace(R_b, R_t, k_blade);
for k = 2:size(xs,2)
    x_prev = xs(k-1);
    x = xs(k);
    for elem = 1:size(sliceElements,1)
        cylElemk = squeeze(sliceElements(elem, :, :));
        cylElemk(:, 1) = x_prev;
        cylElemk1 = squeeze(sliceElements(elem, :, :));
        cylElemk1(:, 1) = x;
        element = [cylElemk; cylElemk1];

        checkLeftHanded(element);

        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx,2);
            if tag !=  1
                elements(end+1, :, :) = element;
                if tag ==  3 || tag == 5|| tag == 9
                    boundaries(end+1, :) = [size(elements,1); 2; 2];
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
k = 1:k_downstream;
r = (mult_stream.^k-mult_stream)/(mult_stream^k_downstream-mult_stream);
xs = (R_downstream-R_b) * r + R_b;
for k = 2:size(xs,2)
    x_prev = xs(k-1);
    x = xs(k);
    for elem = 1:size(sliceElements,1)
        cylElemk = squeeze(sliceElements(elem, :, :));
        cylElemk(:, 1) = x_prev;
        cylElemk1 = squeeze(sliceElements(elem, :, :));
        cylElemk1(:, 1) = x;
        element = [cylElemk1; cylElemk];

        checkLeftHanded(element);

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, sliceBoundaries(:, 1));
        if isBoundary
            tag = sliceBoundaries(idx,2);
            if tag ==  1 && k == 2
                boundaries(end+1, :) = [size(elements,1); 6; 5];
            end
            if tag ==  3 || tag == 5|| tag == 9
                boundaries(end+1, :) = [size(elements,1); 2; 2];
            end
            if tag == 7
                boundaries(end+1, :) = [size(elements,1); 3; 1];
            end
        end
        if k == size(xs,2)
            boundaries(end+1, :) = [size(elements,1); 5; 4];
        end
    end
end
size(elements)
size(boundaries)
% plotElements3D(elements)

tmp = elements(:, :, 1);
elements(:, :, 1) = -elements(:, :, 2);
elements(:, :, 2) = tmp;


exportREA("output.rea", elements, boundaries)
% plotBC(elements, boundaries)
