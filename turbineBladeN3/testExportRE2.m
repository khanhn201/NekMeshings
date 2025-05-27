elements = [];
boundaries = [];
zs = squeeze([0, 2, 5, 10, 15]);
zs = zs(:);
% xs = [xs; R_end_caps'];
[cylElements, cylBoundaries] = gridBlade(zs);
config
ys = linspace(0, 5, 5)(:);
for k = 2:size(ys,1)
    y_prev = ys(k-1);
    y = ys(k);
    for elem = 1:size(cylElements,1)
        cylElemk = squeeze(cylElements(elem, :, :));
        cylElemk(:, 2) = y_prev;
        cylElemk1 = squeeze(cylElements(elem, :, :));
        cylElemk1(:, 2) = y;
        element = [cylElemk1; cylElemk];

        elements(end+1, :, :) = element;
        [isBoundary, idx] = ismember(elem, cylBoundaries(:, 1));
        if isBoundary
            boundaries(end+1, :) = [size(elements,1); 1; 3];
        end
        if k == 2 && R_downstream == 0
            boundaries(end+1, :) = [size(elements,1); 6; 4];
        end
        if k == size(ys,1)
            boundaries(end+1, :) = [size(elements,1); 5; 2];
        end
    end
end
size(elements)
size(boundaries)
exportRE2('testa', elements, boundaries)
plotBC(elements, boundaries)
