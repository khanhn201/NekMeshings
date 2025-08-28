function slicesCoord=readSlices(filename);
    load(filename);
    slice_coords = lofted_shape;
    tmp = slice_coords(:, :, 1);
    slice_coords(:, :, 1) = -slice_coords(:, :, 2);
    slice_coords(:, :, 2) = -tmp;
    slicesCoord = slice_coords(:, 1:end-1, :); % opened airfoil
end
