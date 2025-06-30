function [elements,boundaries] = meshTurbineOuterDouble(front, back, shift)
    config
    elements = [];
    boundaries = [];
    [elements,boundaries] = meshSubDomain(front, back, left, right, shift)
    elements(:,:, 2) -= shift;
    elementsLeft = elements;
    elementsLeft(:, :, 1) += R_a;
    elementsRight = elements;
    elementsRight(:, :, 1) -= R_a;
end
