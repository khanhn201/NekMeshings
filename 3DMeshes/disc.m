function faces = disc(r1, r2, n, axis)
    theta = linspace(0, 2*pi, n+1);
    faces = zeros(n, 4, 3);

    innerRing = [r1 * cos(theta)', r1 * sin(theta)', zeros(n+1, 1)];
    outerRing = [r2 * cos(theta)', r2 * sin(theta)', zeros(n+1, 1)];

    switch axis
        case 'x'
            innerRing = innerRing(:, [3, 1, 2]);
            outerRing = outerRing(:, [3, 1, 2]);
        case 'y'
            innerRing = innerRing(:, [2, 3, 1]);
            outerRing = outerRing(:, [2, 3, 1]);
    end

    for i = 1:n
        faces(i, :, :) = [innerRing(i, :);
                          outerRing(i, :);
                          outerRing(i+1, :);
                          innerRing(i+1, :)];
    end
end
