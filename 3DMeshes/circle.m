function faces = circle(r, height, n)
    center = [0, height, 0];
    faces = zeros(n, 4, 3);
    theta_step = 2 * pi / n;
    for i = 1:n
        theta1 = (i - 1) * theta_step;
        theta3 = (i) * theta_step;
        theta2 = (i - 1) * theta_step + theta_step/2;
        
        p1 = [r * cos(theta1), height, r * sin(theta1)];
        p2 = [r * cos(theta2), height, r * sin(theta2)];
        p3 = [r * cos(theta3), height, r * sin(theta3)];
        faces(i, :, :) = [center; p3; p2; p1];
    end
end

