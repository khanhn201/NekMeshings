function lines = square(R, n, refineYPlus)
    if nargin < 3
        refineYPlus = 0; % Default value
    end


    lines = [];
    idx = 0;

    x_corners = [R, -R, -R, R];
    y_corners = [R, R, -R, -R];
    for side = 1:4
        x_start = x_corners(side);
        y_start = y_corners(side);
        x_end = x_corners(mod(side, 4) + 1);
        y_end = y_corners(mod(side, 4) + 1);

        refine = n;
        if side == 1
            % Y+ facing
            refine += refineYPlus;
        end
        x = linspace(x_start, x_end, refine + 1);
        y = linspace(y_start, y_end, refine + 1);

        for i = 1:refine
            idx += 1;
            lines(idx, 1, 1) = x(i);
            lines(idx, 1, 2) = y(i);
            lines(idx, 2, 1) = x(i+1);
            lines(idx, 2, 2) = y(i+1);
        end
    end
    
end
