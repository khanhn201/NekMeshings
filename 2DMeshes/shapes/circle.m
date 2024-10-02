function lines = circle(r, n, refineYPlus)
    if nargin < 3
        refineYPlus = 0;
    end

    lines = [];
    idx = 0;
 
    for side = 1:4
        refine = n;
        if side == 1
            % Y+ facing
            refine += refineYPlus;
        end
        theta = linspace(0, pi/2, refine+1) + side*pi/2 - pi/4;

        x = r * cos(theta);
        y = r * sin(theta);
        
        for i = 1:refine
            idx += 1;
            lines(idx, 1, 1) = x(i);
            lines(idx, 1, 2) = y(i);
            lines(idx, 2, 1) = x(i+1);
            lines(idx, 2, 2) = y(i+1);
        end
    end
end
