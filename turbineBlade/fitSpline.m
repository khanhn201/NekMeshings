function [pp, arc_length, arc_length_at_max_y]= fitSpline(slice, flipped);
    [min_y, min_idx] = min(slice(:, 2));
    slice = circshift(slice, -min_idx + 1);
    [max_y, max_idx] = max(slice(:, 2));
    slice_c = [slice; slice(1, :)];

    differences = diff(slice_c);
    arc_lengths = sqrt(sum(differences.^2, 2));
    cumulative_arc_lengths = [0; cumsum(arc_lengths)];
    x = cumulative_arc_lengths;
    arc_length = cumulative_arc_lengths(end);
    arc_length_at_max_y = cumulative_arc_lengths(max_idx);

    if flipped == true
        x = x(1:end-1);
        x = circshift(x, -max_idx + 1);
        x(x >= arc_length_at_max_y) = x(x >= arc_length_at_max_y) - arc_length;
        x = [x; arc_length_at_max_y];
        y = circshift(slice, -max_idx + 1); 
        y = [y; y(1, :)];
    else
        y = [slice; slice(1, :)];
    end
    y = [[0;0;0]';y;[0;0;0]'];
    pp = spline(x, y');
    % pp = splinefit(x, y', x, "periodic", true);
end
