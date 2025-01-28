function [pp, arc_length, arc_length_at_max_y]= fitSpline(slice);
    [min_y, min_idx] = min(slice(:, 2));
    slice = circshift(slice, -min_idx + 1);
    slice = [slice; slice(1, :)];
    [max_y, max_idx] = max(slice(:, 2));

    differences = diff(slice);
    arc_lengths = sqrt(sum(differences.^2, 2));
    cumulative_arc_lengths = [0; cumsum(arc_lengths)];
    arc_length = cumulative_arc_lengths(end);
    arc_length_at_max_y = cumulative_arc_lengths(max_idx);

    cumulative_arc_lengths = [-arc_lengths(end); cumulative_arc_lengths; cumulative_arc_lengths(end) + arc_lengths(1)];
    pp = spline(cumulative_arc_lengths, [slice(end-1,:); slice; slice(2,:)]);
end
