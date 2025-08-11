% function slice=splineSlices(slicesCoord, skip=1);
%     nSlice = size(slicesCoord)(1)
%     for i = 1:skip:nSlice
%         slicesCoord(i, :, :)
%         break
%     end
% end


function [pp, arc_length, arc_length_at_max_y]= fitSpline(slice);
    [max_y, max_idx] = max(slice(:, 1));
    slice_c = [slice; slice(1, :)];

    differences = diff(slice_c);
    arc_lengths = sqrt(sum(differences.^2, 2));
    cumulative_arc_lengths = [0; cumsum(arc_lengths)];
    x = cumulative_arc_lengths;
    arc_length = cumulative_arc_lengths(end);
    arc_length_at_max_y = cumulative_arc_lengths(max_idx);

    y = [slice; slice(1, :)];
    % y = [[0;0;0]';y;[0;0;0]'];
    y = [
        (y(2,:)-y(1,:))/(x(2)-x(1)); 
        y; 
        (y(end,:)-y(end-1,:))/(x(end)-x(end-1));];
    pp = spline(x, y');
    % pp = splinefit(x, y', x, "periodic", true);
end
