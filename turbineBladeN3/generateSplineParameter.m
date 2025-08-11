function s_fine=generateSplineParameter(arc_length_at_max_y, arc_length, n);
    % Cheb
    i = 0:(n/2);
    cheb_nodes = cos(i*pi / (n/2));
    cheb_nodes_mapped = (cheb_nodes + 1) / 2;
    s_fine = cheb_nodes_mapped * arc_length_at_max_y;
    s_fine1 = flip(s_fine)(1:end-1);
    s_fine = cheb_nodes_mapped * (arc_length-arc_length_at_max_y) + arc_length_at_max_y;
    s_fine2 = flip(s_fine)(1:end-1);
    s_fine = [s_fine1(1), s_fine1(3:end-1),...
              s_fine2(1), s_fine2(3:end-1)];
    s_fine = [s_fine1, s_fine2];

    % Linear
    s_fine = [linspace(0, arc_length_at_max_y, n/2 + 1)(1:end-1)',
              linspace(arc_length_at_max_y, arc_length, n/2 + 1)(1:end-1)']';

end
