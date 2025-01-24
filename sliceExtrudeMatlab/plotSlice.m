function plotSlice(slice);
    x = slice(:, 1);
    y = slice(:, 2);
    z = slice(:, 3);
    figure;
    scatter3(x, y, z, 'filled');
    title('3D Scatter Plot of Coordinates');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    grid on;
end
