% Read the CSV file
filename = 'turbineBladeWall';
data = csvread(filename, 1, 1); % Skip header row, start from second column

% Extract coordinates
x = data(:, [1, 4, 7, 10]); % Columns x1, x2, x3, x4
y = data(:, [2, 5, 8, 11]); % Columns y1, y2, y3, y4
z = data(:, [3, 6, 9, 12]); % Columns z1, z2, z3, z4

% Plot quadrilaterals
figure;
hold on;
for i = 1:size(x, 1)
    quad_x = [x(i, :) x(i, 1)];
    quad_y = [y(i, :) y(i, 1)];
    quad_z = [z(i, :) z(i, 1)];
    
    plot3(quad_x, quad_y, quad_z, '-k');
end
size(x, 1)

title('3D Quadrilaterals');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
hold off;
