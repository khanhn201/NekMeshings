% three blades
n_segs = 10;
k_inner = 12;
k_outer = 10;
k_upstream = 25;
k_downstream = 35;
k_blade = 8;
k_spacing = 10;

mult = 1.2;
delta_inner = 0.4;
mult_stream = 1.08

Ra = 12; % Inner square width

% Turbine Blade region
R_t = 4.6; % Bounding rect top
R_b = -4.6; % Bounding rect bottom (z from R_b to R_t)
R_x = 65; % End cap position (-R_x and R_x)

second_offset = 20;

% Farfield
R_far_top = 200
R_far_bottom = -90 % Rated hub height
R_far_horizontal = 300
R_downstream = -200
R_upstream = 100
