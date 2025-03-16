n_segs = 4;
k_inner = 9;
k_outer = 6;
k_upstream = 10;
k_downstream = 30;
k_blade = 6;

mult = 1.2;
delta_inner = 0.4;
mult_stream = 1.05

Ra = 1250; % Inner square width (from -1250 to 1250)

% Turbine Blade region
R_t = 550; % Bounding rect top
R_b = -550; % Bounding rect bottom (z from R_b to R_t)
R_x = 5400; % End cap position (-R_x and R_x)

% Farfield
R_far_top = 10000
R_far_bottom = -10000
R_far_horizontal = 10000
R_downstream = -15000
R_upstream = 5000


