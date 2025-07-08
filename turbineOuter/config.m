% two blades
n_segs = 4;
k_inner = 9;
k_outer = 6;
k_upstream = 10;
k_downstream = 30;
k_blade = 8;

mult = 1.2;
delta_inner = 0.4;
mult_stream = 1.08

Ra = 1250; % Inner square width (from -1250 to 1250)

% Turbine Blade region
R_t = 550; % Bounding rect top
R_b = -550; % Bounding rect bottom (z from R_b to R_t)
R_x = 5400; % End cap position (-R_x and R_x)

% Farfield
R_far_top = 10000
R_far_bottom = -10000
R_far_horizontal = 10000
R_downstream = -10000
R_upstream = 5000


% three blades
n_segs = 20;
k_inner = 48;
k_outer = 36;
k_upstream = 40;
k_downstream = 55;
k_blade = 32;

mult = 1.01;
delta_inner = 0.4;
mult_stream = 1.08

Ra = 12; % Inner square width (from -1250 to 1250)

% Turbine Blade region
R_t = 4.6; % Bounding rect top
R_b = -4.6; % Bounding rect bottom (z from R_b to R_t)
R_x = 65; % End cap position (-R_x and R_x)

% Farfield
R_far_top = 150
R_far_bottom = -90
R_far_horizontal = 150
R_downstream = -200
R_upstream = 80

% three blades rod
%n_segs = 8;
%k_inner = 18;
%k_outer = 12;
%k_upstream = 12;
%k_downstream = 20;
%k_blade = 20;
%
%mult = 1.2;
%delta_inner = 0.4;
%mult_stream = 1.08
%
%Ra = 12; % Inner square width (from -1250 to 1250)
%
%% Turbine Blade region
%R_t = 4.6; % Bounding rect top
%R_b = -78; % Bounding rect bottom (z from R_b to R_t)
%R_b = -78; % Bounding rect bottom (z from R_b to R_t)
%R_x = 73; % End cap position (-R_x and R_x)
%
%% Farfield
%R_far_top = 150
%R_far_bottom = -150
%R_far_horizontal = 150
%R_downstream = -300
%R_upstream = 100
