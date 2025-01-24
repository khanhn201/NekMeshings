function smoothMesh()
    % Smooth the layers
    n_smoothing = 10;
    for i = 1:n_smoothing
        for k = 2:k_max - 1
            for t = 1:n_segments
                if t == 1 || t == n_segments / 2
                    continue;
                end
                next_t = mod(t, n_segments) + 1;
                prev_t = t - 1;
                s1 = [all_layers_y{k - 1}(prev_t) - all_layers_y{k - 1}(t); ...
                      all_layers_z{k - 1}(prev_t) - all_layers_z{k - 1}(t)];
                s2 = [all_layers_y{k - 1}(next_t) - all_layers_y{k - 1}(t); ...
                      all_layers_z{k - 1}(next_t) - all_layers_z{k - 1}(t)];
                s3 = [all_layers_y{k}(t) - all_layers_y{k - 1}(t); ...
                      all_layers_z{k}(t) - all_layers_z{k - 1}(t)];

                force_on_s3 = angular_force(s1, s2, s3);
                torque = cross([s3; 0], [force_on_s3; 0]);
                delta_theta = torque(3) * 0.00005;
                rotation_matrix = [cos(delta_theta), -sin(delta_theta); ...
                                   sin(delta_theta), cos(delta_theta)];
                s3_rotated = rotation_matrix * s3;
                s3_new = dot(s3, s3_rotated) / dot(s3_rotated, s3_rotated) * s3_rotated;
                all_layers_y{k}(t) = all_layers_y{k - 1}(t) + s3_new(1);
                all_layers_z{k}(t) = all_layers_z{k - 1}(t) + s3_new(2);
            end
        end
    end
end

function force = angular_force(s1, s2, s3, k_theta)
    if nargin < 4
        k_theta = 1.0;
    end
    s1_norm = norm(s1);
    s2_norm = norm(s2);
    s3_norm = norm(s3);

    cos_theta_1 = dot(s3, s1) / (s3_norm * s1_norm);
    cos_theta_2 = dot(s3, s2) / (s3_norm * s2_norm);

    cos_theta_1 = max(-1.0, min(1.0, cos_theta_1));
    cos_theta_2 = max(-1.0, min(1.0, cos_theta_2));

    theta_1 = acos(cos_theta_1);
    theta_2 = acos(cos_theta_2);

    force_spring_1 = k_theta * theta_1 * (s1 / s1_norm);
    force_spring_2 = k_theta * theta_2 * (s2 / s2_norm);

    force = force_spring_1 + force_spring_2;
end
