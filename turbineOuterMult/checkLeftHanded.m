function checkLeftHanded(element)
    corner_list = [
        1, 2, 4, 5;
        2, 3, 1, 6;
        3, 4, 2, 7;
        4, 1, 3, 8;
        5, 8, 6, 1;
        6, 5, 7, 2;
        7, 6, 8, 3;
        8, 7, 5, 4
    ];

    for i = 1:size(corner_list, 1)
        face = corner_list(i, :);
        node1 = element(face(1), :);
        node2 = element(face(2), :);
        node4 = element(face(3), :);
        node5 = element(face(4), :);
        
        vec12 = node2 - node1;
        vec14 = node4 - node1;
        vec15 = node5 - node1;

        cross_vec = cross(vec12, vec14);
        dot_prod = dot(cross_vec, vec15);

        if dot_prod <= 0
            disp('Left-handed/ill-shaped element detected!');
            return;
        end
    end
end

