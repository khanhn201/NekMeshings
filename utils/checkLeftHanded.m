function checkLeftHanded(elements)
  Ne = size(elements,1);
  for e = 1:Ne
    
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
        node1 = elements(e,face(1), :);
        node2 = elements(e,face(2), :);
        node4 = elements(e,face(3), :);
        node5 = elements(e,face(4), :);
        
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
end

