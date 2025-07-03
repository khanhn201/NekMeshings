config
 [elements,boundaries] = meshSubDomain(R_upstream, -12, -R_far_horizontal, R_far_horizontal, 0, 0, 
        k_upstream, 7,
        boundariesMap=[2, 3, 6, 6]
);
exportRE2("outer1", elements, boundaries);

N = size(elements,1);
X = permute(elements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk('outer1.vtk', X,Hexes, boundaries,'');
size(elements)


 [elements1,boundaries1] = meshSubDomain(-8, R_downstream, 0, R_far_horizontal, -20, 80, 
        7, k_downstream,
        boundariesMap=[3, 4, 0, 6]
);
 [elements2,boundaries2] = meshSubDomain(-8, R_downstream, -R_far_horizontal, 0, -20, -80, 
        7, k_downstream,
        boundariesMap=[3, 4, 6, 0]
);
elements = [elements1; elements2];
boundaries2(:, 1) += size(elements1, 1);
boundaries = [boundaries1; boundaries2];



exportRE2("outer2", elements, boundaries);

N = size(elements,1);
X = permute(elements, [2, 1, 3]);
X = reshape(X, [], 3); 
Hexes = reshape(1:(N*8), 8, N)';
draw_Hexes_vtk('outer2.vtk', X,Hexes, boundaries,'');
size(elements)

