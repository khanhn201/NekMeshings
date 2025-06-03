function exportToVTK(filename, hexMesh)
    % hexMesh is (N, 8, 3)
    N = size(hexMesh, 1);

    % Open file
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file for writing.');
    end

    % Write header
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Hex mesh exported from MATLAB\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n');

    % Write points
    numPoints = size(uniquePoints, 1);
    fprintf(fid, 'POINTS %d float\n', numPoints);
    fprintf(fid, '%f %f %f\n', uniquePoints');

    % Write cells
    totalIndices = N * 9; % 8 points per hex + 1 for count
    fprintf(fid, '\nCELLS %d %d\n', N, totalIndices);
    for i = 1:N
        fprintf(fid, '8 %d %d %d %d %d %d %d %d\n', connectivity(i, :));
    end

    % Write cell types (VTK_HEXAHEDRON = 12)
    fprintf(fid, '\nCELL_TYPES %d\n', N);
    fprintf(fid, '%d\n', repmat(12, N, 1));

    fclose(fid);
    fprintf('VTK file "%s" written successfully.\n', filename);
end
