function data=readSliceFile(filename);
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Could not open the file. Check the file path and name.');
    end

    fgetl(fileID);
    data = fscanf(fileID, '%f %f %f', [3 Inf]);
    fclose(fileID);

    data = data';
end
