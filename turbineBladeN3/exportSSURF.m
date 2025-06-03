function exportSSURF(filename, surfaces)
  fid = fopen(strcat(filename, '.ssurf'), 'w');
  test = single(6.54321);
  num_dim = 3;
  header = sprintf('#v001%16d%3d hdr', length(surfaces), num_dim);
  hdr = blanks(80);
  hdr(1:length(header)) = header;
  fwrite(fid, hdr, 'char');
  fwrite(fid, test, 'single'); % endian thingy
  
  fwrite(fid, length(surfaces), 'double'); % N_surf
  for i=1:length(surfaces)
      fwrite(fid, surfaces(i).elem, 'double');
      fwrite(fid, surfaces(i).face, 'double');
      for j=1:4
        for dim=1:3
          fwrite(fid, surfaces(i).splineCoeffs(j,dim,:), 'double');
        end
        fwrite(fid, surfaces(i).splineEnds(j,:), 'double');
      end
  end
  fclose(fid);
end
function doub=strToDouble(str)
    bytes = zeros(1, 8, 'uint8');
    for i=1:length(str)
        bytes(i) = str(i);
    end
    doub=typecast(bytes, 'double');
end
