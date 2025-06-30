function exportRE2(filename, elements, boundaries)
  % Rea file
  fileID = fopen(strcat(filename, '.rea'), 'w');

% Add header
  headerID = fopen('header.rea', 'r');
  while ~feof(headerID)
      line = fgets(headerID);
      fprintf(fileID, '%s', line);
  end
  fclose(headerID);

% Mesh
  fprintf(fileID,'**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n');
  fprintf(fileID,'      %8i   %8i   %8i NELT,NDIM,NELV\n',-size(elements,1),3,size(elements,1));
% Add footer
  footerID = fopen('footer.rea', 'r');
  while ~feof(footerID)
      line = fgets(footerID);
      fprintf(fileID, '%s', line);
  end
  fclose(footerID);
  fprintf(fileID, '\n');
  fclose(fileID);

  % Re2 file
  fid = fopen(strcat(filename, '.re2'), 'w');
  test = single(6.54321);
  num_dim = 3;
  num_elem = size(elements, 1);
  nBC = 1;
  header = sprintf('#v004%16d%3d%16d%4d hdr', num_elem, num_dim, num_elem, nBC);
  hdr = blanks(80);
  hdr(1:length(header)) = header;
  fwrite(fid, hdr, 'char');
  fwrite(fid, test, 'single'); % endian thingy

  % xyz
  for i = 1:size(elements,1)
    fwrite(fid, 0, 'double'); % igroup
    fwrite(fid, elements(i, 1:8, 1), 'double'); %x
    fwrite(fid, elements(i, 1:8, 2), 'double'); %y
    fwrite(fid, elements(i, 1:8, 3), 'double'); %z
  end
  % curve
  fwrite(fid, 0, 'double'); % 0 curve
  % boundaries
  fwrite(fid, size(boundaries, 1), 'double'); % bc

  for i = 1:size(boundaries, 1)
    elem = boundaries(i, 1);
    face = boundaries(i, 2);
    tag = boundaries(i, 3);
    buf2 = zeros(1, 8);
    buf2(1) = double(elem);
    buf2(2) = double(face);
    buf2(3) = 0.0;
    buf2(4) = 0.0;
    buf2(5) = 0.0;
    buf2(6) = 0.0;
    buf2(7) = 0.0;
    if tag == 1
      buf2(8) = strToDouble("W  ");
    elseif tag == 2
      buf2(8) = strToDouble("v  ");
    elseif tag == 3
      buf2(8) = strToDouble("int");
    elseif tag == 4
      buf2(8) = strToDouble("O  ");
    elseif tag == 5
      buf2(8) = strToDouble("s  ");
    elseif tag == 6
      buf2(8) = strToDouble("SYM");
    else
      tag
    end
    fwrite(fid, buf2, 'double'); %x
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
