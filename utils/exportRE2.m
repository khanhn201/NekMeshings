function exportRE2(filename, elements, boundaries, num_dim=3, nScalar=0)

  % Re2 file
  fid = fopen(strcat(filename, '.re2'), 'w');
  test = single(6.54321);
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
    for j = 1:num_dim
      fwrite(fid, elements(i, 1:2^num_dim, j), 'double');
    end
  end
  % curve
  fwrite(fid, 0, 'double'); % 0 curve
  fwrite(fid, size(boundaries, 1), 'double'); % bc


  % boundaries fluid
  bcMap = getBCMap();
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
    if bcMap.id2str.isKey(tag)
        buf2(8) = strToDouble(bcMap.id2str(tag));
    else
        tag
        buf2(8) = 0; % unknown
    end

    fwrite(fid, buf2, 'double'); %x
  end
  % boundaries scalar
  for i = 1:nScalar
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
        if bcMap.id2strS.isKey(tag)
            buf2(8) = strToDouble(bcMap.id2strS(tag));
        else
            tag
            buf2(8) = 0; % unknown
        end

        fwrite(fid, buf2, 'double'); %x
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
