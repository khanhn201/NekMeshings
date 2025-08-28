function [elements, boundaries] = importRE2(filename)
  fid = fopen(filename, 'r');
  if fid < 0
    error('Could not open file %s.re2', filename);
  end

  % --- Read header (80 chars)
  hdr = fread(fid, 80, 'char=>char')';
  nums = sscanf(hdr, '%*s %d %d %d %d %*s');
  num_elem = nums(1)
  num_dim = nums(2)

  % --- endian test
  test = fread(fid, 1, 'single'); %#ok<NASGU> % usually 6.54321

  % --- read elements
  elements = zeros(num_elem, 8, num_dim);
  for i = 1:num_elem
    igroup = fread(fid, 1, 'double'); %#ok<NASGU> % not used
    x = fread(fid, 8, 'double');
    y = fread(fid, 8, 'double');
    z = fread(fid, 8, 'double');
    elements(i,:,1) = x;
    elements(i,:,2) = y;
    elements(i,:,3) = z;
  end

  % --- read curve (ignored here, always zero)
  curveFlag = fread(fid, 1, 'double'); %#ok<NASGU>

  % --- read boundary conditions
  nBCread = fread(fid, 1, 'double');
  nBCread = int32(nBCread);
  boundaries = zeros(nBCread, 3);

  bcMap = getBCMap();
  for i = 1:nBCread
    buf2 = fread(fid, 8, 'double');
    elem = int32(buf2(1));
    face = int32(buf2(2));
    tagD = buf2(8);

    tagStr = sprintf('%-3s', strtrim(doubleToStr(tagD)));

    % default -1 if unknown
    if bcMap.str2id.isKey(tagStr)
        tag = bcMap.str2id(tagStr);
    else
        tagStr
        tag = 0;
    end

    boundaries(i,:) = [elem, face, tag];
  end

  fclose(fid);
end

function str = doubleToStr(doub)
  bytes = typecast(doub, 'uint8');
  % strip trailing nulls and turn into chars
  str = char(bytes(bytes > 0));
end
