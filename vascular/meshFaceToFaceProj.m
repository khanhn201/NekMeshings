function [elements, boundaries] = meshFaceToFaceProj(bottom, top, p, sliceBoundaries, boundaryMap)
  Nz = length(p);
  Nf = size(top, 1);
  layers = zeros(Nz, Nf, 4, 3);
  for i = 1:Nz
      layers(i,:,:,:) = (1-p(i))*bottom + p(i)*top;
  end

  Ne = Nf*(Nz-1);
  elements = zeros(Ne, 8, 3);
  boundaries = [];
  e = 1;
  for i = 1:Nz-1
      for j = 1:Nf
          [isBoundary, idx] = ismember(j, sliceBoundaries(:, 1));
          if isBoundary && ~isempty(boundaryMap)
              tag = sliceBoundaries(idx,3);
              face = sliceBoundaries(idx,2);
              [found, mapIdx] = ismember(tag, boundaryMap(:,1));
              if found
                newtag = boundaryMap(mapIdx, 2);
                boundaries(end+1, :) = [e; face; newtag];
              end
          end
          elements(e,1:4,:) = layers(i, j, :, :);
          elements(e,5:8,:) = layers(i+1, j, :, :);
          if i == 1
            boundaries(end+1, :) = [e; 5; 5];
          end
          if i == Nz-1
            boundaries(end+1, :) = [e; 6; 6];
          end


          e = e + 1;
      end
  end
end
