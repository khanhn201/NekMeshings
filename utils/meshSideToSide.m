function [elements, boundaries] = meshSideToSide(bottom, top, p)
  Nx = size(p, 2);
  Ny = size(top, 1);
  layers = zeros(Nx, Ny, 3);
  for i = 1:Nx
      layers(i,:,:) = (1-p(i))*bottom + p(i)*top;
  end

  Ne = (Ny-1)*(Nx-1);
  elements = zeros(Ne, 4, 3);
  boundaries = [];
  e = 1;
  for i = 1:Nx-1
      for j = 1:Ny-1
          elements(e,1,:) = squeeze(layers(i,   j, :))';
          elements(e,2,:) = squeeze(layers(i,   j+1,   :))';
          elements(e,3,:) = squeeze(layers(i+1, j+1,   :))';
          elements(e,4,:) = squeeze(layers(i+1, j, :))';
          if i == 1
            boundaries(end+1, :) = [e; 1; 1];
          end
          if i == Nx-1
            boundaries(end+1, :) = [e; 3; 3];
          end
          if j == 1
            boundaries(end+1, :) = [e; 2; 2];
          end
          if j == Ny-1
            boundaries(end+1, :) = [e; 4; 4];
          end
          e = e + 1;
      end
  end
end
