function exportREA(filename, elements, boundaries)
  fileID = fopen(filename, 'w');

% Add header
  headerID = fopen('header.rea', 'r');
  while ~feof(headerID)
      line = fgets(headerID);
      fprintf(fileID, '%s', line);
  end
  fclose(headerID);

% Mesh
  fprintf(fileID,'**MESH DATA**  1st line is X of corner 1,2,3,4. 2nd line is Y.\n');
  fprintf(fileID,'%8i   %8i   %8i NEL NDIM NEL\n',size(elements,1),2,size(elements,1));
  for i = 1:size(elements,1)
      fprintf(fileID,'            ELEMENT%12i [    1A]    GROUP     0\n',i);
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,1:4,1));
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,1:4,2));
  end

  fprintf(fileID,'  ***** BOUNDARY CONDITIONS ***** \n');
  fprintf(fileID,' ***** FLUID   BOUNDARY CONDITIONS ***** \n');
  for i = 1:size(boundaries,1)
    for j = 1:4
      if boundaries(i,j, 1) == 1
        fprintf(fileID,' W    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 2
        fprintf(fileID,' v    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 3
        fprintf(fileID,' O    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 4
        fprintf(fileID,' SYM  %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 5
        fprintf(fileID,' v    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      else
        fprintf(fileID,' E    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      end
    end
  end

  % fprintf(fileID,'  ***** THERMAL BOUNDARY CONDITIONS *****\n');
  % for i = 1:size(boundaries,1)
  %   for j = 1:4
  %     if boundaries(i,j, 1) == 1
  %       fprintf(fileID,' t    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
  %     elseif boundaries(i,j, 1) == 2
  %       fprintf(fileID,' t    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
  %     elseif boundaries(i,j, 1) == 3
  %       fprintf(fileID,' t    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
  %     else
  %       fprintf(fileID,' E    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
  %     end
  %   end
  % end
% Add footer
  footerID = fopen('footer.rea', 'r');
  while ~feof(footerID)
      line = fgets(footerID);
      fprintf(fileID, '%s', line);
  end
  fclose(footerID);
  fprintf(fileID, '\n');

  fclose(fileID);
end
