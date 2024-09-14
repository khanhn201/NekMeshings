function exportREA(filename, elements, boundaries)
  fileID = fopen(filename, 'w');
  fprintf(fileID,'%8i   %8i   %8i NEL NDIM NEL\n',size(elements,1),3,size(elements,1));
  for i = 1:size(elements,1)
      fprintf(fileID,'            ELEMENT%12i [    1A]    GROUP     0\n',i);
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,1:4,1));
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,1:4,2));
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,1:4,3));
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,5:8,1));
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,5:8,2));
      fprintf(fileID,'%15.7g  %15.7g  %15.7g  %15.7g\n',elements(i,5:8,3));
  end

  fprintf(fileID,' ***** FLUID   BOUNDARY CONDITIONS ***** \n');
  for i = 1:size(boundaries,1)
    for j = 1:6
      if boundaries(i,j, 1) == 1
        fprintf(fileID,' W    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 2
        fprintf(fileID,' v    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 3
        fprintf(fileID,' O    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      else
        fprintf(fileID,' E    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      end
    end
  end

  fprintf(fileID,'  ***** THERMAL BOUNDARY CONDITIONS *****\n');
  for i = 1:size(boundaries,1)
    for j = 1:6
      if boundaries(i,j, 1) == 1
        fprintf(fileID,' t    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 2
        fprintf(fileID,' t    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      elseif boundaries(i,j, 1) == 3
        fprintf(fileID,' t    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      else
        fprintf(fileID,' E    %3d %13.5e %13.5e %13.5e %13.5e %13.5e\n', i, 0, 0, 0, 0, 0);
      end
    end
  end
  fclose(fileID);
end
