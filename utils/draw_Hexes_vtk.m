function draw_Hexes_vtk(fname, X,Hexes,BC,str)
t0=tic; iftoiv=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8]; 

if(isinf(X(end,1))==1); X=X(1:end-1,:);end
nX=size(X,1); nH=size(Hexes,1); 

ide=(1:nH)'; 
nF=length(ide);

vtk_title='Hexes'; vtk_title = vtk_title(1:min(length(vtk_title),256));

write_a = @(fid,str) [fprintf(fid,str); fprintf(fid,'\n');];
write_b = @(fid,dat,prec) [fwrite(fid,dat,prec); fprintf(fid,'\n');];
write_b2= @(fid,dat,prec) [fwrite(fid,dat,prec)];

fid=fopen(fname,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,[vtk_title '\n']);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'\n');

fprintf(fid,'POINTS %d float\n',nX);
s='%f %f %f \n'; fprintf(fid,s,X');fprintf(fid,'\n');


fprintf(fid,'CELLS %d %d\n',nH,nH*9);
dat=[8*ones(nH,1),Hexes-1];
fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',dat');
fprintf(fid,'\n');

fprintf(fid,'CELL_TYPES %d\n',nH);
fprintf(fid,'%d\n',12*ones(nH,1));
fprintf(fid,'\n');

tmp = zeros(8,nH);
for i=1:size(BC,1)
for j=iftoiv(BC(i,2), :)
    tmp(j,BC(i,1)) = BC(i,3);
end
end
fprintf(fid,'POINT_DATA %d\n',nF*8);
fprintf(fid,'SCALARS bc_id int 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%d\n',reshape(tmp, 8*nH, []));
fprintf(fid,'\n');
fclose(fid);
