function draw_Hexes_vtk(X,Hexes,BC,str)
t0=tic; iftoiv=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8]; 

fname = 'Hexes.vtk';

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
% fid=fopen(fname,'wb','ieee-be');
%
% write_a(fid,'# vtk DataFile Version 2.0');
%
% write_a(fid,vtk_title);
% write_a(fid,'BINARY\n');
% write_a(fid,'DATASET UNSTRUCTURED_GRID');
%
% write_a(fid,['POINTS ' num2str(nX) ' float']);
% write_b(fid,X','float32');
% write_a(fid,'');
%
% write_a(fid,['CELLS ' num2str(nH) ' ' num2str(9*nH)]);
% dat=uint32([8*ones(nH,1),Hexes-1]);
% write_b2(fid,dat','uint32');
% write_a(fid,''); write_a(fid,'');
%
% write_a(fid,['CELL_TYPES ' num2str(nH)]);
% write_b(fid,uint32(12*ones(1,nH)),'uint32');
% write_a(fid,'');
%
% write_a(fid,['POINT_DATA ' num2str(nH*8)]);
% write_a(fid,'SCALARS bc int 1');
% write_a(fid,'LOOKUP_TABLE default');
% write_b(fid,uint32(reshape(tmp,[1, nH*8])),'uint32');
% write_a(fid,'');
%
% write_a(fid,['CELL_DATA ' num2str(nH)]);
% write_a(fid,'SCALARS cell_id int 1');
% write_a(fid,'LOOKUP_TABLE default');
% % tmp = 0*ones(1,nH);
% % for i=1:size(BC,1)
% %     tmp(BC(i,1)) = BC(i,3);
% % end
% write_b(fid,uint32(1:nH),'uint32');
% write_a(fid,'');
%
%
% % write_a(fid,['CELLS ' num2str(nF) ' ' num2str(5*nF)]);
% % for i=1:size(BC,1);ie=BC(i, 1);iF=BC(i,2); 
% %     dat=uint32([4,Hexes(ie,iftoiv(iF,:))-1]);
% %     write_b2(fid,dat,'uint32'); % no break line
% % end;
% % write_a(fid,''); write_a(fid,'');
% %
% % write_a(fid,['CELL_TYPES ' num2str(nF)]);
% % write_b(fid,uint32(7*ones(1,nF)),'uint32');
% % write_a(fid,'');
% %
% % write_a(fid,['CELL_DATA ' num2str(nF)]);
% % write_a(fid,'SCALARS cell_id int 1');
% % write_a(fid,'LOOKUP_TABLE default');
% % in1=sub2ind(size(CBC),ide,idf);
% % write_b(fid,uint32(CBC(in1)),'uint32');
% % write_a(fid,'');
% %
fclose(fid);
% %
