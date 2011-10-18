function [mat_A, mat_B] = TriMeshTKBinaryRead(filename)
% read binary Nystrom sampling matrix

% file IO
finfo = fopen(strcat(filename,'.info'),'r');
N = fscanf(finfo, '#size %d', 1);
sampleidx = fscanf(finfo, '%d', [1, inf])' +1;
unsampleidx = 1:N;
unsampleidx(sampleidx) =[];
fclose(finfo);
fdata = fopen(filename,'r');
mat_origin = fread(fdata, [N, length(sampleidx)], 'double');
fclose(fdata);
% permutation
mat_A = mat_origin(sampleidx, :);
mat_B = mat_origin(unsampleidx, :);
end