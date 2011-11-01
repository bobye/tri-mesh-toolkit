function [L D] = TriMeshTKFBaseRead( matrix_filename, eigenvalue_filename )
D = load(eigenvalue_filename);
n = length(D);
fdata = fopen(matrix_filename,'r');
L = fread(fdata, [n-1, n-1], 'double');
end

