function TriMeshBiHDMSVD( name )

[A, B] = TriMeshTKBinaryRead(strcat(name,'.bihdmat'));
S = TriMeshTKNystromSVD(A, B);

dlmwrite(strcat(name,'.bihdmat.svd'), S);
end

