function TriMeshBiHDMSVD( name, svdname )

[A, B] = TriMeshTKBinaryRead(strcat(name,'.bihdmat'));
S = TriMeshTKNystromSVD(A, B);

dlmwrite(strcat(svdname,'.bihdmat.svd'), S);
end

