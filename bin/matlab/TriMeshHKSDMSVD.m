function TriMeshBiHDMSVD( name, svdname )

[A, B] = TriMeshTKBinaryRead(strcat(name,'.hksdmat'));
S = TriMeshTKNystromSVD(A, B);

dlmwrite(strcat(svdname,'.hksdmat.svd'), S);
end

