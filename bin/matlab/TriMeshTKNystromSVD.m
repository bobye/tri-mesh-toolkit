function S = TriMeshTKNystromSVD( mat_A, mat_B )
% return Nystrom approximated SVD

[U, D, H] = svd(mat_A);
ZU = [U*sqrt(D); mat_B*H*sqrt(pinv(D))];
ZH = [H*sqrt(D); mat_B*U*sqrt(pinv(D))];

S=abs(eigs(ZH'*ZU, 100));  
end

