function S = TriMeshTKNystromSVD( mat_A, mat_B )
% return Nystrom approximated SVD

[U, D, H] = svd(mat_A);
ZU = [U*sqrt(D); mat_B*H*sqrt(inv(D))];
ZH = [H*sqrt(D); mat_B*U*sqrt(inv(D))];

S=svd(ZH'*ZU);
end

