function M = TriMeshTKFBase2M( f, L, D )
% assemble Fourier simplifed matrix

n = length(D);
M=zeros(n,n);

K= L'*f(D(2:n));
M(2:n,1) = K;
M(1,2:n) = K';
M(1,1) = 2*sum(f(D(2:n)));
for i=2:n
    M(i,i) = -2* f(D(i));
end

end

