function eigen_solver(varargs)

num_eigval = 5;
output_boolean = false;

if (nargin == 1)   
    L = PetscBinaryRead(varargs{1});
    M = [];
elseif (nargin == 2)
    if (isscalar(varargs{2}))
        L = PetscBinaryRead(varargs{1});
        M = [];
        num_eigval = varargs{2};
    else
        L = PetscBinaryRead(varargs{1});
        M = PetscBinaryRead(varargs{2});
    end
elseif (nargin == 3)
    L = PetscBinaryRead(varargs{1});
    M = PetscBinaryRead(varargs{2});    
    num_eigval = varargs{3};
    output_boolean = true;    
end

%L = PetscBinaryRead('~/data/meshtk_workshop/sample.stiff');
%D = PetscBinaryRead('~/data/meshtk_workshop/sample.mass');

opts.issym = 1;
opts.isreal = 1;
tic;
if (M)
    [V, D, flag] = eigs(L,M,eigval,'sm',opts);
else
    [V, D, flag] = eigs(L,num_eigval,'sm',opts);
end
toc;

if (flag)
    fprintf('not all eigenvalues are converged\n');
end

[D, IDX] = sort(diag(D));
V = V(:,IDX);

fprintf('%.12G\n',D);

if (output_boolean)
    output_filename = sprintf('%s/_ev.ascii', varargs{4});
    export(output_filename, D);
    
    n = length(D);
    for i=1:n           
        output_filename = sprintf('%s/_%d.petsc', varargs{4},i-1);
        PetscBinaryWrite(output_filename, V(i));
    end
end
