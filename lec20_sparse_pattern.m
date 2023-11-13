% sample run: lec17_jacobi(100);

function error=lec20_sparse_pattern(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: initialize the grid and initialize u^0_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the grid
x_left = 0.0;
x_right = 1.0;
x = linspace(x_left,x_right,N+1);
% x is a row vector, we want u to be a column vetor so take transpose here
x = x';
x_interior = x(2:end-1);

dx = 1/N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Assemble spatial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N,N);
for j = 1:N
    if (j>1)
        A(j,j-1) = -1;
    end

    if (j<N)
        A(j,j+1) = -1;
    end

    A(j,j) = 2;
end
%A = sparse(A);

close all
figure(1)
spy(inv(A))

figure(2)
spy(A)

end