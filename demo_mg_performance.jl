 
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using AlgebraicMultigrid 
using Random
Random.seed!(10086)


 function test_solver(N,solver_option)   
    #########################################
    # Step 1: set up the solver and the mesh
    #########################################
    tol = 1e-7
    
    h = 1.0/(N+1)
    
    ###################################################
    # Step 2: set up the discrete Laplacian operator
    ##################################################
    # central difference
    e_vec = ones(N)
    A = spdiagm(0 => (-2*e_vec))+spdiagm(1 => (e_vec[1:end-1]))+spdiagm(-1 => (e_vec[1:end-1]))
    A .*= (1.0/h^2)
    
    # 3D
    A = kron(kron(A,sparse(I, N, N)),sparse(I, N, N))+
            kron(sparse(I, N, N),kron(A,sparse(I, N, N)))+
            kron(sparse(I, N, N),kron(sparse(I, N, N),A))
    M = N*N*N

    # RHS
    x_exact = rand(M);
    rhs = A*x_exact;

    #ml = ruge_stuben(A,strength = Classical(0.5))
    # Define AMG solver
    ml = ruge_stuben(A)
   

    if (solver_option==1)
        p = aspreconditioner(ml,AlgebraicMultigrid.V());
        time = @elapsed sol, cg_log = cg(A,rhs, Pl=p, log=true, reltol=1e-8);
    elseif (solver_option==2)
        p = aspreconditioner(ml,AlgebraicMultigrid.W());
        time = @elapsed sol, cg_log = cg(A,rhs, Pl=p, log=true, reltol=1e-8);
    else
        time = @elapsed sol, cg_log = cg(A,rhs, log=true, reltol=1e-8);
    end
    return time,cg_log.iters,norm(sol-x_exact,Inf)
 end


N_array = [50;100;200]
solver_option = 3
println("8*log(8) = ",8*log(8)," 8*sqrt(8) = ",8*sqrt(8))


for i = 1:length(N_array)
    N = N_array[i];

    time,iter,error = test_solver(N,solver_option)

    if (solver_option == 1)
        println("V cycle AMG preconditioner")
        println(" DOFs = ",N^3," Time = ",time," Iter = ",iter," l-2 error = ",error)
    elseif (solver_option == 2)
        println("W cycle AMG preconditioner")
        println(" DOFs = ",N^3," Time = ",time," Iter = ",iter," l-2 error = ",error)
    else
        println("CG")
        println(" DOFs = ",N^3," Time = ",time," Iter = ",iter," l-2 error = ",error)
    end
    println()
end