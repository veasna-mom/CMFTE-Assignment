% ------------------------------------------------------------
% Gauss-Seidel and SOR Method 
% Solves Ax = b using two iterative methods
% ------------------------------------------------------------

% Define system of equations
A = [4 -1  0;
    -1  4 -1;
     0 -1  3];   % This is my good way to write matrix
b = [15;
     10; 
     10];

% Solver parameters
tol = 1e-5;           % Convergence tolerance
maxIter = 1000;       % Maximum allowed iterations
x0 = zeros(length(b),1);   % Initial guess (start at zero) you can start with any input

% ------------------- Run Gauss-Seidel -----------------------
[x_gs, iter_gs, res_gs] = GaussSeidel(A,b,x0,tol,maxIter);
disp(" Gauss-Seidel Method")
disp("   Solution vector x = ")
disp(x_gs)
disp("   Converged in " + iter_gs + " iterations")

% ------------------- Run SOR -------------------------------
omega = 1.11;  % Relaxation factor (1 < omega < 2 for over-relaxation), you may try with different omega see what happen? 
[x_sor, iter_sor, res_sor] = MySOR(A,b,x0,tol,maxIter,omega);
disp(" SOR Method (omega = " + omega + ")")
disp("   Solution vector x = ")
disp(x_sor)
disp("   Converged in " + iter_sor + " iterations")

% ------------------- Plot results ---------------------------
figure;
semilogy(res_gs,'b-o','LineWidth',1.2); hold on;
semilogy(res_sor,'r-s','LineWidth',1.2);
grid on;
xlabel('Iteration'); ylabel('Residual (||Ax-b||_\infty)');
legend('Gauss-Seidel','SOR','Location','northeast');
title('Convergence Comparison: Gauss-Seidel vs SOR');

% ============================================================
% Function: Gauss-Seidel
% ============================================================
function [x,Iter,residuals] = GaussSeidel(A,b,x0,tol,maxIter)
    n = length(b);
    x = x0;
    residuals = zeros(maxIter,1);
    
    for k = 1:maxIter
        x_old = x;
        
        % Update each variable in sequence
        for i = 1:n
            sum1 = A(i,1:i-1)*x(1:i-1);
            sum2 = A(i,i+1:n)*x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2)/A(i,i);
        end
        
        % Compute residual (largest imbalance in Ax=b)
        residuals(k) = norm(A*x-b,inf);
        
        % Convergence check (largest change in x)
        if norm(x - x_old,inf) < tol
            Iter = k;
            residuals = residuals(1:k); % trim unused entries
            return;
        end
    end
    
    % If maxIter reached without convergence
    Iter = maxIter;
    residuals = residuals(1:maxIter);
    warning("Gauss-Seidel did not converge within max iterations");
end

% ============================================================
% Function: SOR (Successive Over-Relaxation)
 
  % SOR can help you to reach convergence faster
  % SOR can also take longer to converge
  % SOR use omega= Relaxation Factor
  % If omega>=1 it is called over relaxation
  % If 0<omega <1 it is called under-relaxation
  % Imagine you are simulating compressible flow, the density change
  % suddenly this lead to fluctuting floating point error in CFD code we
  % use under relaxation to ensure that code work well, I use OpenFOAM for
  % High Speed Flow it often be a good practice
  % Next time when you see CFD code, when dealing with numerical oscilation
  % this is a good way to implement this technique
  % At IIT Bombay, I finds this course is very very challenging




% ============================================================
function [x,Iter,residuals] = MySOR(A,b,x0,tol,maxIter,omega)
    n = length(b);
    x = x0;
    residuals = zeros(maxIter,1);
    
    for k = 1:maxIter
        x_old = x;
        
        % Update each variable in sequence
        for i = 1:n
            sum1 = A(i,1:i-1)*x(1:i-1);
            sum2 = A(i,i+1:n)*x_old(i+1:n);
            x_gs = (b(i) - sum1 - sum2)/A(i,i); % GS update
            x(i) = (1-omega)*x_old(i) + omega*x_gs; % Relaxed update
        end
        
        % Compute residual
        residuals(k) = norm(A*x-b,inf); % Why I use norm here?  because our sol is in vector representation
        
        % Convergence check
        if norm(x - x_old,inf) < tol 
            Iter = k;
            residuals = residuals(1:k);
            return;
        end
    end
    
    Iter = maxIter;
    residuals = residuals(1:maxIter);
    warning("SOR did not converge within max iterations");
end
