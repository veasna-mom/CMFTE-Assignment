%% -------- Assignment 5 --Coded by Veasna-25M17722
clear;
clc;

%% Parameters
w = 1.0;              % length of rod [m]
k = 400;              % thermal conductivity [W/mK]
cp = 385;             % specific heat [J/kgK]
rho = 8000;           % density [kg/m^3]
T0 = 25;              % initial temperature [°C]
dx = 0.1;             % spatial step [m]
dt = 1;               % time step [s]

alpha = k/(rho*cp);   % thermal diffusivity [m^2/s]
N = round(w/dx) + 1;  % number of nodes
T_base = 400;         % Bc at x=0
T_end = 400;          % BC at x=L
tol = 1e-4;           % tolerance 
maxIter = 900;        %iterations
mid = round((N+1)/2); % midpoint index




%% a Run all schemes
[T_mid_ftcs, time_ftcs,t_reach] = FTCS(dx, N, dt, tol, T0, alpha, T_base, T_end, maxIter);
[T_mid_btcs, time_btcs] = BTCS(dx, N, dt, tol, T0, alpha, T_base, T_end, maxIter);
[T_mid_dufort, time_dufort] = DuFortFrankel(dx, N, dt, tol, T0, alpha, T_base, T_end, maxIter);
 %% b)
    fprintf('>> b. Midpoint reaches 200 °C at t = %.2f s (FTCS)\n', t_reach);
%% Plot comparison
figure;
plot(time_ftcs, T_mid_ftcs, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Explicit FTCS');
hold on;
plot(time_btcs, T_mid_btcs, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Implicit BTCS');
plot(time_dufort, T_mid_dufort, 'c--', 'LineWidth', 1.5, 'DisplayName', 'Implicit DuFort-Frankel');
xlabel('Time [s]');
ylabel('Midpoint Temperature [°C]');
title('Midpoint Temperature vs Time');
grid on;
legend;
hold off;

%% 1. Explicit FTCS
function [T_mid, time,t_reach] = FTCS(dx, N, dt, tol, T0, alpha, T_base, T_end, maxIter)
    Fo = alpha*dt/(dx^2); % Fourier number
    % Initial condition
    T = ones(N,1)*T0;
    T_mid = zeros(maxIter,1); % Preallocate memory to store midpoint temperature
    % Apply boundary conditions
    T(1) = T_base;
    T(end) = T_end;
    t_reach = NaN;  % store when midpoint reaches 200 °C
    for n = 1:maxIter
        T_old = T;
        % Update interior nodes
        for i = 2:N-1
            T(i) = T_old(i) + Fo*(T_old(i+1) - 2*T_old(i) + T_old(i-1));
        end
        % Store midpoint temperature
        T_mid(n) = T(round((N+1)/2));
        
     
        % Check target temperature (200 °C) For question b
        if isnan(t_reach) && T_mid(n) >= 200
            t_reach = n*dt; % record first hitting time
        end
        
        
        
        
        % Check convergence (steady state)
        if max(abs(T - T_old)) < tol
            fprintf('FTCS converged at iteration %d\n', n);
            T_mid = T_mid(1:n); % Trim unused zeros
            break
        end
    end
    time = (0:length(T_mid)-1)*dt; % Time array starting from t=0
end

%% 2.Implicit BTCS Scheme
function [T_mid, time] = BTCS(dx, N, dt, tol, T0, alpha, T_base, T_end, maxIter)
    a = alpha*dt/(dx^2); % Coefficient for matrix
    % Initial condition
    T = ones(N,1)*T0;
    T(1) = T_base; % Boundary condition at x=0
    T(end) = T_end; % Boundary condition at x=L
    T_mid = zeros(maxIter,1); % Preallocate memory to store midpoint temperature

    % Construct tridiagonal matrix A
    A = zeros(N,N);
    A(1,1) = 1; % Boundary condition at x=0
    A(N,N) = 1; % Boundary condition at x=L
    % Assign coefficients for interior nodes
    for i = 2:N-1
        A(i,i-1) = -a;  % Here I write  -aT_(i-1)+(1+2a)T_i+a*T_i+1 = Ti_old
        A(i,i) = 1 + 2*a;
        A(i,i+1) = -a;
    end

    % Time-stepping loop to solve the system of equations
    for n = 1:maxIter
        T_old = T;
        b = T_old; % Right-hand side(T at previous time step)
        b(1) = T_base; % Dirichlet boundary condition at x=0
        b(end) = T_end; % Dirichlet boundary condition at x=L
        T = A\b; % Solve linear system A*T = b for new temperatures
        % Store midpoint temperature
        T_mid(n) = T(round((N+1)/2));
        % Check convergence
        if max(abs(T - T_old)) < tol
            fprintf('BTCS converged at iteration %d\n', n);
            T_mid = T_mid(1:n); % Trim unused zeros
            break;
        end
    end
    time = (0:length(T_mid)-1)*dt; % Time array starting from t=0
end

function [T_mid, time, T_final] = DuFortFrankel(dx, N, dt, tol, T0, alpha, T_base, T_end, maxIter)
    

    Fo = alpha*dt/dx^2;   % Fourier number

    % Initial condition
    T = ones(N,1)*T0;
    T(1) = T_base;
    T(end) = T_end;

    % Storage
    T_mid = zeros(maxIter,1);
    time  = zeros(maxIter,1);

    % First step: FTCS
    T_ftcs = T;
    for i = 2:N-1
        T_ftcs(i) = T(i) + Fo*(T(i+1) - 2*T(i) + T(i-1));
    end
    T_ftcs(1) = T_base;
    T_ftcs(end) = T_end;

    % Time levels
    T_prev = T;       % u^0
    T_curr = T_ftcs;  % u^1
    T_new  = T_curr;

    % Time stepping
    for n = 1:maxIter
        for i = 2:N-1
            T_new(i) = ((1 - 2*Fo)/(1 + 2*Fo))*T_prev(i) +(2*Fo/(1 + 2*Fo))*(T_curr(i+1) + T_curr(i-1));
        end
        T_new(1)   = T_base;
        T_new(end) = T_end;

        % Save midpoint & time
        T_mid(n) = T_curr(round(N/2));
        time(n)  = n*dt;

        % Convergence check
        if max(abs(T_new - T_curr)) < tol
            fprintf('DuFortFrankel converged at iteration %d\n', n);
            T_mid = T_mid(1:n);
            time  = time(1:n);
            break;
        end

        % Shift time levels
        T_prev = T_curr;
        T_curr = T_new;
    end

    T_final = T_curr; % return last profile
end
