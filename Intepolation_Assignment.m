%% ----Assigment4--Interpolation
%% coded by Veasna Mom-25M17722


%% Answer a

clear; clc;

% Define the function to interpolate
f = @(x) 1./(1+ 25*x.^2);

% Define query points for smooth plotting
x = linspace(-1, 1, 200);

%
N = input('Enter number of interpolation nodes: ');  % This is dynamic input
% --- Perform interpolation ---
Pn_LG = LG_interpolate(f, x, N);           % Lagrange
Pn_NDD = Newton_Divided_Diff(f, x, N);     % Newton Divided Difference

m = input('Enter degree of least squares polynomial: ');
Pn_LS = Least_Squares(f, x, N, m);         % Least Squares


%%  --- Plot results ---
figure;
plot(x, f(x), 'k', 'LineWidth', 2); hold on;
plot(x, Pn_LG, 'r--', 'LineWidth', 1.5);
plot(x, Pn_NDD, 'b:', 'LineWidth', 1.5);
plot(x, Pn_LS, 'g-', 'LineWidth', 1.5);
legend('Exact f(x)', 'Lagrange', 'Newton', 'Least Squares', 'Location', 'Best');
xlabel('x');
ylabel('f(x)');
title(['Interpolation & Approximation with N = ', num2str(N), ', LS degree = ', num2str(m)]);
grid on;

%% Answer b --- Calculate relative error for Least Squares at interpolation nodes ---
x_nodes = linspace(-1, 1, N);
y_true = f(x_nodes);
y_pred = Least_Squares(f, x_nodes, N, m);

relative_error = abs(y_true - y_pred)./abs(y_true);

% --- Write relative errors to ASCII file ---
filename = "D:\M.TECH\CMTFE\Assignment4\relative_error.txt";
fileID = fopen(filename, 'w');

fprintf(fileID, 'x\tTrueValue\tPredictedValue\tRelativeError\n');
for i = 1:N
    fprintf(fileID, '%.6f\t%.6f\t%.6f\t%.6f\n', x_nodes(i), y_true(i), y_pred(i), relative_error(i));
end
fclose(fileID);

disp(['Relative error data saved to ', filename]);



%% --------- FUNCTION: Lagrange Interpolation -------------
% xq is simple x that we want to interpolate
function Pn = LG_interpolate(f, xq, N)
    x_nodes = linspace(-1, 1, N);
    y_nodes = f(x_nodes); % This give value of all y_i we want to do L*y(x)
    Pn = zeros(size(xq));
    for k = 1:length(xq) % For how many x we want to do interpolation
        P = 0;
        for i = 1:N
            Li = 1;
            for j = 1:N
                if j ~= i
                    Li = Li * (xq(k) - x_nodes(j)) / (x_nodes(i) - x_nodes(j));
                end
            end
            P = P + y_nodes(i) * Li;
        end
        Pn(k) = P;
    end
end

%% --------- FUNCTION: Newton Divided Difference -------------
function Pn_ndd = Newton_Divided_Diff(f, xq, N)
    x_nodes = linspace(-1, 1, N);
    y_nodes = f(x_nodes);

    % Compute divided difference coefficients
    a = y_nodes;
    for j = 2:N
        for i = N:-1:j
            a(i) = (a(i) - a(i-1)) / (x_nodes(i) - x_nodes(i-j+1));
        end
    end

    % Evaluate polynomial at query points
    Pn_ndd = zeros(size(xq));
    for k = 1:length(xq)
        value = a(N);
        for j = N-1:-1:1
            value = a(j) + (xq(k) - x_nodes(j)) * value;
        end
        Pn_ndd(k) = value;
    end
end

%% --------- FUNCTION: Least Squares Polynomial -------------
function Pn_ls = Least_Squares(f, xq, N, m)
    x_nodes = linspace(-1, 1, N);
    y_nodes = f(x_nodes);

    % Build the normal equations for polynomial coefficients
    A = zeros(m+1, m+1);
    b = zeros(m+1, 1);

    for j = 0:m
        for k = 0:m
            A(j+1, k+1) = sum(x_nodes.^(j+k));
        end
        b(j+1) = sum(y_nodes .* (x_nodes.^j));
    end

    % Solve for coefficients
    coeff = A \ b; % Here I used builtin matrix solving for easiness

    % Evaluate polynomial at query points
    Pn_ls = zeros(size(xq));
    for k = 1:length(xq)
        value = 0;
        for j = 0:m
            value = value + coeff(j+1) * xq(k)^j;
        end
        Pn_ls(k) = value;
    end
end
