%% MA342 Project 3 Wave Equation Disc
clear variables;
close all;
clc;
            
% Definitions
r = 1;          % m; radius of disc
T = 120;         % s; total time duration
alpha = 2.0;    % ul; laplacian constant this is [H / rho]
dx = 0.05;      % m; position displacement increment
dy = dx;
dt = dx * dy * sqrt(alpha);
p_val = 10;

x = -r:dx:r;
y = x;

[X, Y] = meshgrid(x, y);
Time = unique([0:dt:T, T]);     % define positions to look at
Sol = zeros(length(X), length(Y), length(Time));

% Initial conditions
for i = 1:length(X)
    for j = 1:length(Y)
        Sol(i, j, 1) = norm([X(i, j) Y(i, j)], p_val) / 10;
    end
end
Sol(:, :, 2) = Sol(:, :, 1);

% how close the boundary conditions
% should be calculated
boundary_precision = 0.1;
boundary_condition = 0.1;

% Boundary conditions
for i = 1:length(X)
    for j = 1:length(Y)

        d = norm([X(i, j) Y(i, j)], p_val);
        diff = abs(1 - d);


        if diff < boundary_precision
            Sol(i, j, :) = boundary_condition;
        elseif d > 1
            Sol(i, j, :) = boundary_condition;
        end
    end
end

% Expanded Laplacian Constant
r_x = (alpha * dt^2) / dx^2;
r_y = (alpha * dt^2) / dy^2;

% Index array
indx = 2:length(X)-1;

for j = 3:length(Time)
    Sol(indx, indx, j) = ...
        2*Sol(indx, indx, j-1) ...
        - Sol(indx, indx, j-2) ...
        + r_x*(Sol(indx + 1, indx, j-1) - 2*Sol(indx, indx, j-1) + Sol(indx - 1, indx, j-1)) ...
        + r_y*(Sol(indx, indx + 1, j-1) - 2*Sol(indx, indx, j-1) + Sol(indx, indx - 1, j-1));
end

% Movie Creation
loops = size(Sol, 1);

% Initialize frames
F(loops) = struct('cdata', [], 'colormap', []);

axis tight;
ax = gca;
ax.NextPlot = 'replaceChildren';

for j = 1:loops
    surf(X, Y, Sol(:, :, j))
    F(j) = getframe(gcf);
end

