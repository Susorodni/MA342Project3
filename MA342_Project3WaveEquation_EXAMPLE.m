%% MA342 Project 3 Wave Equation EXAMPLE
clear variables;
close all;
clc;
            
% Definitions
L = 2;          % m; total length of rope
T = 30;         % s; total time duration
alpha = 0.1;    % ul; laplacian constant this is [H / rho]
dt = 0.028;     % s; time displacement increment
dx = 0.01;      % m; position displacement increment

% Iteratives
PosX = unique([0:dx:L, L]);     % define times to look at
Time = unique([0:dt:T, T]);     % define positions to look at
Sol = zeros(length(PosX), length(Time));

% Initial Conditions
Sol(:, 2) = (exp(-10*(PosX-(L/2)).^2) - exp(-10*(L/2).^2))*dt;

% NORMALLY we would define boundary conditions here. However, since
% they are zero they're already defined in the Sol definition, so
% we don't have to do that again.

% Expanded Laplacian Constant
r = (alpha*(dt^2))/(dx^2);

% Index array
indx = 2:length(PosX)-1; % don't want ends included as they are boundaries

% This iterates through time FOR EVERY position of X
for j = 3:length(Time)
    Sol(indx, j) = ...
        2*Sol(indx, j-1) ...
        - Sol(indx, j-2) ...
        + r*(Sol(indx+1, j-1) - 2*Sol(indx, j-1) + Sol(indx-1, j-1));
end

% In Sol, columns represent the position of x
% rows are for each iteration of time

% Movie Creation
loops = size(Sol, 1);

% Initialize frames
F(loops) = struct('cdata', [], 'colormap', []);

axis tight equal;
ax = gca;
ax.NextPlot = 'replaceChildren';

for j = 1:loops
    plot(PosX, Sol(:, j));
    F(j) = getframe(gcf);
end