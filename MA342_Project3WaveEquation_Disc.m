%% MA342 Project 3 Wave Equation Disc
clear variables;
close all;
clc;

%% Disc Inital Formation
% Universal Constants

% The Laplacian constant, the unitless ratio between the 
% Horizontal forces and density of the increment
%       Unit: ul
alpha = 2.0;

% The increment in both the x and y directions for the
% mesh grid. Defines the "fineness"
%       Unit: m
d = 0.04;
            
% XY Mesh Definitions

% The radius of the unit disc
%       Unit: m
r = 1;

% The positions on both the x and y axis
%       Unit: ([] x 1) m
mesh_axis = transpose(-r:d:r);

% The number of position increments in both the x and y axis
%       Unit: ul
num_axis_inc = length(mesh_axis);

% Return the position for every index in the position matrix
%       Unit: ([] x []) m
[X, Y] = meshgrid(mesh_axis, mesh_axis);

% Time Definitions

% The total time duration
%       Unit: s
T = 5;

% The increment in time for every computed frame
%       Unit: s
dt = d^2 * sqrt(alpha);

% Every time to compute the disc position
%       Unit: ([] x 1)
Time = transpose(unique([0:dt:T, T]));

% The number of time increments
%       Unit: ul
num_time_inc = length(Time);

% The solution matrix that defines the z axis displacement at every
% position on the XY grid at every frame
%       Unit: ([] x [] x []) m
Sol = zeros(num_axis_inc, num_axis_inc, num_time_inc);

%% Initial Conditions

% Define the initial z axis displacement for every position
Sol(:, :, 1) = (X.^2 + Y.^2) / 10;

% The next time interval retains the same displacement
Sol(:, :, 2) = Sol(:, :, 1);

%% Boundary Conditions
% The unit disc has positions defined in a rectangular matrix. Because
% of this, the positions of the boundary of the outer circle are not
% perfect and have to be estimated.

% This defines the margin of error to determine if a location on the
% XY matrix is a boundary of the unit circle
%       Unit: m
boundary_precision = 0.05;

% This defines the constant z displacement that all unit circle boundary
% positions both on and outside it to be for all time.
%       Unit: m
boundary_condition = 0.1;

% To determine which positions will be treated with a constant position
% and which will not, the "ban list" makes a logical matrix for
% MATLAB to check both for computing the central difference when 
% doing the displacement calculations and movie plotting. This allows
% Proper computation to occur while ensuring that the boundary
% conditions do not change
ban_list = zeros(length(X), length(Y));

for i = 1:num_axis_inc
    for j = 1:num_axis_inc

        % Grab the distance each point is from the center, if
        % it's within the precision to r, or greater than one
        % outright, it's put on the ban list
        p_dist = X(i, j)^2 + Y(i, j)^2;
        diff = abs(r - p_dist);
        
        if diff < boundary_precision || p_dist > r
            % All positions here will be at the 
            % boundary condition for all time
            Sol(i, j, :) = boundary_condition;
            ban_list(i, j) = 1;
        end
    end
end

% Converts the ban list to a logical matrix
ban_list = logical(ban_list);

%% Computation

% These are the "Expanded Laplacian Constant." These are simply
% constants within the computation equation.
r_x = (alpha * dt^2) / d^2;
r_y = (alpha * dt^2) / d^2;

% Define the array of indexes to compute the central difference
% equation at in the for loop. The first and last positions will not
% be inspected since they would go beyond the scope of the central
% difference method
indx = 2:num_axis_inc-1;

% The for loop starts at 3 since the first two "frames" are already
% defined to be used for following frames
for k = 3:num_time_inc
    for i = indx
        for j = indx

            % We only compute the operation if the position
            % is not on the ban list
            if ~ban_list(i, j)
                Sol(i, j, k) = ...
                    2*Sol(i, j, k-1) ...
                    - Sol(i, j, k-2) ...
                    + r_x*(Sol(i + 1, j, k-1) - 2*Sol(i, j, k-1) + Sol(i - 1, j, k-1)) ...
                    + r_y*(Sol(i, j + 1, k-1) - 2*Sol(i, j, k-1) + Sol(i, j - 1, k-1));
            end
        end
    end
end

%% Plotting
% Once our matrix is fully computed, we can set all the positions
% on the "ban list" to be NaN so that they are invisible during
% any animation
for i = 1:num_axis_inc
    for j = 1:num_axis_inc
        if ban_list(i, j)
            Sol(i, j, :) = NaN;
        end
    end
end


% We have to define all the frames of the animation in the form
% of a struct by defining its size with the number of frames
F(num_time_inc) = struct('cdata', [], 'colormap', []);

% Gives good view of the plots
axis tight manual;
grid on;

% Keeps the limits of the XY meshgrid within the radius of the
% unit circle
xlim([-r, r]);
ylim([-r, r]);

% Finds the absolute highest and lowest points the unit circle goes
% in order to constrain the height to within this range
z_max = max(Sol, [], 'all');
z_min = min(Sol, [], 'all');

% Obtain the properties of the plot
ax = gca;
ax.NextPlot = 'replaceChildren';

% Plot the 3D matrix
for j = 1:num_time_inc
    surf(X, Y, Sol(:, :, j));
    zlim([z_min, z_max]);
    F(j) = getframe(gcf);
    sgtitle(append("Disc at t = ", sprintf('%5.4f', Time(j))))
end

