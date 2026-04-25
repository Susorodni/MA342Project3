%% MA342 Project 3 Wave Equation Disc
clear variables;
close all;
clc;
            
% Definitions
r = 1;          % m; radius of disc
T = 5;         % s; total time duration
alpha = 2.0;    % ul; laplacian constant this is [H / rho]
dx = 0.05;      % m; position displacement increment
dy = dx;
dt = dx * dy * sqrt(alpha);
p_vals = [1; 3; 5; 10];

x = -r:dx:r;
y = x;

[X, Y] = meshgrid(x, y);
Time = unique([0:dt:T, T]);     % define positions to look at
Sol = zeros(length(X), length(Y), length(Time), length(p_vals));

% Movie stuff
fps = 60;
total_movie_frames = T * fps;
movie_snapshot_times = floor(linspace(0, length(Time), total_movie_frames));

% Initial conditions
for i = 1:length(X)
    for j = 1:length(Y)
        for v = 1:length(p_vals)
            Sol(i, j, 1, v) = norm([X(i, j) Y(i, j)], p_vals(v)) / 10;
        end
    end
end
Sol(:, :, 2, :) = Sol(:, :, 1, :);

% how close the boundary conditions
% should be calculated
boundary_precision = 0.05;
boundary_condition = 0.1;

ban_list = zeros(length(X), length(Y));

% Boundary conditions
for i = 1:length(X)
    for j = 1:length(Y)

        d = X(i, j)^2 + Y(i, j)^2;
        diff = abs(1 - d);


        if diff < boundary_precision || d > 1
            Sol(i, j, :, :) = boundary_condition;
            ban_list(i, j) = 1;
        end
    end
end

ban_list = logical(ban_list);

% Expanded Laplacian Constant
r_x = (alpha * dt^2) / dx^2;
r_y = (alpha * dt^2) / dy^2;

% Index array
indx = 2:length(X)-1;

for k = 3:length(Time)
    for i = 2:length(X)-1
        for j = 2:length(Y)-1
            if ~ban_list(i, j)
                Sol(i, j, k, :) = ...
                    2*Sol(i, j, k-1, :) ...
                    - Sol(i, j, k-2, :) ...
                    + r_x*(Sol(i + 1, j, k-1, :) - 2*Sol(i, j, k-1, :) + Sol(i - 1, j, k-1, :)) ...
                    + r_y*(Sol(i, j + 1, k-1, :) - 2*Sol(i, j, k-1, :) + Sol(i, j - 1, k-1, :));
            end
        end
    end
end

% Add NaNs to make it look goofdf
for i = 1:length(X)
    for j = 1:length(Y)
        if ban_list(i, j)
            Sol(i, j, :, :) = NaN;
        end
    end
end


% Movie Creation
loops = size(Sol, 3);

% Initialize frames
F(loops) = struct('cdata', [], 'colormap', []);

axis tight manual;
grid on;
xlim([-r, r]);
ylim([-r, r]);

z_max = max(Sol, [], 'all');
z_min = min(Sol, [], 'all');

zlim([z_min, z_max]);

ax = gca;
ax.NextPlot = 'replaceChildren';

num_p_vals = length(p_vals);

for j = 1:loops

    for v = 1:num_p_vals
        figure(1);

        subplot(num_p_vals / 2, num_p_vals / 2, v);
        surf(X, Y, Sol(:, :, j, v));
        zlim([z_min, z_max]);
        title(append("P-Val = ", sprintf('%.0f', p_vals(v))));
        F(j) = getframe(gcf);
   
        % exportgraphics(gcf, 'simple_disc.gif', 'Append', true);

        % figure(2);
        % subplot(num_p_vals / 2, num_p_vals / 2, v);
        % imagesc(x, y, Sol(:, :, j, v));
        % drawnow;
    end

    sgtitle(append("Disc at t = ", sprintf('%5.4f', Time(j))))

    % if any(j == movie_snapshot_times)
    %     exportgraphics(gcf, 'disc_p_val.gif', 'Append', true);
    % end
end

