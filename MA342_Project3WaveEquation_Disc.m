%% MA342 Project 3 Wave Equation Disc
clear variables;
close all;
clc;
            
% Definitions
r = 1;          % m; radius of disc
T = 10;         % s; total time duration
alpha = 2.0;    % ul; laplacian constant this is [H / rho]
dx = 0.05;      % m; position displacement increment
dy = dx;
dt = dx * dy * sqrt(alpha);

x = -r:dx:r;
y = x;

[X, Y] = meshgrid(x, y);
Time = unique([0:dt:T, T]);     % define positions to look at
% Sol = zeros(length(X), length(Y), length(Time));
movie_fps = 60;
total_frames = movie_fps*T;
movie_2d = true;
movie_3d = true;
snapshot_times = floor(linspace(0, length(Time), total_frames));

% Initial conditions
prev_iter = (X.^2 + Y.^2) / 10;
select_iter = prev_iter;

% initialize some movie stuff

if movie_2d || movie_3d
    h = figure;
    h.Visible = 'off';
end
if movie_3d
    surf(X, Y, select_iter, 'EdgeColor','none');
    M3d(total_frames) = struct('cdata',[],'colormap',[]);
    M3d(1) = getframe;
end
if movie_2d
    imagesc(x, y, select_iter);
    M2d(total_frames) = struct('cdata',[],'colormap',[]);
    M2d(1) = getframe;
end
if movie_2d || movie_3d
    colorbar;
    axis tight manual;
    xlim([-r, r]);
    ylim([-r, r]);
end
snap_index = 2;

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
            % Sol(i, j, :) = boundary_condition;
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

% Movie creation
loops = size(length(Time));

% Initialize frames
F(loops) = struct('cdata', [], 'colormap', []);

axis tight manual;
grid on;
xlim([-r, r]);
ylim([-r, r]);

ax = gca;
ax.NextPlot = 'replaceChildren';

for k = 1:total_frames

    if k >= 3
        % tic
        for i = 2:length(X)-1
            for j = 2:length(Y)-1
                if ~ban_list(i, j)
                    next_iter(i, j) = ...
                        2*select_iter(i, j) ...
                        - prev_iter(i, j) ...
                        + r_x*(select_iter(i + 1, j) - 2*select_iter(i, j) + select_iter(i - 1, j)) ...
                        + r_y*(select_iter(i, j + 1) - 2*select_iter(i, j) + select_iter(i, j - 1));
                end
            end
        end
        % toc

        prev_iter = select_iter;
        select_iter = next_iter;

    else
        switch k
            case 1
                next_iter = prev_iter;
            case 2
                next_iter = select_iter;
        end
    end

    if k == snapshot_times(snap_index)
        plot_iter = next_iter;

        % Remove outside circle pieces
        for i = 1:length(X)
            for j = 1:length(Y)
                if ban_list(i, j)
                    plot_iter(i, j, :) = NaN;
                end
            end
        end

        if movie_3d
            surf(X, Y, plot_iter, 'EdgeColor','none');
            drawnow
            M3d(snap_index) = getframe;
        end
        if movie_2d
            imagesc(x, y, plot_iter);
            drawnow
            M2d(snap_index) = getframe;
        end
        snap_index = snap_index + 1;
    end
end

% movie(s)
if movie_3d || movie_2d
    h.Visible = 'on';
end
if movie_3d
    movie(M3d, 1, movie_fps);
end
if movie_2d
    if movie_3d
        figure
    end
    movie(M2d, 1, movie_fps);
end