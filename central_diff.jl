using LinearAlgebra;

r = 1;      # radius of disk
T = 1000;   # total time analysis
alpha = 2.0;
dx = 0.05;
dy = dx;
dt = dx * dy * sqrt(alpha);

x = -r:dx:r;
y = x;

X = x' .* ones(length(x));
Y = ones(length(y)) .* y';

Time = unique(transpose(0:dt:T));
Time = vcat(Time, T);
# Sol = zeros(length(X), length(Y), length(Time));

#=
function central_diff(Sol, indx, Time, r_x, r_y)
    for j in 3:length(Time)
        Sol[indx, indx, j] = 
            2*Sol[indx, indx, j-1]
            - Sol[indx, indx, j-2]
            + r_x*(Sol[indx + 1, indx, j-1] - 2*Sol[indx, indx, j-1] + Sol[indx - 1, indx, j-1])
            + r_y*(Sol[indx, indx + 1, j-1] - 2*Sol[indx, indx, j-1] + Sol[indx, indx - 1, j-1])
    end
end
=#