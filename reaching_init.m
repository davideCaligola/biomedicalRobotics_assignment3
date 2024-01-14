clear all;

% fixed step time
time_step = 0.05;

% threshold for reaching target / home center
th = 0.005;

% define home and target coordinates
n_target = 8;
coord_home = [0,0];
coord_targets = zeros(n_target,2);
step = 2*pi/n_target;
for i=0:n_target-1
    coord_targets(i+1,:) = [0.1*cos(i*step), 0.1*sin(i*step)];
end

% randomize coordiantes order
rand_idx = randperm(n_target);
coord_targets(rand_idx,:) = coord_targets;
