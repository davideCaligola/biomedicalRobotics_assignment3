clear all;

% fixed step time
time_step = 0.05;
Ts = time_step;

% threshold for reaching target / home center
th = 0.005;

% define home and target coordinates
n_target = 8;
radius = 0.1;
coord_home = [0,0];
coord_targets = zeros(n_target+2,2);
step = 2*pi/n_target;
for i=0:n_target-1
    coord_targets(i+1,:) = [radius*cos(i*step), radius*sin(i*step)];
end

% randomize coordiantes order
rand_idx = randi(n_target,[1,10]);
coord_targets = coord_targets(rand_idx,:);
