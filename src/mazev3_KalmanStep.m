function [TotalDW, net]=mazev3_KalmanStep(net, objMaze, DeltaJ, E, ERRORS, MAZES)
%% C++ version

adjusted_alpha = -1*net.alpha;  % 49 by 1

net.R  = eye(size(net.R)) .* net.inline_r(E); % mazev3_r(E, 0, 0);

% now regular Kalman filter
GAMMA   = net.C * net.K * net.C'+ net.R; % 49 by 49
G       = net.K * net.C' * inv(GAMMA);   % 86 by 49
TotalDW = G * adjusted_alpha;  % this is the weight update % 86 by 1
net.K = net.K - G * net.C * net.K + net.Q; % 86 by 86
