%%
% Generates a force field to attract the manipulator end-effector to a goal
% position proportional to the distance and an anisotorpic viscousity,
% stronger on the orthogonal direction to the distance.
% If the end-effector position is outside a threshold with respect to the
% ideal linear trajectory, a resistive force  proportional with respect to 
% the distance from the ideal trajectory is applied to bring back the
% end-effector closer to the ideal trajectory.
% Home position is considered the origin at (0,0)
%
% NOTE:
% position in meter [m]
% velocity in meter/second [m/s+
%
% Notation:
%   p_ee = end-effector position (vector, 3x1)
%   v_ee = end-effector velocity (vector, 3x1)
%   p_g = goal position (vector, 3x1)
%   v_g = goal velocity (vector, 3x1)
%   a_pp = position error gain parallel (scalar)
%   a_po = position error gain orthogonal (scalar)
%   a_vp = viscosity gain parallel (scalar)
%   a_vo = viscosity gain orthogonal (scalar)
%   
%   p_margin = margin before apply direction compensation
%
% The attractive force is considered only on the plane, thus the following
% projector is considered:
%        
%   Rxy = [1 0 0;
%          0 1 0;
%          0 0 0]
%
% The coordinates will be projected in the plane using such an operator.
% Define distance frame with x-axis along the direction of the distance
% vector between the goal and home, y-axis orthogonal to it.
%
%   cos(theta) = p_g(1) / norm(p_g)
%   sin(theta) = p_g(2) / norm(p_g)
%
% Rotation matrix to move from world frame to distance frame and world
% frame:
%
%   dRw = [ cos(theta) sin(theta) 0
%          -sin(theta) cos(theta) 0
%                0          0     1]
%
% From now on it is considered to work in distance frame
% The attractive potential is defined as follows:
%
%   U_att = (1/2)*(p_g-p_ee)'*A_p*(p_g-p_ee) + (1/2)*(v_g-v_ee)'*A_v*(v_g-v_ee)
%
% so that the defined Uatt has a minimum at the goal position at zero velocity.
% The attractive force will be the opposite of the attractive potential
% gradient with respect to pee and vee, thus
%
%   F_att = - ∇p_ee(U_att) - ∇v_ee(U_att) = a_p(p_g - p_ee) + a_v(v_g - v_ee)
%
% By hypothesis the goal velocity v_g is null => v_g = 0.
% The resulting actractive force results as follows:
%
%   F_att = (a_p*(p_g-p_ee) - a_v*v_ee)
%
% Thresholds (p_th, v_th) are also considered, after which both
% the position and the velocity errors are considered small enough and the
% target reached.
% The found force needs to be converted into the world frame
%
%   wF_att = dRw'*F_att
%

function wFatt = forceField(p_ee, v_ee, p_g, a_pp, a_po, a_vp, a_vo, p_th, v_th, p_margin)
    
    % check if input vectors are rows. In case convert them to column
    if isrow(p_ee)
        p_ee = p_ee';
    end

    if isrow(v_ee)
        v_ee = v_ee';
    end

    if isrow(p_g)
        p_g = p_g';
    end


    % planar projector
    Rxy = [1 0 0;
           0 1 0;
           0 0 0];

    % project goal and end-effector coordinates on the plane
    p_g = Rxy*p_g;
    p_ee = Rxy*p_ee;
    v_ee = Rxy*v_ee;

    % planar rotation angle between world and distance frame:
    cos_theta = p_g(1) / norm(p_g);
    sin_theta = p_g(2) / norm(p_g);

    % rotation frame from world to distance frame
    dRw = [ cos_theta sin_theta 0;
           -sin_theta cos_theta 0;
                0         0     1];

    % move from world to distance frame
    dp_g = dRw * p_g;
    dp_ee = dRw * p_ee;
    dv_ee = dRw * v_ee;

    % attractive force
    wFatt = zeros(3,1);

    % position error matrix gain
    A_p = zeros(3);
    A_p(1,1) = a_pp;
    % check if distance from ideal trajectory is greater than the margin
    if (abs(dp_ee(2)) > p_margin)
        A_p(2,2) = a_po;
    end

    % viscous matrix gain
    A_v = zeros(3);
    A_v(1,1) = a_vp;
    A_v(2,2) = a_vo;

    if (norm(Rxy*(dp_g - dp_ee)) > p_th) || (norm(Rxy*dv_ee) > v_th)
        % errors greater than thresholds calculate the attractive force
        % in world reference
        wFatt = dRw' * (A_p*(dp_g - dp_ee) - A_v*(dv_ee));
    end
end
