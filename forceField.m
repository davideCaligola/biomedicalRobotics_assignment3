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
% velocity in meter/second [m/s]
%
% Notation:
%   p_ee = end-effector position (vector, 3x1)
%   v_ee = end-effector velocity (vector, 3x1)
%   p_end = goal position (vector, 3x1)
%   v_end = goal velocity (vector, 3x1)
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
% vector between the target point (p_end) and the starting point (p_start), 
% y-axis orthogonal right-handed to it.
%
%   cos(theta) = (p_end(1) - p_start(1)) / norm(p_end - p_start)
%   sin(theta) = (p_end(2) - p_start(2)) / norm(p_end - p_start)
%
% Rotation matrix to move from world frame to distance frame:
%
%   dRw = [ cos(theta) sin(theta) 0
%          -sin(theta) cos(theta) 0
%                0          0     1]
%
% From now on it is considered to work in distance frame
% The attractive potential is defined as follows:
%
%   U_att = (1/2)*(p_end-p_ee)'*A_p*(p_end-p_ee) + (1/2)*(v_end-v_ee)'*A_v*(v_end-v_ee)
%
% so that the defined Uatt has a minimum at the goal position at zero velocity.
% The attractive force will be the opposite of the attractive potential
% gradient with respect to pee and vee, thus
%
%   F_att = - ∇p_ee(U_att) - ∇v_ee(U_att) = A_p(p_end - p_ee) + A_v(v_end - v_ee)
%
% By hypothesis the goal velocity v_end is null => v_end = 0.
% The resulting actractive force results as follows:
%
%   F_att = (A_p*(p_end-p_ee) - A_v*v_ee)
%
% Thresholds (p_th, v_th) are also considered, after which both
% the position and the velocity errors are considered small enough and the
% target reached.
% The found force needs to be converted into the world frame
%
%   wF_att = dRw'*F_att
%

function wF_att = forceField_aniso(enable, p_ee, v_ee, p_start, p_end, a_pp, a_po, a_vp, a_vo, p_th, v_th, p_margin)
    
    if enable

        % check if input vectors are rows. In case convert them to column
        if isrow(p_ee)
            p_ee = p_ee';
        end
    
        if isrow(v_ee)
            v_ee = v_ee';
        end
    
        if isrow(p_start)
            p_start = p_start';
        end
    
        if isrow(p_end)
            p_end = p_end';
        end
    
    
        % planar projector
        Rxy = [1 0 0;
               0 1 0;
               0 0 0];
    
        % project goal and end-effector coordinates on the plane
        p_end = Rxy*p_end;
        p_ee = Rxy*p_ee;
        v_ee = Rxy*v_ee;
    
        % planar rotation angle between world and distance frame:
        cos_theta = (p_end(1) - p_start(1)) / norm(p_end - p_start);
        sin_theta = (p_end(2) - p_start(2)) / norm(p_end - p_start);
    
        % rotation frame from world to distance frame
        dRw = [ cos_theta sin_theta 0;
               -sin_theta cos_theta 0;
                    0         0     1];
    
        % move from world to distance frame
        dp_end = dRw * p_end;
        dp_ee = dRw * p_ee;
        dv_ee = dRw * v_ee;
    
        % attractive force
        wF_att = zeros(3,1);
    
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
    
        if (norm(Rxy*(dp_end - dp_ee)) > p_th) || (norm(Rxy*dv_ee) > v_th)
            % errors greater than thresholds
            wF_att = dRw' * (A_p*(dp_end - dp_ee) - A_v*(dv_ee));
        end
    
    else % force field disabled
        wF_att = [0,0,0]';
    end
