function omega_bar = compute_omega_bar(qp, enabled, method)
% Compute omega_bar given all motors are working OR motor 1 is disabled etc.

    kF = qp.kF;
    kM = qp.kM;
    g  = qp.g;

    %   Nonlinear motor mixing matrix: [F_T, tau_x, tau_y, tau_z] = M_mix_NL * [omega_1^2; omega_2^2; ...; omega_6^2]
    M_mix_NL = [ ...
        kF * enabled; ...
        kF * enabled .* qp.y; ...
       -kF * enabled .* qp.x; ...
       -kM * enabled .* qp.dirs; ...
    ];
    
    if strcmp(method, 'least_squares')
        % Solve for the omega_bar s.t. [F_T, tau_x, tau_y, tau_z] = [mg, 0, 0, 0]
        %  - multiple solutions possible
        %  - use of pinv() selects the minimum Euclidean norm soln
        %  - sometimes pinv() gives solutions with trivial imaginary components - safely ignore
        req = [qp.m*g; 0; 0; 0];
        W = pinv(M_mix_NL) * req;
        omega_bar = real(sqrt(W))'; % rad/s

    elseif strcmp(method, 'keep_opposite_motor_at_nominal_RPM')
        
        % This is hexacopter-specific
        % Solve for the omega_bar s.t. [F_T, tau_x, tau_y, tau_z] = [mg, 0, 0, 0]
        %  - multiple solutions possible
        %  - use of pinv() selects the minimum Euclidean norm soln
        %       - consequence: if M1 is offline, then the opposite motor M4 is offline - the equilibrium is something like omega_bar = [0 100 100 0 100 100] rad/s
        %       - in this approach, we mandate that M4 is not taken offline - equilibrium is something like [0 ? ? w_4 ? ?] rad/s

        % Find the disabled motor index, and the motor opposite it.
        disabled_motor_idx = find(enabled ~= 1);
        
        if isempty(disabled_motor_idx)
            % disp("keep_opposite_motor_at_nominal_RPM method used despite no motor failure - reverting to least_squares")
            omega_bar = compute_omega_bar(qp, enabled, 'least_squares');
            return
        elseif numel(disabled_motor_idx) > 1
            warning("Invalid usage: keep_opposite_motor_at_nominal_RPM assumes only one motor is offline. Reverting to least squares.")
            omega_bar = compute_omega_bar(qp, enabled, 'least_squares');
            return
        end

        opposite_motor_idx = mod(disabled_motor_idx+2,6)+1;
        i = opposite_motor_idx;
        
        % Design assuming Motor i (opposite to the offline motor) maintains its normal hover RPM.
        nominal_omega_bar = compute_omega_bar(qp, [1 1 1 1 1 1], 'least_squares');
        omega_i_bar = nominal_omega_bar(opposite_motor_idx);

        % What force/torques does Motor i yield at this RPM?
        F_i =  qp.kF               * omega_i_bar^2;
        Tx_i = qp.kF *  qp.y(i)    * omega_i_bar^2;
        Ty_i = qp.kF * -qp.x(i)    * omega_i_bar^2;
        Tz_i = qp.kM * -qp.dirs(i) * omega_i_bar^2;
        
        % For equilibrium, the force and torques must be [F_t; tau_x; tau_y; tau_z] = [qp.m*9.81; 0; 0; 0]
        % With Motor i prescribed, the net force and torque required from the other 5 motors is:
        req = [qp.m*9.81;0;0;0] - [F_i;Tx_i;Ty_i;Tz_i];

        % Find the least-squares solution to this problem
        to_solve_for = 1:6;
        to_solve_for(i) = [];

        M_mix_NL = M_mix_NL(:,to_solve_for);
        W = pinv(M_mix_NL)*req;
        omega_bar_solved_for = real(sqrt(W));
        
        % Concatenate solved for RPMs + prescribed RPM
        omega_bar = zeros(1,6);
        omega_bar(i) = omega_i_bar;
        omega_bar(to_solve_for) = omega_bar_solved_for;

    end

end 
