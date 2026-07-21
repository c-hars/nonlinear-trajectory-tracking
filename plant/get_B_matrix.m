function B = get_B_matrix(qp, opts)
    arguments
        qp
        opts.enabled = []
    end

    if isempty(opts.enabled)
        enabled = ones(1, qp.n_rotors); % all motors working nominally
        omega_bar = compute_omega_bar(qp, enabled, 'least_squares');
    else
        enabled = opts.enabled;
        omega_bar = compute_omega_bar(qp, enabled, 'keep_opposite_motor_at_nominal_RPM');
    end

    M = compute_M(qp, omega_bar, enabled);

    B = zeros(12, qp.n_rotors);
    B(6,:)  = (1/qp.m)    * M(1,:);  % z_c_ddot = F_t / m
    B(10,:) = (1/qp.I_xx) * M(2,:);  % ω_x_dot  = τ_x / I_xx
    B(11,:) = (1/qp.I_yy) * M(3,:);  % ω_y_dot  = τ_y / I_yy
    B(12,:) = (1/qp.I_zz) * M(4,:);  % ω_z_dot  = τ_z / I_zz
end

function M = compute_M(qp, omega_bar, enabled)
    % Linearised motor mixing matrix about the trim speed omega_bar:
    %   [delta_F_T, delta_tau_x, delta_tau_y, delta_tau_z] ~= M * [delta_omega_1; delta_omega_2; ...; delta_omega_6]
    M = [ 2*qp.kF*omega_bar          .* enabled; ...  % F_t  = Σ kF·ω_i²
          2*qp.kF*omega_bar.*qp.y    .* enabled; ...  % τ_φ  = Σ  y_i·(kF·ω_i²)
         -2*qp.kF*omega_bar.*qp.x    .* enabled; ...  % τ_θ  = Σ –x_i·(kF·ω_i²)
         -2*qp.kM*omega_bar.*qp.dirs .* enabled; ...  % τ_ψ  = Σ (±1)·(kM·ω_i²)
    ];
end