function [Q, R, C, Qy, Qyf, Q_euler, sel] = get_weights(qp, AttRep, opts)
% Recompute the cost matrices for a given attitude representation and sample rate.
% Uses Q_euler (defined via Bryson's rule) as the canonical reference.
% Re-call whenever AttRep or Ts change to ensure a consistent cost function.
%
% Euler-coordinate-centric, DARE terminal cost, scaling for invariance across attitude representations.
%   NB: attitude scaling is first-order equivalence only: exact at identity attitude, degrades with tilt.
%
%   Q       : full-state stage weight, rep-scaled
%   R       : input weight (rep-independent)
%   C       : output selection matrix
%   Qy      : stage output weight, rep-scaled
%   Qyf     : terminal output weight, DARE-based, rep-scaled
%   Q_euler : full-state stage weight in Euler coordinates (physical angles, unscaled)


    arguments
        qp
        AttRep {mustBeMember(AttRep, ["Euler","Quaternion","MRP","FRA", ...
                                      "euler","quaternion","mrp","fra"])} = 'Euler'
        opts.OutputSelection = [1, 2, 3, 9, 10, 11]
    end

    sel = opts.OutputSelection;

    % Q matrix via Bryson (for the Euler coordinate plant)
    max_dx = [[1 1 1]*10*0.01, ...        % | 10 | 5   | [cm]
              [1 1 1]*5*0.01, ...         % | 5  | 5   | [cm/s]
              deg2rad([1 1 1]*5), ...     % | 5  | 2.5 | [deg]
              deg2rad([1 1 1]*50)];       % | 50 | 50  | [deg/s]
    Q_euler = diag(1 ./ max_dx.^2) * qp.Ts;

    % R matrix: based on copter's RPM headroom
    R = diag(1 ./ (qp.max_du).^2) * qp.Ts;
    R = R * 10;

    % Attitude representation scaling
    % Physical angle ~= s * (rep coordinate)
    % Euler/FRA: s=1,  Quaternion: s=2,  MRP: s=4
    switch lower(AttRep)
        case {'euler', 'fra'}, s = 1;
        case 'quaternion',     s = 2;
        case 'mrp',            s = 4;
    end

    % AttRep-scaled Q (attitude block gets scaled by s^2)
    Dx = diag([ones(1,6), s*ones(1,3), ones(1,3)]);
    Q  = Dx * Q_euler * Dx;

    % Output selection
    C = zeros(numel(sel), 12);
    for i = 1:numel(sel), C(i, sel(i)) = 1; end

    is_att = ismember(sel, 7:9);     % attitude states in 12-state layout
    D = diag(1 + (s - 1) * is_att);  % s on attitude outputs, 1 elsewhere

    % Stage output weight
    Qy_eul = Q_euler(sel, sel);
    Qy     = D * Qy_eul * D;

    % Terminal weight: DARE in Euler coords, project, clean, transform
    [A0, B0] = c2d_zoh_expm(get_A_matrix(), get_B_matrix(qp), qp.Ts);
    P   = idare(A0, B0, C' * Qy_eul * C, R);
    Qyf = P(sel, sel);

    % Mask numerical noise (odd tiny negatives from idare)
    mask_thresh = 1e-9 * max(abs(Qyf), [], 'all');
    Qyf(abs(Qyf) < mask_thresh) = 0;

    % Drop cross terms (couplings are small, and this makes the rep transform a pure diagonal scaling)
    Qyf = diag(diag(Qyf)); 

    % Euler -> AttRep coordinates
    Qyf = D * Qyf * D;
end