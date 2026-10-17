function A_fcn = get_SDC_A_function(rep)
    switch lower(rep)
        case 'euler',      A_fcn = @(xk,qp) get_A_matrix_SDRE_EulerAttitude(xk,qp);
        case 'mrp',        A_fcn = @(xk,qp) get_A_matrix_SDRE_MRPAttitude(xk,qp);
        case 'fra',        A_fcn = @(xk,qp) get_A_matrix_SDRE_FRAAttitude(xk,qp);
        case 'quaternion', A_fcn = @(xk,qp) get_A_matrix_SDRE_QuaternionAttitude(xk,qp);
    end
end