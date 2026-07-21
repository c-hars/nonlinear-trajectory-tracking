function xe = to_euler_state(x, rep)
% Any-representation state -> Euler-equivalent [p; v_c; euler; omega].
    switch rep
        case "Euler"
            xe = x;
        case "FRA"
            xe = [x(1:6); quat2euler(fra2quat(x(7:9))); x(10:12)];
        case "MRP"
            xe = [x(1:6); quat2euler(mrp2quat(x(7:9))); x(10:12)];
        case "Quaternion"
            q  = x(7:10) / norm(x(7:10));
            xe = [x(1:6); quat2euler(q); x(11:13)];
    end
end