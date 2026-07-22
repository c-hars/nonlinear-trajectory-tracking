function xe = to_euler_state(x, rep)
% Any-representation state -> Euler-equivalent [p; v_c; euler; omega].

    if iscolumn(x) && (length(x)==12 || length(x)==13)
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

    elseif size(x,2)==12 || size(x,2)==13
        nr = size(x,1);
        xe = zeros(nr,12);
        xe(:,1:6) = x(:,1:6);
        switch rep
            case 'Euler'
                xe = x;
            case "FRA"
                for k=1:nr, xe(k,7:9) = quat2euler(fra2quat(x(k,7:9))); end
                xe(:,10:12) = x(:,10:12);
            case "MRP"
                for k=1:nr, xe(k,7:9) = quat2euler(mrp2quat(x(k,7:9))); end
                xe(:,10:12) = x(:,10:12);
            case "Quaternion"
                for k=1:nr, xe(k,7:9) = quat2euler(x(k,7:10)); end
                xe(:,10:12) = x(:,11:13);
        end

    end

end