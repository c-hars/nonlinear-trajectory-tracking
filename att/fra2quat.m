function q = fra2quat(th)
    th = th(:);
    ang = norm(th);
    if ang < 1e-8
        q = [1; th/2];
    else
        q = [cos(ang/2); sin(ang/2)*(th/ang)];
    end
    q = q / norm(q);
end