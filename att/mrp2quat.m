function q = mrp2quat(p)
    p = p(:);
    ps = p'*p;
    q  = [ (1 - ps); 2*p ] / (1 + ps);
end