function q = mrp2quat(p)
    ps = p'*p;
    q  = [ (1 - ps); 2*p ] / (1 + ps);
end