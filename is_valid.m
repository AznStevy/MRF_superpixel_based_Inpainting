function valid = is_valid(x, y, Omega)
    valid = true;
    [M, N] = size(Omega);
    if (x > 0 && y > 0 && x <= M && y <= N)
       valid = true; 
    else
       valid = false; 
    end
end

