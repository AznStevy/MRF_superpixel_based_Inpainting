function Z_n = find_Z(Q, n)
    % define the angles
    if n == 1
        s_a = 337.5; e_a = 22.5;
        i_1 = {1, 2:3, 3:6, 3:6, 5:12};
        i_2 = {1, 10:11, 19:22, 19:22, 37:44};
    elseif n == 2
        s_a = 22.5; e_a = 67.5;
        i_1 = {1, 7:8, 15:18, 15:18, 29:36};
        i_2 = {1, [1 12], [31, 32, 1 2], [31, 32, 1 2], ...
            [61:64, 1:4]};
    elseif n == 3
        s_a = 67.5; e_a = 112.5;
        i_1 = {1, 6:7, 11:14, 11:14, 21:28};
        i_2 = {1, 14:15, 27:30, 27:30, 53:60};
    elseif n == 4    
        s_a = 112.5; e_a = 157.5;
        i_1 = {1, 4:5, 7:10, 7:10, 5:12};
        i_2 = {1, 12:13, 23:26, 23:26, 45:52};
    end
    
    Z_n = [];
    for scale = 2:5
        indices = i_1{scale};
        coefs = Q{scale}{indices};
        Z_n = [Z_n coefs(:)']; % put it on the side
        
%         indices = i_2{scale};
%         Q{scale}
%         coefs = Q{scale}{indices};
%         Z_n = [Z_n coefs(:)']; % put it on the side
    end
end

