function rec_image = get_match(image, S)
    % S is an offset matrix
    rec_image = zeros(size(image));
    for u = 1:size(image, 1)
        for v = 1:size(image, 2)
            try
                % new_image(u, v, :) = image(NNF(u, v, 1), NNF(u, v, 2), :);
                rec_image(u, v, :) = image(S(u,v,1) + u, S(u,v,2) + v, :);
            catch; 
            % disp('error');
            end;
        end
    end

end

