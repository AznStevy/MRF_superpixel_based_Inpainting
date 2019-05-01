function S = get_offset(NNF, image)
	new_NNF = cat(3, NNF(:, :, 2), NNF(:, :, 1), NNF(:,:,3)); 
    NNF = new_NNF;  
    
    S = zeros(size(NNF));

    for u = 1:size(image, 1)
        for v = 1:size(image, 2)
            S(u, v, 1) = NNF(u, v, 1)-u; 
            S(u, v, 2) = NNF(u, v, 2)-v;
        end
    end
end
