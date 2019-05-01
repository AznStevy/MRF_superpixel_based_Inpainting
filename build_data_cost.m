% function data_cost = build_data_cost(L_map, mask)
%     [M, N] = size(E_d); MN = M*N;
%     data_cost = ones(MN) * 10000; % max cost for everything else
%     
%     for k = 1:MN
%        data_cost(1, k) = k;
%        data_cost(2, k) = data_fn([u, v], L_map, mask);
%     end
% end

function data_cost = build_data_cost(Label_map, Labels, mask)
    K = size(Labels, 2);
    max_value = 1E6;
  
    [x_coords, y_coords] = find(mask == 0);
    s_idx_x = min(x_coords); e_idx_x = max(x_coords);
    s_idx_y = min(y_coords); e_idx_y = max(y_coords);
    M_patch = e_idx_x - s_idx_x; N_patch = e_idx_y - s_idx_y;
    data_cost = zeros(K, M_patch*N_patch);
    
    for u = 1:size(mask, 1)
        for v = 1:size(mask, 2)
            if mask(u, v) == 0 
                
                E_d_r = u - s_idx_x + 1; % on the mask, first is 1.
                E_d_c = v - s_idx_y + 1;

                try; k = sub2ind([M_patch, N_patch], E_d_r, E_d_c); catch; end;
                
                shifted_r = u + Label_map(u, v, 1);
                shifted_c = v + Label_map(u, v, 2);
                
                [idx, l_idx] = ismember([Label_map(u, v, 1), Label_map(u, v, 2)], ...
                    Labels, 'rows');

                if is_valid(shifted_r, shifted_c, mask) && ...
                            mask(shifted_r, shifted_c) ~= 0
                    data_cost(l_idx, k) = 1; % max_value;
                    % fprintf('%d, %d \n',);
                end    
            end
        end
    end
    
    data_cost = (1 - data_cost)*max_value;
end