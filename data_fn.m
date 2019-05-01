function data_val = data_fn(point, L_map, mask)
    max_val = 1000000; % basically infinity
    
    r = point(1); c = point(2);
    try
        if isequal(L_map(r, c), [0,0])
            data_val = max_val; return;
        end

        new_r = r + L_map(r, c, 1);
        new_c = c + L_map(r, c, 2);


        if mask(new_r, new_c) == 1
            data_val = 0; return;
        else
            data_val = max_val; return;
        end
    catch
        data_val = max_val; return;
    end
end


% function data_val = data_fn(point, L_map, mask)
%     max_val = 1000000; % basically infinity
%     r = point(1); c = point(2);
%     try
% 
%        % see if we are in the mask
%        if ~is_hole(r, c, mask)
%            if L_map(r,c,1) == 0 && L_map(r,c,2) == 0
%               data_val = 0; return;
%            else
%               data_val = max_val; return;
%            end
%        end
% 
%        % test to see if new is valid
%        new_r = r + L_map(r, c, 1);
%        new_c = c + L_map(r, c, 2);
% 
%        if ~is_hole(new_r, new_c, mask)
%           data_val = 0; return;
%        end
%    catch
%         data_val = max_val; return;
%    end
%    data_val = max_val; return;
% end

% function data_val = data_fn(point, L_map, mask)
%     max_val = 1000000; % basically infinity
%     r = point(1); c = point(2);
%     try
% 
%        % see if we are in the mask
%        if is_hole(r, c, mask)
%            if L_map(r,c,1) == 0 && L_map(r,c,2) == 0
%               data_val = 0; return;
%            else
%               data_val = max_val; return;
%            end
%        end
% 
%        % test to see if new is valid
%        new_r = r + L_map(r, c, 1);
%        new_c = c + L_map(r, c, 2);
% 
%        if is_hole(r, c, mask)
%           data_val = 0; return;
%        end
%    catch
%         data_val = max_val; return;
%    end
%    data_val = max_val; return;
% end



