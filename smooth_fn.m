% function smooth_val = smooth_fn(p_1, p_2, l_1, l_2, image, mask)
%     image = rescale(image, 0, 1);
%     if (l_1 == l_2)
%         return; 
%     end
%      
%     %% if all of these are ok points
%  
%     % value of the image, might be RGB or gray
%     a = imtranslate(image, l_1);
%     b = imtranslate(image, l_2);
% 
%     smooth_val = sum((a(:)-b(:)).^2);
% end

function smooth_val = smooth_fn(x_1, x_2, s_a, s_b, image, mask)
    max_val = 1000000; % basically infinity
    smooth_val = 0;
    if isequal(s_a, s_b)
        return; 
    end
    
    x1_s_a = s_a + x_1;
    x1_s_b = s_b + x_1;
    x2_s_a = s_a + x_2;
    x2_s_b = s_b + x_2;
    
    %% if all of these are ok points
    if is_valid(x1_s_a(1), x1_s_a(2), mask) && ...
            is_valid(x1_s_b(1), x1_s_b(2), mask) && ...
            is_valid(x2_s_a(1), x2_s_a(2), mask) && ...
            is_valid(x2_s_b(1), x2_s_b(2), mask)
        % disp('Exists');
        
        % value of the image, might be RGB or gray
        v1_a = image(x1_s_a(1), x1_s_a(2));
        v1_b = image(x1_s_b(1), x1_s_b(2));
        v2_a = image(x2_s_a(1), x2_s_a(2));
        v2_b = image(x2_s_b(1), x2_s_b(2));
        
        smooth_val = smooth_val + ...
            (norm(v1_a(:)-v1_b(:)).^2 + norm(v2_a(:)-v2_b(:)).^2);
        % smooth_val = smooth_val;
        
        return;
    end
    smooth_val = max_val;
end

% function smooth_val = smooth_fn(p_1, p_2, l_1, l_2, image, mask)
%     max_val = 1000000; % basically infinity
%     smooth_val = 0;
%     if (l_1 == l_2)
%         return; 
%     end
%     
%     x1_s_a = l_1 + p_1;
%     x1_s_b = l_2 + p_1;
%     x2_s_a = l_1 + p_2;
%     x2_s_b = l_2 + p_2;
%     
%     %% if all of these are ok points
%     if is_valid(x1_s_a(1), x1_s_a(2), mask) && ...
%             is_valid(x1_s_b(1), x1_s_b(2), mask) && ...
%             is_valid(x2_s_a(1), x2_s_a(2), mask) && ...
%             is_valid(x2_s_b(1), x2_s_b(2), mask)
%         
%         % value of the image, might be RGB or gray
%         v1_a = image(x1_s_a(1), x1_s_a(2));
%         v1_b = image(x1_s_b(1), x1_s_b(2));
%         v2_a = image(x2_s_a(1), x2_s_a(2));
%         v2_b = image(x2_s_b(1), x2_s_b(2));
%         
%         for i = 1:length(v1_a(:))
%             smooth_val = smooth_val + ...
%                 (v1_a(i)-v1_b(i)).^2 + ...
%                 (v2_a(i)-v2_b(i)).^2;
%         end
%         
%         return;
%     end
%     smooth_val = max_val;
% end

