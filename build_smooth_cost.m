function [energy] = build_smooth_cost(L_map, Labels, h_image, mask)
    % for image only
    reg = 50; % h_image = rescale(h_image, 0, 255);
    [M, N] = size(mask);
    
    % for patch only
    [x_coords, y_coords] = find(mask == 0);
    s_idx_x = min(x_coords); e_idx_x = max(x_coords);
    s_idx_y = min(y_coords); e_idx_y = max(y_coords);
    M_patch = e_idx_x - s_idx_x; N_patch = e_idx_y - s_idx_y;
    energy = zeros(M_patch*N_patch); % we need a xlen*ylen by xlen*ylen matrix

%     shift_energy = zeros(size(h_image));
%     for i = 1:size(Labels, 1)
%         shift_energy = shift_energy + (h_image - imtranslate(h_image, Labels(i, :))).^2;
%     end
    
    for r = 1:M
        for c = 1:N
            % r,c
            if mask(r, c) == 0
                E_s_r = r - s_idx_x + 1; E_s_c = c - s_idx_y + 1;
                % [M_patch N_patch], E_s_i, E_s_j
              
                try
                    k = sub2ind([M_patch N_patch], E_s_r, E_s_c); % where in the patch
                catch
                    continue;
                end
                
                % boundary conditions
%                 if E_s_i < boundary || E_s_i < boundary || ...
%                         E_s_i > M_patch - boundary || E_s_j > N_patch - boundary
%                     A(sub2ind([M_patch N_patch], E_s_i, E_s_j),k) = 0;
%                     continue;
%                 end

                % define p_1 and l_1 with respect to the IMAGE
                p_1 = [r, c];        
                l_1 = get_offset(r, c, h_image, Labels, L_map); % offset
                if E_s_r > 1
                    rr=r-1; cc=c; E_s_ii = E_s_r - 1; E_s_jj = E_s_c; % up
                    p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);

                    % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2; % 
                    smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
                    energy(k, sub2ind([M_patch N_patch], E_s_ii, E_s_jj)) = smooth_val;
                    % energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
                end
                if E_s_r < M_patch
                    rr=r+1; cc=c; E_s_ii = E_s_r + 1; E_s_jj = E_s_c; % down
                    p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);

                    % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2;
                    smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
                    energy(k, sub2ind([M_patch N_patch], E_s_ii, E_s_jj)) = smooth_val;
                    % energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
                end
                if E_s_c > 1
                    rr=r; cc=c-1; E_s_ii = E_s_r; E_s_jj = E_s_c-1; % left
                    p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);
                    % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2;
                    smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
                    energy(k, sub2ind([M_patch N_patch], E_s_ii, E_s_jj)) = smooth_val;
                    % energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
                end
                if E_s_c < N_patch
                    rr=r; cc=c+1; E_s_ii = E_s_r; E_s_jj = E_s_c+1; % right
                    p_2 = [rr, cc]; % l_2 = get_offset(rr, cc, h_image, Labels, L_map);
                    % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2;
                    smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
                    energy(k, sub2ind([M_patch N_patch], E_s_ii, E_s_jj)) = smooth_val;
                    % energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
                end
            end
               
        end
    end
end

function offset = get_offset(r, c, image, Labels, L_map)
    offset = L_map(r, c, :);
%     if idx_val == 0
%         offset = [0, 0];
%     else
%         offset = Labels(idx_val, :);
%     end
end

% function [energy] = build_smooth_cost_he(L_map, Labels, h_image, mask)
%     % for image only
%     reg = 1; h_image = rescale(h_image, 0, 1);
%     [M, N] = size(h_image);
%     
% 
%     shift_energy = zeros(size(image));
%     for i = 1:size(Labels, 1)
%         shift_energy = shift_energy + imtranslate(h_image, Labels(i, :));
%     end
% 
%     % for patch only
%     [x_coords, y_coords] = find(mask == 0);
%     s_idx_x = min(x_coords); e_idx_x = max(x_coords);
%     s_idx_y = min(y_coords); e_idx_y = max(y_coords);
%     M_patch = e_idx_x - s_idx_x; N_patch = e_idx_y - s_idx_y;
%     energy = zeros(M_patch*N_patch); % we need a xlen*ylen by xlen*ylen matrix
% 
%     for r = 1:M
%         for c = 1:N            
%             if mask(r, c) == 0
%                 E_s_r = r - s_idx_x + 1; E_s_c = c - s_idx_y + 1;
%                 % [M_patch N_patch], E_s_i, E_s_j
%               
%                 try
%                     k = sub2ind([M_patch N_patch], E_s_r, E_s_c); % where in the patch
%                 catch
%                     continue;
%                 end
%                 
%                 % boundary conditions
% %                 if E_s_i < boundary || E_s_i < boundary || ...
% %                         E_s_i > M_patch - boundary || E_s_j > N_patch - boundary
% %                     A(sub2ind([M_patch N_patch], E_s_i, E_s_j),k) = 0;
% %                     continue;
% %                 end
% 
%                 % define p_1 and l_1 with respect to the IMAGE
%                 p_1 = [r, c];        
%                 l_1 = get_offset(r, c, h_image, Labels, L_map); % offset
%                 
%                 if E_s_r > 1
%                     rr=r-1; cc=c; E_s_ii = E_s_r - 1; E_s_jj = E_s_c; % up
%                     p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);
% 
%                     % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2; % 
%                     smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
%                     energy(k, sub2ind([M_patch N_patch], E_s_ii, E_s_jj)) = smooth_val;
%                     energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
%                 end
%                 if E_s_r < M_patch
%                     rr=r+1; cc=c; E_s_ii = E_s_r + 1; E_s_jj = E_s_c; % down
%                     p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);
% 
%                     % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2;
%                     smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
%                     energy(k, sub2ind([M_patch N_patch],E_s_ii,E_s_jj)) = smooth_val;
%                     energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
%                 end
%                 if E_s_c > 1
%                     rr=r; cc=c-1; E_s_ii = E_s_r; E_s_jj = E_s_c-1; % left
%                     p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);
%                     % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2;
%                     smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
%                     energy(k, sub2ind([M_patch N_patch],E_s_ii,E_s_jj)) = smooth_val;
%                     energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
%                 end
%                 if E_s_c < N_patch
%                     rr=r; cc=c+1; E_s_ii = E_s_r; E_s_jj = E_s_c+1; % right
%                     p_2 = [rr, cc]; l_2 = get_offset(rr, cc, h_image, Labels, L_map);
%                     % smooth_val = (shift_energy(rr, cc) - shift_energy(r, c)).^2;
%                     smooth_val = smooth_fn(p_1, p_2, l_1, l_2, h_image, mask) + reg;
%                     energy(k, sub2ind([M_patch N_patch],E_s_ii,E_s_jj)) = smooth_val;
%                     energy(sub2ind([M_patch N_patch], E_s_ii, E_s_jj), k) = smooth_val;
%                 end
%             end
%                
%         end
%     end
% end


% function [A] = build_smooth_cost(E_s, reg)
%     [M,N] = size(E_s{1}); bi_dir = 1;
%     A = zeros(N*M); % we need a xlen*ylen by xlen*ylen matrix
% 
%     % defined by graph cut paper table
%     n_link = @(u, v, n) E_s{n}(u, v);
%     
%     for i = 1:M
%         for j = 1:N
%             
%             k = sub2ind([M N],i,j);
%             if i > 1
%                 ii=i-1; jj=j; % up
%                 A(k,sub2ind([M N],ii,jj)) = n_link(ii,jj,2) + reg;
%                 if bi_dir
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj)); 
%                 end 
%             end
%             if i < M
%                 ii=i+1; jj=j; % down
%                 A(k,sub2ind([M N],ii,jj)) = n_link(ii,jj,1) + reg;
%                 if bi_dir
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                 end
%             end
%             if j > 1
%                 ii=i; jj=j-1; % left
%                 A(k,sub2ind([M N],ii,jj)) = n_link(ii,jj,4) + reg; 
%                 if bi_dir
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                 end
%             end
%             if j < N
%                 ii=i; jj=j+1; % right
%                 A(k,sub2ind([M N],ii,jj)) = n_link(ii,jj,3) + reg;
%                 if bi_dir
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                 end
%             end
%                
%         end
%     end
% end




% function [A] = build_adj_smooth(E_s, Omega)
%     [x_coords, y_coords] = find(Omega == 0);
%     s_idx_x = min(x_coords); e_idx_x = max(x_coords);
%     s_idx_y = min(y_coords); e_idx_y = max(y_coords);
%     N = e_idx_x - s_idx_x; M = e_idx_y - s_idx_y;
%     
%     A = []; % zeros(N*M); % we need a xlen*ylen by xlen*ylen matrix
% 
%     % defined by graph cut paper table
%     n_link = @(u, v, n) E_s{n}(u, v);
% 
%     for i = 1:N
%         for j = 1:M
%             omega_x_coord = i + s_idx_x;
%             omega_y_coord = j + s_idx_y;
%             % st_link = 0;
%             
%             if Omega(omega_x_coord, omega_y_coord) == 0
%                 % fprintf('%d, %d\n', omega_x_coord, omega_y_coord);
%                 
%                 k = sub2ind([M N],i,j);
%                 if i > 1
%                     ii=i-1; jj=j; % up
%                     ii_o=omega_x_coord-1; jj_o=omega_y_coord; % up
%                     A(k,sub2ind([M N],ii,jj)) = n_link(ii_o,jj_o,2);
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                     % st_link = st_link + n_link(ii_o,jj_o,2);
%                 end
%                 if i < N
%                     ii=i+1; jj=j; % down
%                     ii_o=omega_x_coord+1; jj_o=omega_y_coord; % down
%                     A(k,sub2ind([M N],ii,jj)) = n_link(ii_o,jj_o,1);
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                     % st_link = st_link + n_link(ii_o,jj_o,1);
%                 end
%                 if j > 1
%                     ii=i; jj=j-1; % left
%                     ii_o=omega_x_coord; jj_o=omega_y_coord-1; % left
%                     A(k,sub2ind([M N],ii,jj)) = n_link(ii_o,jj_o,4); 
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                     % st_link = st_link + n_link(ii_o,jj_o,4);
%                 end
%                 if j < M
%                     ii=i; jj=j+1; % right
%                     ii_o=omega_x_coord; jj_o=omega_y_coord+1; % right
%                     A(k,sub2ind([M N],ii,jj)) = n_link(ii_o,jj_o,3);
%                     A(sub2ind([M N],ii,jj),k) = A(k,sub2ind([M N],ii,jj));
%                     % st_link = st_link + n_link(ii_o,jj_o,3);
%                 end
%                 
%                 % create labeling node
% 
%                 % connect st links, sum of all 4 edges
%                 % fprintf('%d %d %d\n',s_ind, sub2ind([M N],i,j),st_link);
%                 % fprintf('%d %d %d\n',t_ind, sub2ind([M N],i,j),st_link);
% %                 A(s_ind,sub2ind([M N],i,j)) = st_link; 
% %                 A(sub2ind([M N],i,j),s_ind) = st_link; 
% %                 A(t_ind,sub2ind([M N],i,j)) = st_link; 
% %                 A(sub2ind([M N],i,j),t_ind) = st_link; 
%             end
%         end
%     end
%     
%     
%     %% Disconnect any extraneous nodes
%     % A(s_ind, s_ind) = 0; A(t_ind, t_ind) = 0;
%     % A(t_ind, s_ind) = 0; A(s_ind, t_ind) = 0; % disconnect s,t
% end

