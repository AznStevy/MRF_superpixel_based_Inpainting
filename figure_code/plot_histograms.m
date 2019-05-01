function plot_histograms(dds_indices, h_top)
    figure(6); clf;
    rs = 1; cs = length(dds_indices); i = 1;
    for idx = dds_indices
        % determine appropriate bounds
        rel_idx = find(h_top{idx}(:) > 0); % single number
        first_idx = rel_idx(1); last_idx = rel_idx(end);
        [f_r,f_c] = ind2sub(size(h_top{idx}), first_idx);
        [l_r,l_c] = ind2sub(size(h_top{idx}), last_idx);
        sqr_size = max(l_r - f_r, l_c - f_c);
        
        subplot(rs, cs, i);
        try
        imagesc(h_top{idx}(f_r:f_r + sqr_size, f_c:f_c + sqr_size));
        catch
            imagesc(h_top{idx});
        end
        axis equal tight; title(idx); colormap(parula(256)); colorbar; 
        i = i + 1;
        
    end
end

