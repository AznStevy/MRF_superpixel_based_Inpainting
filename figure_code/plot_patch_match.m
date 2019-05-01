function plot_patch_match(image,dds_indices,I_ns,S_n)
figure(5); clf;
rs = 3; cs = length(dds_indices) ; i = 1; % for making plots
for idx = dds_indices
    patch = get_match(image, S_n{idx});
    subplot(rs, cs, i); imagesc(patch); 
    axis equal tight; colormap(gray); title(sprintf('Patchmatch idx=%d', idx))

    val_range = 1:230;
    subplot(rs, cs, i + cs); imagesc(S_n{idx}(val_range,val_range,1)); 
    axis equal tight; title(sprintf('x_{%d}', idx))

    subplot(rs, cs, i + 2*cs); imagesc(S_n{idx}(val_range,val_range,2)); 
    axis equal tight; title(sprintf('y_{%d}', idx))
    
    i = i + 1;
end
end

