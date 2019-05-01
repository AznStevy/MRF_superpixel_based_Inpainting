function plot_ilabels(L_mask, L_map)
figure(7); clf;
rs = 1; cs = 3; 
subplot(rs, cs, 1);
imagesc(L_mask); 
axis equal tight; colormap(gray); title('Labeling'); colorbar;
subplot(rs, cs, 2);
imagesc(L_map(:, :, 1)); 
axis equal tight; colormap(gray); title('Label Offset X');
subplot(rs, cs, 3);
imagesc(L_map(:, :, 2)); 
axis equal tight; colormap(gray); title('Label Offset Y');
end

