function plot_inputs(og_image, mask)
figure(1); clf; rs = 1; cs = 2;
subplot(rs, cs, 1);
imagesc(og_image); axis equal tight;
subplot(rs, cs, 2);
imagesc(mask); axis equal tight; title('Mask'); colormap(gray);
end

