function plot_variance(image,I_ns,struct_union,I_s_n)
figure(4); clf;
rs = 2; cs = 4;
subplot(rs, cs, 1);
imagesc(image); axis equal tight; title('Original');
subplot(rs, cs, 2); 
imagesc(I_ns); axis equal tight; colormap(gray); title('I_{ns}');
subplot(rs, cs, 3); 
imagesc(struct_union); axis equal tight; colormap(gray); title('Union');
subplot(rs, cs, 5); 
imagesc(I_s_n{1}); axis equal tight; colormap(gray); title('I_s^n_{1}');
subplot(rs, cs, 6); 
imagesc(I_s_n{2}); axis equal tight; colormap(gray); title('I_s^n_{2}');
subplot(rs, cs, 7); 
imagesc(I_s_n{3}); axis equal tight; colormap(gray); title('I_s^n_{3}');
subplot(rs, cs, 8); 
imagesc(I_s_n{4}); axis equal tight; colormap(gray); title('I_s^n_{4}');
end

