function plot_curvelet_transforms(I_g_n,I_d_n,I_gs_d_n,I_gs_d_n_var)
figure(3); clf;
rows = 3; cols = 4;
subplot(rows,cols,1);
imagesc(I_g_n{1}); axis equal tight; colormap(gray); title('0 degrees');
subplot(rows,cols,2);
imagesc(I_g_n{2}); axis equal tight; colormap(gray); title('45 degrees');
subplot(rows,cols,3);
imagesc(I_g_n{3}); axis equal tight; colormap(gray); title('90 degrees');
subplot(rows,cols,4);
imagesc(I_g_n{4}); axis equal tight; colormap(gray); title('135 degrees');
subplot(rows,cols,5);
imagesc(I_d_n{1}); axis equal tight; colormap(gray);
subplot(rows,cols,6);
imagesc(I_d_n{2}); axis equal tight; colormap(gray);
subplot(rows,cols,7);
imagesc(I_d_n{3}); axis equal tight; colormap(gray);
subplot(rows,cols,8);
imagesc(I_d_n{4}); axis equal tight; colormap(gray);
subplot(rows,cols,9); % varience
imagesc(I_gs_d_n{1}); axis equal tight; colormap(gray); title(I_gs_d_n_var{1});
subplot(rows,cols,10);
imagesc(I_gs_d_n{2}); axis equal tight; colormap(gray); title(I_gs_d_n_var{2});
subplot(rows,cols,11);
imagesc(I_gs_d_n{3}); axis equal tight; colormap(gray); title(I_gs_d_n_var{3});
subplot(rows,cols,12);
imagesc(I_gs_d_n{4}); axis equal tight; colormap(gray); title(I_gs_d_n_var{4});
end

