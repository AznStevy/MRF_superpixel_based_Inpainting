function plot_energy(data_cost, smooth_cost, M, N)
figure(8); clf;
rs = 1; cs = 3;
label_avg = mean(data_cost, 1); label_avg = reshape(label_avg, [M N]);
subplot(rs,cs,1); imagesc(label_avg); title('data single label'); colormap(gray); axis equal tight;
subplot(rs,cs,2); imagesc(data_cost); title('data cost'); colormap(gray);
subplot(rs,cs,3); imagesc(smooth_cost); axis equal tight; title('neigh cost'); axis equal tight;

end

