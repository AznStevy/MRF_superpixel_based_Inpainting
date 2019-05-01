function plot_result(L_map_op, L_optimized_x,L_optimized_y,image,new_image,mask)
figure(9); clf;
rs = 2; cs = 3; % down, up, right, left
subplot(rs, cs, 1); imagesc(L_map_op); axis equal tight; colormap(gray); title('Label Map');
subplot(rs, cs, 2); imagesc(L_optimized_x); axis equal tight; colormap(gray); title('Shift row');
subplot(rs, cs, 3); imagesc(L_optimized_y); axis equal tight; colormap(gray); title('Shift col');
subplot(rs, cs, 5); imagesc(image); axis equal tight; colormap(gray); title('Original Image');
subplot(rs, cs, 6); imagesc(new_image); axis equal tight; colormap(gray); title('Infilled image');
hold on; visboundaries(~mask); hold off; %// draw rectangle on image; %// draw rectangle on image

figure(10); clf; psnr_ = 20*log(psnr(new_image, image)); ssim_ = ssim(new_image, image);
imagesc(new_image); axis equal tight; colormap(gray); 
title(sprintf('PSNR: %f dB, SSIM: %f',psnr_, ssim_));
hold on; visboundaries(~mask); hold off; %// draw rectangle on image; %// draw rectangle on image
end

