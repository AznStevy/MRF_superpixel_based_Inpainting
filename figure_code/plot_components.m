function plot_components(Omega,O_missing, pOmega, Omega_L, M, F_L, I_d_n, I_gs_n, I_gs_d_n)
    figure(2); clf;
    rs = 3; cs = 3;
    subplot(rs, cs, 1); 
    imagesc(~Omega); axis equal tight; colormap(gray); title('Original Hole');
    subplot(rs, cs, 2); 
    imagesc(O_missing); axis equal tight; colormap(gray); title('Hole');
    subplot(rs, cs, 3); 
    imagesc(pOmega); axis equal tight; colormap(gray); title('\partial\Omega');
    subplot(rs, cs, 4); 
    imagesc(Omega_L); axis equal tight; colormap(gray); title('<\delta'); 
    subplot(rs, cs, 5); 
    imagesc(M); axis equal tight; colormap(gray); title('M'); 
    subplot(rs, cs, 6); 
    imagesc(F_L); axis equal tight; colormap(gray); title('F_L');
    subplot(rs, cs, 7); 
    imagesc(I_d_n); axis equal tight; colormap(gray); title('I_{d}^n');
    subplot(rs, cs, 8); 
    imagesc(I_gs_n); axis equal tight; colormap(gray); title('I_{gs}^n');
    subplot(rs, cs, 9); 
    imagesc(I_gs_d_n); axis equal tight; colormap(gray); title('I_{gs\_d}^n');
end

