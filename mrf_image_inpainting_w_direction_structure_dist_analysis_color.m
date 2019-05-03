% Implementation of the algorithm developped by Jixiang Cheng and Zhidan
% Li for Markov random field-based image inpainting with direction
% structure distribution analysis for maintaining structure coherence.

clear all; close all;
%% Path
addpath(genpath(pwd));

%% Tunable parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
image_name = 'lena'; % must be 256 by 256 (PM limit), images folder, png
alpha = 2; beta = .2; % (12)
delta = 10; % implementation by paper 
super_pixel = true;
boundary = 2;

% Put hole in image; row col width height
hole = [80 80 100 100]; % [50 160 180 70] [80 80 100 100]; [20 100 50 50]

%% Load the image
im_scale = 1;
image = imread(sprintf('%s.png', image_name)); 
image = im2double(image);
image = imresize(image, im_scale); og_image = image; image = rescale(image, 0, 1);
image_gray = image; try; image_gray = rgb2gray(image_gray); catch; end;
try; image= gray2rgb(image); catch; end;
r_image = size(image, 1); c_image = size(image, 2);


%% Start Process
mask = ones(size(image_gray));
mask_og = ones(size(mask)); % define the hole
mask(hole(1)-boundary: hole(1)+hole(3)+boundary, ...
    hole(2)-boundary:hole(2)+hole(4)+boundary) = 0;
mask_og(hole(1):hole(1)+hole(3), hole(2):hole(2)+hole(4)) = 0;
h_image = mask_og .* image;

%% plot inputs
plot_inputs(og_image, mask);

%% This is to determine the DDS
% We input an image I_Y and we perform a forward transform and obtain
% the coefficient matrix.
tic;
% This operation here is what allows us to recover structural information
% from the image.
% Tuning parameters
finest = 1; % curvelets at finest scale
nbscales = log2(r_image) - 3; % evals to 5
nbangles_coarse = 8; %8;
h_image_gray = h_image; try; h_image_gray = rgb2gray(h_image); catch; end;
Q = fdct_wrapping(image_gray,0,finest,nbscales,nbangles_coarse); % (3)
Z = {};
A = {};
for n = 1:4
    Z{n} = find_Z(Q, n);
    % we need to modify Q according to (5)
    H_n = Q; % varies with n
    for s = 1:size(Q, 2)
        for d = 1:size(Q{s}, 2)
            % for the coeffecients in Q, check they are in Z_n
            coefs = Q{s}{d};
            member_index = ~ismember(coefs, Z{n});
            coefs(member_index) = 0;
            H_n{s}{d} = coefs;
        end
    end
    A{n} = real(ifdct_wrapping(H_n, 0, r_image, c_image));
end
disp('Determined coefficient matrix.'); toc; 

%% 3.1 Do edge detection and evaluate the gradients of the image.
% for each of the values of A, obtain an x and y gradient.
tic;
I_g_n = {}; I_d_n = {}; I_gs_n = {};
I_gs_d_n = {}; I_gs_d_n_var = {};
grad_op = [-0.5, 0, 0.5];

for n = 1:4

    G_x = conv2(A{n}, grad_op, 'same');
    G_y = conv2(A{n}, grad_op', 'same');
    
    %% for 3.2 Structure feature selection
    I_g_n{n} = sqrt((G_x).^2 + (G_y).^2); % full image gradients
    % canny edge and expansion operators
    im_edge = round(edge(A{n}, 'Canny'));
    % dilate the edges slightly
    se = strel('ball', 2, 2);
    I_d_n{n} = imdilate(im_edge, se);
    I_d_n{n} = rescale(I_d_n{n}, 0, 1) > 0.1; % normalize
    % Find I_gs_n
    O_missing = edge(mask, 'Sobel'); %~Omega; % ~edge(Omega, 'Canny'); % the region edge
    
    % find shortest distance from point to edge of Omega, done by bwdist
    pOmega = bwdist(O_missing); 
    Omega_L = pOmega < delta; 
    M = mask ~= 1; % generate mask
    F_L = M .* Omega_L;% good feature.
    I_gs_n{n} = I_g_n{n} .* F_L;
    I_gs_d_n{n} = I_gs_n{n} .* I_d_n{n};
    I_gs_d_n_var{n} = var(I_gs_d_n{n}(:));
    
    %% Plot each of those
    plot_components(mask,O_missing, pOmega, ...
        Omega_L, M, F_L, I_d_n{n}, I_gs_n{n}, I_gs_d_n{n});
end
disp('Selected features.'); toc;

%% display the gradient magnitudes
plot_curvelet_transforms(I_g_n,I_d_n,I_gs_d_n,I_gs_d_n_var);

%% Selection based on varience
tic;
I_gs_d_n_var = cell2mat(I_gs_d_n_var);
min_var = min(I_gs_d_n_var);
max_var = max(I_gs_d_n_var);
alphas = I_gs_d_n_var./min_var > alpha;
betas = I_gs_d_n_var./max_var > beta;

intersect_a_b = alphas .* betas;
dds_indices = find(intersect_a_b==1); % good indices
num_dds = length(dds_indices);
disp('Selected features based on varience.'); toc;

%% 3.3 calculate offsets in selected direction structure 
% and non-structure image for the good DDS
I_s_n = {}; S_n = {}; i = 1; % compute patch region
% initialize for figures
for n = 1:4
    I_s_n{n} = zeros(size(image));
end

% for rerunning
dds_indices = unique(dds_indices(dds_indices ~= 5));
for idx = 1:4
    % direction edge image
    I_s_n{idx} = double(image).*I_d_n{idx}; % inverse
    I_s_n{idx} = I_s_n{idx}./max(I_s_n{idx}(:));
    
    % we want to invert this because important information is 1, other 0
    if idx == 1
       struct_union = I_s_n{idx};
    end
    struct_union = max(struct_union, I_s_n{idx}); 
end
I_ns = image./struct_union; 
I_ns = 1-I_ns./max(I_ns(:)); 
I_ns(isnan(I_ns))=0;

%% Show the intersection and I_ns
plot_variance(og_image, I_ns, struct_union, I_s_n);

%% Calculate using patch match
% uses the intersections to align structures on edges.

patchsize = 9; iters = 5; cores = 2; tau = 8;

dds_indices = dds_indices(dds_indices ~= 5);
num_dds = length(dds_indices);

bw_dist = bwdist(~mask); 
valid_mask = bw_dist > tau; valid_mask = ~valid_mask;
valid_mask = double(gray2rgb(valid_mask) > 0);

%% create superpixel boundary
if super_pixel
    %% Super pixels
    N = 200;
    [L, NumLabels] = superpixels(image,N);
    N = NumLabels;

    BW = boundarymask(L);
    idx = label2idx(L);
    
    %% Figure for boundaries
    boundary_mask = zeros(size(image),'like',image);
    for labelVal = 1:N
        Idx = idx{labelVal};
        if any(h_image(Idx) == 0)
            boundary_mask(Idx) = 1;
        end
    end
    
    valid_mask = valid_mask .* boundary_mask;
end

%% get NNF

image = rescale(image, 0, 1);

tic; NNF = cell(1,5);
for idx = dds_indices
    % compute the offset according to for n = 1 to 4
    % this requires another steepest descent algorithm using expectation
    % maximization of nearest neighbors.
    % for every part of the image, compute the patch match of it.

    % nnmex(1A, 2B, 3[algo='cpu'], 4[patch_w=7], 5[nn_iters=5], 
    % 6[rs_max=100000], 7[rs_min=1], 8[rs_ratio=0.5], 9[rs_iters=1.0], 10[cores=2], 11[bmask=NULL]...
    % 12[win_size=[INT_MAX INT_MAX]], 13[nnfield_prev=[]], 14[nnfield_prior=[]], 15[prior_winsize=[]], 16[knn=1], 17[scalerange=4])
    temp = I_s_n{idx}; try; temp = gray2rgb(temp); catch; end; % temp = image
    temp_2 = I_s_n{idx}; try; temp_2 = gray2rgb(temp_2); catch; end; % temp = image
    NNF{idx} = nnmex(temp, temp_2, 'cputiled', patchsize, iters, [], [], [], [], cores, valid_mask);
    % image, h_image
    fprintf('Finished PatchMatch for n=%d ', idx); toc;
end
% according to the section, n is 1-5 where 5th uses I_ns.
% struct_union
temp = struct_union; try; temp = gray2rgb(temp); catch; end;
temp_2 = struct_union; try; temp_2 = gray2rgb(temp_2); catch; end; % temp = image
NNF{5} = nnmex(temp, temp_2, 'cputiled', patchsize, iters, [], [], [], [], cores, valid_mask);
fprintf('Finished PatchMatch for n=5 '); toc;

%% Convert to offsets + distance constraint
dds_indices = unique([dds_indices 5]);
num_dds = length(dds_indices);

tic; disp('Converting to offsets.');
create_str_mask = @(image) (rescale(im2bw(image), 0 ,1) >.5);
for idx = dds_indices
    S_n{idx} = get_offset(NNF{idx}, image);
    if idx == 5
        structure = create_str_mask(~mask .* I_ns);    
    else
        structure = create_str_mask(~mask .* I_s_n{idx});
    end
    S_n{idx} = S_n{idx} .* structure; % (rescale(im2bw(I_gs_n{idx}), 0 ,1) >.5);
end
toc;

%% Visualize Patch match results
plot_patch_match(image, dds_indices, I_ns, S_n);

%% calculate 2-d direction histogram (14)
disp('Starting histogram calculation.');
tic;

r2_image = 2*r_image; c2_image = 2*c_image;
offsets = {}; h = {};
% valid_dds_offsets = S_n;
for idx = dds_indices
    offsets{idx} = [];
    for u = 1:r_image
        for v = 1:c_image
            if F_L(u, v) == 1 || mask(u, v) == 0
                if S_n{idx}(u,v,1) ~= 0 || S_n{idx}(u,v,2) ~= 0
                    new_point = [S_n{idx}(u,v,1) S_n{idx}(u,v,2)];
                    offsets{idx} = [offsets{idx}; new_point];
                end
            end
        end 
    end
end

for idx = dds_indices
    h{idx} = zeros(r2_image, c2_image);
    for i = 1:size(offsets{idx},1)
        o_set = offsets{idx}(i,:); % first offset pair
        vert_off = o_set(1) + r_image;
        horz_off = o_set(2) + c_image;
        h{idx}(vert_off, horz_off) = h{idx}(vert_off, horz_off) + 1;
    end
end

disp('Histograms computed.'); toc;

%% select number of peaks, based on size of DDS
tic;

num_S_1_4 = num_dds - 1; % because num_dds includes the 5th.
if num_S_1_4 == 1
    K1 = 40; K2 = 20;
elseif num_S_1_4 == 2
    K1 = 20; K2 = 20;
elseif num_S_1_4 == 3
    K1 = 15; K2 = 15;
elseif num_S_1_4 == 4
    K1 = 10; K2 = 20;        
else
    K2 = 20;
end
% if super_pixel
%     K1 = 200; K2 = 200;
% else
%     K1 = 200; K2 = 200;
% end

h_top = {};
% Pick out K highest peak values for all of them
for idx = dds_indices
    sorted_h = sort(h{idx}(:), 'descend');
    if idx ~= 5
        sorted_h = sorted_h(1:K1);
    else
        sorted_h = sorted_h(1:K2);
    end
    h_top{idx} = h{idx}.*ismember(h{idx}, sorted_h);
end

% pick out the top K offsets
top_K_offsets = [];
for idx = dds_indices
    for u = 1:r2_image
        for v = 1:c2_image
            if h_top{idx}(u, v) ~= 0
                x_offset = u-r_image; y_offset = v-c_image;
                top_K_offsets = [top_K_offsets; x_offset, y_offset];
            end
        end
    end
end

disp('Top offsets from histograms picked.'); toc;


%% Display top histograms + labels in S
plot_histograms(dds_indices, h_top);

%% Compute labeling map L for Omega
tic;
L_map = zeros(r_image, c_image, 2);
L_mask = zeros(r_image, c_image);
% just to stay consistent with the other file
for idx = dds_indices
    disp(sprintf('Created labeling mask %d.', idx));
    for u = 1:r_image
        for v = 1:c_image
            valid = 0;
            if mask_og(u, v) == 0
                for i = 1:size(top_K_offsets, 1)
                    if S_n{idx}(u, v, 1) == top_K_offsets(i, 1) ...
                            && S_n{idx}(u, v, 2) == top_K_offsets(i, 2)
                        L_mask(u, v) = i;
                        L_map(u, v, 1) = S_n{idx}(u, v, 1);
                        L_map(u, v, 2) = S_n{idx}(u, v, 2);
                        break;
                    end
                end
            end
        end
    end
end

% compute the top labeling masks
disp('Created labeling masks.'); toc;

% based on He's work https://ieeexplore.ieee.org/document/6853394
% With this step, we have computed the labeling map L(x) which contains the
% necessary offsets to optimize in the next step.

%% view histogram and L
plot_ilabels(L_mask, L_map);

%% Optimize the energy function using offsets.
% This is the part that actually does the infilling by optimizing the
% Markov Random Field energy function (which is the labeling/MAP function)
% Here, we use the same graph-based image completion method using the
% statistics of patch offsets.
% http://imagine.enpc.fr/~marletr/enseignement/mathimage/lecture5_slides.pdf
% http://vision.stanford.edu/teaching/cs231a_autumn1112/lecture/lecture6_clustering_and_seg_p2_cs231a.pdf
% What I want is the energy distribution of omega...

%% Use alpha expansion defined by the original paper.
% This was implemented by Olga Veksler and Andrew Delong and is available
% at https://vision.cs.uwaterloo.ca/code/
[r_coords, c_coords] = find(mask == 0);
s_idx_r = min(r_coords); e_idx_r = max(r_coords);
s_idx_c = min(c_coords); e_idx_c = max(c_coords);
M = e_idx_r - s_idx_r; N = e_idx_c - s_idx_c; MN = M*N;
tic; disp('Calculating energies');
% create data cost
neighbors = build_grid(M, N); % build a NM by NM matrix
data_cost = build_data_cost(L_map, top_K_offsets, mask); % MN x K
smooth_cost = build_smooth_cost(L_map, top_K_offsets, h_image, mask); % compute smoothing term
disp('Finished Calculating energies.'); toc;

%% display energies
plot_energy(data_cost, smooth_cost, M, N);

%% Do multilabel graph cuts
tic; disp('Starting graph cut optimization algorithm.');

K = size(data_cost, 1); % For graph optimization
H = GCO_Create(MN, K); % GCO_create(sites, labels)
GCO_SetDataCost(H, data_cost); % unary matrix
GCO_SetNeighbors(H, smooth_cost); %smooth_cost);
GCO_Expansion(H);

labeling = GCO_GetLabeling(H);
[E_a_op E_d_op E_s_op l_op] = GCO_ComputeEnergy(H);
GCO_Delete(H);

% GCO_Delete(H);
disp('Completed graph cut optimization algorithm.'); toc;

%% convert label to offset
tic;
raw_offsets = reshape(labeling, [M N]);

L_optimized = zeros([M N, 2]);
L_map_optimized = zeros([M N]);
for k = 1:size(labeling)
    label = labeling(k);
    [x, y] = ind2sub([M N],k); % position
    
    L_map_optimized(x, y) = label;
    L_optimized(x, y, 1) = top_K_offsets(label, 1);
    L_optimized(x, y, 2) = top_K_offsets(label, 2);
end

% interpolate it by a factor of scale
L_optimized_x = L_optimized(:,:,1);
L_optimized_y = L_optimized(:,:,2);
L_optimized = cat(3, L_optimized_x, L_optimized_y);
disp('Converted to offsets.'); toc;

%% Fill in image
tic;
new_image = image;
for u = 1:r_image
    for v = 1:c_image
        mask_x = u - s_idx_r + 1; mask_y = v - s_idx_c + 1;
        if mask_x < M && mask_y < N && mask_x > 0 && mask_y > 0
            try
                if mask_og(u,v) == 0
                    offset_x = L_optimized(mask_x, mask_y, 1);
                    offset_y = L_optimized(mask_x, mask_y, 2);
                    new_image(u, v, :) = image(u + offset_x, v+offset_y, :);
                end
            catch d
                new_image(u, v, :) = image(u, v, :);
            end
        end
        
    end
end

disp('Completed infilling.'); toc;

%% Poisson image blending
% As done in both Cheng et al. and He et al. papers.
tic; disp('Starting poisson blending.');
try
    source = new_image;
    target = image;
    se = strel('ball', 50, 50);
    p_mask = im2bw(rescale(imdilate(double(~mask), se), 0, 1) > .5);

    bdry = 50;
    [nonzero_rows, nonzero_cols] = find(p_mask);
    min_row = min(nonzero_rows) - bdry; max_row = max(nonzero_rows) + bdry;
    min_col = min(nonzero_cols) - bdry; max_col = max(nonzero_cols) + bdry;

    source_region = source(min_row:max_row, min_col:max_col, :);
    mask_region = p_mask(min_row:max_row, min_col:max_col);
    target_region = target(min_row:max_row, min_col:max_col, :);

    result = poisson_blending(source_region, mask_region, target_region);
    new_image_blended = target;
    new_image_blended(min_row:max_row, min_col:max_col, :) = result(:,:,:);
    % new_image_2 = rgb2gray(new_image_2);
catch
    new_image_blended = new_image;
end
toc; disp('Finished poisson blending.');
%% Display Filled Image
plot_result(L_map_optimized,L_optimized_x,L_optimized_y, image, new_image_blended, mask);

%% Output images to output file
if super_pixel
    sp = '_sp';
else
    sp = '';
end

imwrite(h_image, sprintf('output/%s_masked.png', image_name))
imwrite(new_image, sprintf('output/%s_reconstructed%s.png', image_name, sp))
imwrite(new_image_blended, sprintf('output/%s_blended%s.png', image_name, sp))
