function result = poisson_blending(source, mask, target)
source = padarray(source, [1,1], 'symmetric');
target = padarray(target, [1,1], 'symmetric');
mask = padarray(mask, [1,1]);

[t_rows, t_cols, ~] = size(target);

s = reshape(source, t_rows*t_cols, []);
t = reshape(target, t_rows*t_cols, []);

%% Allocate the RHS vector b.
b = zeros(t_rows*t_cols, 3);

row_vec = zeros(t_rows*t_cols, 1);
col_vec = zeros(t_rows*t_cols, 1);
value_vec = zeros(t_rows*t_cols, 1);

equation_num = 1;

for index = 1:t_rows*t_cols
    if mask(index)
        % 
        b(index,:) = 4*s(index,:) - s(index-1,:) - s(index+1,:) - s(index+t_rows,:) - s(index-t_rows,:);
        
        %% Insert a 4 into A at the index of the current central pixel.
        row_vec(equation_num) = index;
        col_vec(equation_num) = index;
        value_vec(equation_num) = 4;
        equation_num = equation_num + 1;
        
        %% Insert a -1 for the pixel below the current pixel:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index + 1;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;

        %% Insert a -1 for the pixel above the current pixel:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index - 1;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;

        %% Insert a -1 for the pixel to the left of the current pixel:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index - t_rows;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;
        
        %% Insert a -1 for the pixel to the right of the current pixel:    
        row_vec(equation_num) = index;
        col_vec(equation_num) = index + t_rows;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;
    else
        % If the current pixel location 
        row_vec(equation_num) = index;
        col_vec(equation_num) = index;
        value_vec(equation_num) = 1;
        equation_num = equation_num + 1;
        
        b(index,:) = t(index,:);
    end
end

%% We create the sparse matrix efficiently:
A = sparse(row_vec, col_vec, value_vec, t_rows*t_cols, t_rows*t_cols);

%% Solve for each color channel:
f_red = A \ b(:,1);
f_green = A \ b(:,2);
f_blue = A \ b(:,3);

%% Reshape to the original size:
f_red = reshape(f_red, [t_rows, t_cols]);
f_green = reshape(f_green, [t_rows, t_cols]);
f_blue = reshape(f_blue, [t_rows, t_cols]);

%% Stack the channels back together:
result = zeros(t_rows, t_cols, 3);
result(:,:,1) = f_red;
result(:,:,2) = f_green;
result(:,:,3) = f_blue;

%% Chop off the border:
result = result(2:t_rows-1, 2:t_cols-1, :);