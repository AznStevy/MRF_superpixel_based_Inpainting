function [A] = build_grid(M, N, weight)
    if nargin < 3; weight = 1; end;
    A = zeros(N*M); % we need a xlen*ylen by xlen*ylen matrix

    for i = 1:M
        for j = 1:N               
            k = sub2ind([M N],i,j);
            if i > 1
                ii=i-1; jj=j; % up
                A(k,sub2ind([M N],ii,jj)) = weight;
            end
            if i < M
                ii=i+1; jj=j; % down
                A(k,sub2ind([M N],ii,jj)) = weight;
            end
            if j > 1
                ii=i; jj=j-1; % left
                A(k,sub2ind([M N],ii,jj)) = weight; 
            end
            if j < N
                ii=i; jj=j+1; % right
                A(k,sub2ind([M N],ii,jj)) = weight;
            end
        end
    end
end

