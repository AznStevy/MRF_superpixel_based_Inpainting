function [image]=gray2rgb(image)
    %Gives a grayscale image an extra dimension
    %in order to use color within it
    [m n] = size(image);
    rgb = zeros(m,n,3);
    rgb(:,:,1)=image;
    rgb(:,:,2)=rgb(:,:,1);
    rgb(:,:,3)=rgb(:,:,1);
    image=rgb/255;
end
