function [new_img, interp_img_size] = complex_img_interpolate( img, step )
% interpolate a complex image with arbitrary sampling rate
% input:
% img is a img_size by img_size complex image
% step is an arbitrary sampling period
%
% output:
% new_img is a interp_img_size by interp_img_size complex image after interpolation

    img_size = size(img,1);
    interp_img_size = numel(1:step:img_size);
    new_img = fftshift(fft2(ifft2(ifftshift(img)),interp_img_size, interp_img_size));
end
