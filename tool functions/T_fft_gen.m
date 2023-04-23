function [T_fft_H, T_fft_V, xy_pos_in_FOV, F_pos_in_NA] = T_fft_gen(img_digital_size, img_physical_size, H, V)
% generate inverse fft matrix that map frequency to space
%
% outputs:
% T_fft_H and T_fft_V are numel(modeimg_F) by # of output spatial chs
% xy_pos_in_FOV is index of selected spatial chs of the img_size by img_size output image
% F_pos_in_NA is a numel(modeimg_F) by 1 list specifying the frequency chs 
%
% inputs:
% img_digital_size  is the arbitrary number of sampling points along an output dimension
% img_physical_size is the arbitrary output physical dimension (unit: um)
% H and V are structures including all experimental variables

    [x, y] = meshgrid(1:img_digital_size);
    temp = sqrt( (x-round(img_digital_size/2)).^2 + (y-round(img_digital_size/2)).^2 ) < img_digital_size/0.5;
    xy_pos_in_FOV = find(temp == 1);
    
    [T_fft_H, F_pos_in_NA.H] = T_fft_assemble(H, xy_pos_in_FOV, img_digital_size, img_physical_size);
    [T_fft_V, F_pos_in_NA.V] = T_fft_assemble(V, xy_pos_in_FOV, img_digital_size, img_physical_size);
end

%% self-defined functions
function [T_fft, F_pos_in_NA] = T_fft_assemble(P, xy_pos_in_FOV, img_digital_size, img_physical_size)
% input:
% P is a structure of experimental variables
% xy_pos_in_FOV is index of selected spatial chs at the output
% img_digital_size  is the number of sampling points along an output dimension
% img_physical_size is the output physical dimension (unit: um)
% 
% output:
% T_fft is a fft matrix transforming between output spatial chs to freq. chs
% F_pos_in_NA is a list of counted freq. chs

    delta = size(P.dis_SF,1)-1;
    n_dis_pos = numel(xy_pos_in_FOV);
    DC_aperture = circshift(P.dis_SF, [P.dis_kxshift, P.dis_kyshift]);      % shift the F aperture to DC
    
    step = delta/(2*img_physical_size);
    if step == 1 % the original physical FOV
        F_pos_in_NA = find( DC_aperture == 1 );
        F_size = delta+1;
    else % for larger physical FOV
        [interp_DC_aperture, F_size] = DC_aperture_interp( DC_aperture, step );
        F_pos_in_NA = find( interp_DC_aperture >= 0.5 ); % for DC aperture
    end
    n_F = numel(F_pos_in_NA);

    T_fft = complex(zeros(n_F, n_dis_pos));
    for ii = 1:n_dis_pos
        img_out = zeros(img_digital_size);
        img_out( xy_pos_in_FOV(ii) ) = 1;
        temp = fftshift(fft2( imresize(img_out, [F_size, F_size]) ));       % interpolate in real-space pads 0s in F domain
        
        T_fft(:, ii) = temp(F_pos_in_NA);
    end
end

function [new_DC_aperture, interp_F_size] = DC_aperture_interp( DC_aperture, step )
% input:
% DC_aperture is a F_size by F_size image for masking DC term in F domain
% step is the re-sampling period in the F domain
%
% output:
% new_DC_aperture is the interpolated DC mask in F domain
% interp_F_size is the image size of the interpolated DC mask

    F_size = size(DC_aperture,1); 
    interp_F_size = numel(1:step:F_size);
    temp = fftshift(fft2(fftshift(DC_aperture)));
    n_pad_zeros = floor((interp_F_size - F_size)/2);
    temp_pad = padarray(temp,[n_pad_zeros, n_pad_zeros],0,'both');
    new_DC_aperture = abs(fftshift(fft2(fftshift(temp_pad))));
    % normalize the amplitude
    new_DC_aperture = new_DC_aperture/max(new_DC_aperture(:));
end
