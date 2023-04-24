%% load the data
clear
load('resolution_chart.mat')
MMF_signals_chart = MMF_signals_chart_600;

%% reconstruct confocal images
switch MMF_signals_chart.z
    case 0 % for chart placed at 0 um
        img_physical_size = floor(delta/2); % unit: um 

    case 600 % for chart placed at 600 um
        img_physical_size = floor(1.5*delta/2); % unit: um 

    case 1200 % for chart placed at 1200 um
        img_physical_size = floor(2*delta/2); % unit: um 

    otherwise
end
n_ref = 1.4378;
img_digital_size = 101; % original delta + 1

close all
if ~(exist('T_fft_H','var') && exist('T_fft_V','var'))
    [T_fft_H, T_fft_V, xy_pos_in_FOV, F_pos_in_NA] = T_fft_gen(img_digital_size, img_physical_size, H, V);
end
[ confocal_image_stacks_H, confocal_image_stacks_V ]... 
   = MMF_volumetric_reconstruction( MMF_signals_chart, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, MMF_signals_chart.z/n_ref, T_fft_H, T_fft_V );

%% show the intensity image
norm_image = normalized_image_I(confocal_image_stacks_H, 3);
figure
imagesc(log(norm_image))
colormap('gray'); axis image;

%% self-defined functions
function norm_image = normalized_image_I(field, flip_dim, varargin)
% normalize the image intensity and flip the image if needed
    temp = abs(field).^2;
    max_int = max(temp(:));
    if nargin == 1
        norm_image = temp/max_int;
    else
        switch flip_dim
            case 3
                norm_image = flip(temp/max_int, 1);
                norm_image = flip(norm_image, 2);
            otherwise
                norm_image = flip(temp/max_int, flip_dim);
        end
    end
end
