%% load the data
clear 
load("two_layer_cell_data.mat")

%% reconstruct confocal images at varying z positions
z_positions = [0, 260]; % unit: um
n_ref = 1.44;
img_digital_size = 101; % odd number, maximal is delta + 1
img_physical_size = floor(delta/2); % unit: um

close all
if ~(exist('T_fft_H','var') && exist('T_fft_V','var'))
    ts = tic;
    [T_fft_H, T_fft_V, xy_pos_in_FOV, F_pos_in_NA] = T_fft_gen(img_digital_size, img_physical_size, H, V);
    toc(ts)
end
[ confocal_image_stacks_H, confocal_image_stacks_V ]... 
   = MMF_volumetric_reconstruction( MMF_signals_cell, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V );

%% combine H and V channels by image registration
close all
HV_image_register(confocal_image_stacks_H, confocal_image_stacks_V, z_positions);

%% reconstruct wide-field images
[ widefield_image_stack_H, widefield_image_stack_V ]... 
   = MMF_widefield_reconstruction( MMF_signals_cell, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V );

image_stacks_H = sqrt(widefield_image_stack_H);
image_stacks_V = sqrt(widefield_image_stack_V);

%% combine H and V channels by image registration
close all
HV_image_register(image_stacks_H, image_stacks_V, z_positions);

%% bright-field image
figure
subplot(1,2,1)
clean_image(bf_sample_s1_cell); title('surface 1')
subplot(1,2,2)
clean_image(bf_sample_s2_cell); title('surface 2')

%% show versatile image contrast (one z position at a time)
z_positions = 0;
n_ref = 1.44;

[ confocal_XH, confocal_YV, confocal_bi, confocal_di, confocal_H_scat, confocal_V_scat ]... 
   = MMF_confocal_multicontrasts( MMF_signals_cell, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V );

%% try the correction phase on X and Y polarization channels at the output
% the sample phase image should have reasonable phase variation without wrapping
phase_offset_X = [3.7, -15, -10, 0.5];
phase_offset_Y = [0, -22, -19, 0];
MMF_confocal_multicontrasts_plotting(confocal_XH, confocal_YV, confocal_bi, confocal_di, confocal_H_scat, confocal_V_scat,...
                                    phase_offset_X, phase_offset_Y);

%% apply the phase correction 
[ confocal_XH, confocal_YV, confocal_bi, confocal_di, confocal_H_scat, confocal_V_scat ]... 
   = MMF_confocal_multicontrasts( MMF_signals_cell, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V,...
                                  phase_offset_X, phase_offset_Y);

%% plot the corrected contrasts
phase_zeros = [0, 0, 0, 0];
MMF_confocal_multicontrasts_plotting(confocal_XH, confocal_YV, confocal_bi, confocal_di, confocal_H_scat, confocal_V_scat,...
                                    phase_zeros, phase_zeros);

  
%% self-defined functions
function [z_reflectivity_H, z_reflectivity_V] = HV_image_register(image_stacks_H, image_stacks_V, z_positions)
    figure('Position', [100, 100, 300, 600])
    figure('Position', [400, 100, 300, 600])
    h =  findobj('type','figure');
    optimizer = registration.optimizer.RegularStepGradientDescent;
    metric = registration.metric.MeanSquares;

    n_z_positions = numel(z_positions);
    
%%
    temp = abs(image_stacks_H).^2;
    z_reflectivity_H = sum(sum(temp, 1),2);
    H_max = max(max(max(temp)));
    temp = abs(image_stacks_V).^2;
    z_reflectivity_V = sum(sum(temp, 1),2);
    V_max = max(max(max(temp)));
%% 

    for ii = 1:n_z_positions
%% linear image
        cH = (abs(image_stacks_H(:,:,ii)).^2)/H_max;
        cV = (abs(image_stacks_V(:,:,ii)).^2)/V_max;
        
        cV_reg = imregister(cV,cH, 'translation', optimizer,metric);
        cHV = cH + cV_reg; 
        cHV = cHV/max(cHV(:));
        figure(length(h)-1)
        subplot(3, n_z_positions, ii);   clean_image(cH);      title('H'); caxis([0 1]); axis image
        subplot(3, n_z_positions, ii + n_z_positions);     clean_image(cV);    title('V');     caxis([0 1]); axis image
        subplot(3, n_z_positions, ii + 2*n_z_positions);   clean_image(cHV);   title('H+V');   caxis([0 1]); axis image
        
%% log-scaled image
        figure(length(h))
        subplot(3, n_z_positions, ii);                     clean_image(show_img_log10(cH)); title('H');
        subplot(3, n_z_positions, ii + n_z_positions);     clean_image(show_img_log10(cV));  title('V');
        subplot(3, n_z_positions, ii + 2*n_z_positions);   clean_image(show_img_log10(cHV));   title('H+V');
        
    end

%% integrated lateral reflectivity
    figure
    plot(z_positions, squeeze(z_reflectivity_H)); hold on
    plot(z_positions, squeeze(z_reflectivity_V));
    title('reflectivity plot')
    xlabel('z position (\mum)')
    ylabel('reflectivity (a.u.)')
end

function clean_image(img)
    imagesc(img); colormap('gray'); axis image; set(gca, 'XTick', [], 'YTick', []);
end

function log_img = show_img_log10(img)
    log_img = log10(img/max(img(:)));
    clean_image( log_img )
    caxis([-3 0])
end