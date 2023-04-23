%% load the data
clear 
load('quarter_WP_45deg_rotated_data.mat')

%% reconstruct confocal images
z_positions = 160; % unit: um
n_ref = 1;
img_digital_size = 101; % odd number, maximal is delta + 1
img_physical_size = floor(delta/2); % unit: um

close all
if ~(exist('T_fft_H','var') && exist('T_fft_V','var'))
    ts = tic;
    [T_fft_H, T_fft_V, xy_pos_in_FOV, F_pos_in_NA] = T_fft_gen(img_digital_size, img_physical_size, H, V);
    toc(ts)
end
[ confocal_image_stacks_H, confocal_image_stacks_V ]... 
   = MMF_volumetric_reconstruction( MMF_signals_qwp, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V );

%% combine H and V channels by image registration
close all
HV_image_register(confocal_image_stacks_H, confocal_image_stacks_V, z_positions);

%% reconstruct wide-field images
[ widefield_image_stack_H, widefield_image_stack_V ]... 
   = MMF_widefield_reconstruction( MMF_signals_qwp, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V );

image_stacks_H = sqrt(widefield_image_stack_H);
image_stacks_V = sqrt(widefield_image_stack_V);

%% combine H and V channels by image registration
close all
HV_image_register(image_stacks_H, image_stacks_V, z_positions);

%% bright-field image
figure
clean_image(bf_sample_wpq)
title('brightfield image of Q wave-plate')

%% show versatile image contrast
wave_plate_deg = 45;
z_positions = 160;
n_ref = 1;

close all
if ~(exist('T_fft_H','var') && exist('T_fft_V','var'))
    ts = tic;
    [T_fft_H, T_fft_V, xy_pos_in_FOV, F_pos_in_NA] = T_fft_gen(img_digital_size, img_physical_size, H, V);
    toc(ts)
end
[ confocal_XH, confocal_YV, confocal_ret, confocal_oax, confocal_H_scat, confocal_V_scat ]... 
   = MMF_confocal_multicontrasts( MMF_signals_qwp, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V );

%% try the correction phase on X and Y polarization channels at the output
% the sample phase image should have reasonable phase variation without wrapping
phase_offset_X = [1.8, -48, 1, 0];
phase_offset_Y = [0, -28, -1, 0];
MMF_confocal_multicontrasts_plotting(confocal_XH, confocal_YV, confocal_ret, confocal_oax, confocal_H_scat, confocal_V_scat,...
                                    phase_offset_X, phase_offset_Y);

%% apply the phase correction 
[ confocal_XH, confocal_YV, confocal_ret, confocal_oax, confocal_H_scat, confocal_V_scat ]... 
   = MMF_confocal_multicontrasts( MMF_signals_qwp, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions/n_ref, T_fft_H, T_fft_V,...
                                  phase_offset_X, phase_offset_Y);

%% plot the corrected contrasts
phase_zeros = [0, 0, 0, 0];
MMF_confocal_multicontrasts_plotting(confocal_XH, confocal_YV, confocal_ret, confocal_oax, confocal_H_scat, confocal_V_scat,...
                                    phase_zeros, phase_zeros);

mask = OA_mask_gen(size(confocal_oax,1), [48,54], 42); % for 45 deg
confocal_oax_unwrap = qu_phase(confocal_oax, mask, pi/2); % for 45 deg
   
figure
image(HueOverLum(confocal_oax_unwrap, confocal_ret, colormap(gca, cmap('C6')), [-pi, pi], [0 2])); 
axis image
title('oax(angle) + ret(amplitude)')

figure
draw_OA_color_map([-pi, pi], [0 1]) % OA_phase function should be only for colormap drawing
title('oax color map')


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

function draw_OA_color_map(phase, amp)
    [cmap_x, cmap_y] = meshgrid(linspace(-1, 1, 101),linspace(-1, 1, 101));
    circ_cmap = cmap_x + 1i*cmap_y;
    circ_cmap(sqrt(cmap_x.^2 + cmap_y.^2) > 1) = 0;
    circ_cmap_ang = angle(circ_cmap);
    circ_cmap_2ang = 2*circ_cmap_ang;
    circ_cmap_2ang( circ_cmap_ang>0 & circ_cmap_ang<pi )   = circ_cmap_2ang( circ_cmap_ang>0 & circ_cmap_ang<pi )-pi;
    circ_cmap_2ang( circ_cmap_ang>-pi & circ_cmap_ang<=0 ) = circ_cmap_2ang( circ_cmap_ang>-pi & circ_cmap_ang<=0 )+pi;
    image(HueOverLum(circ_cmap_2ang, abs(circ_cmap), colormap(cmap('C6')), phase, amp));
    set(gca,'XTick',[], 'YTick',[]);    axis image
end

function mask = OA_mask_gen(img_dim, center, r)
    mask = zeros(img_dim);
    [x,y] = meshgrid(1:img_dim);
    mask( sqrt( (x-center(1)).^2 + (y-center(2)).^2)< r ) = 1;
    %figure
    %imagesc(mask)
end

function circ_cmap_2ang = qu_phase(phase, mask, phase_center)
    temp = phase;
    condition_1 = abs(phase - (phase_center - pi))<0.8;
    condition_2 = abs(phase - (phase_center + pi))<0.8;
    temp( condition_1 ) = temp( condition_1 ) + pi;
    temp( condition_2 ) = temp( condition_2 ) - pi;
    temp( mask==0 ) = phase( mask==0 );
    circ_cmap_2ang = temp;
end

