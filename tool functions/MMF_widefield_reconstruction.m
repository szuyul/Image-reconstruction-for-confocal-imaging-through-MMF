function [ widefield_image_stack_H, widefield_image_stack_V ]... 
   = MMF_widefield_reconstruction( sample, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions, T_fft_H, T_fft_V )
%% determine dense confocal gated spatial channels
delta = size(H.SF, 1) - 1;

%% assemble TM at z = 0
n_z_positions = numel(z_positions);
% conversion between numerical defocus coeff. and physical distance
refocus_coeff_X = 7.1*z_positions; 
refocus_coeff_Y = 7.1*(z_positions+10);

[~, IA] = intersect(pre_mode_ind, H.mode_ind);  % mode_ind must be in pre_mode_ind, pre_mode_ind(IA) == mode_ind

step = delta/(2*img_physical_size);
shift_new_FOV = ceil((1-step)*img_digital_size/2);
F_size = delta+1;

% interpolate in F domain for larger physical FOV
if step ~= 1 
    [H.dis_modeimg_F_Hill, H.dis_modeimg_F_Vill] = interpolate_F(H, step, F_pos_in_NA.H);
    [V.dis_modeimg_F_Hill, V.dis_modeimg_F_Vill] = interpolate_F(V, step, F_pos_in_NA.V);
    [~, F_size] = complex_img_interpolate( zeros(delta+1), step );
end

T_fw = [                           H.dis_modeimg_F_Hill(:, IA),                            sample.HV_ill_phase_offset*H.dis_modeimg_F_Vill(:, IA)  ;
        sample.HV_det_phase_offset*V.dis_modeimg_F_Hill(:, IA), sample.HV_ill_phase_offset*sample.HV_det_phase_offset*V.dis_modeimg_F_Vill(:, IA) ];
Tik_par = 0.15;
T_fw_inv  = Tikinv(T_fw, Tik_par); 

figure
widefield_image_stack_H = zeros(img_digital_size, img_digital_size, n_z_positions);
widefield_image_stack_V = zeros(img_digital_size, img_digital_size, n_z_positions);
for ii = 1:n_z_positions
%% calculate the T_fw at refocused plane
Fresnel_F_H = zern_aberr(F_size, [2;0], refocus_coeff_X(ii), F_pos_in_NA.H, [], [], 'vector');     
Fresnel_F_V = zern_aberr(F_size, [2;0], refocus_coeff_Y(ii), F_pos_in_NA.V, [], [], 'vector');     
            
%% calculate the corresponding confocal images
T_fw_refocus_inv = T_fw_inv .* ([Fresnel_F_H.', Fresnel_F_V.']);
widefield_image_H = zeros(img_digital_size);
widefield_image_V = zeros(img_digital_size);

%n = 200; % process n spatial channels at a time
n = img_digital_size;
n_spatial_chs = numel(xy_pos_in_FOV);
n_dis_position_groups = ceil(n_spatial_chs/n);
for jj = 1:n_dis_position_groups
    if jj == n_dis_position_groups
        idx = (n*(jj-1)+1):n_spatial_chs;
    else
        idx = (n*(jj-1)+1):(n*jj);
    end
    % H ch
    in_gating = T_fw_refocus_inv(:, 1:size(T_fft_H, 1))*T_fft_H(:,idx);
    out_gating = in_gating.';
    temp = abs(out_gating*sample.T_2X).^2;
    widefield_image_H(xy_pos_in_FOV(idx)) = sum(temp,2);
    % V ch
    in_gating = T_fw_refocus_inv(:, (size(T_fft_H, 1)+1):end)*T_fft_V(:,idx);
    out_gating = in_gating.';
    temp = abs(out_gating*sample.T_2X).^2;
    widefield_image_V(xy_pos_in_FOV(idx)) = sum(temp,2);
end
widefield_image_stack_H(:,:, ii) = circshift(widefield_image_H, [shift_new_FOV, shift_new_FOV]);
widefield_image_stack_V(:,:, ii) = circshift(widefield_image_V, [shift_new_FOV, shift_new_FOV]);

%% show the results
subplot(2, n_z_positions, ii);               imagesc( widefield_image_stack_H(:,:, ii) ); axis image; colormap('gray')
subplot(2, n_z_positions, ii+n_z_positions); imagesc( widefield_image_stack_V(:,:, ii) ); axis image; colormap('gray') 
end

end

%% self-defined function
function [F_Hill_interp, F_Vill_interp] = interpolate_F(P, step, F_pos_in_NA)
% interpolate the F domain, and output numel(F_pos_in_NA) by n_pre_mode_ind TM
%
% input:
% P is a structure of experimental variables
% step is the re-sampling period in the F domain
% F_pos_in_NA is the list specifying the frequency chs of the interpolated mask in F domain
% 
% output:
% F_Hill_interp and F_Vill_interp are interpolated output frequency chs per input realization

    DC_aperture = circshift(P.dis_SF, [P.dis_kxshift, P.dis_kyshift]);      % shift the F aperture to DC

    temp = zeros(size(DC_aperture));
    n_pre_mode_ind = size(P.dis_modeimg_F_Hill,2);
    F_Hill_interp = zeros(numel(F_pos_in_NA), n_pre_mode_ind);
    F_Vill_interp = zeros(numel(F_pos_in_NA), n_pre_mode_ind);
    for ii = 1:n_pre_mode_ind
        temp(DC_aperture == 1) = P.dis_modeimg_F_Hill(:, ii);
        [temp_interp, ~] = complex_img_interpolate( temp, step );
        F_Hill_interp(:,ii) = temp_interp(F_pos_in_NA);
        
        temp(DC_aperture == 1) = P.dis_modeimg_F_Vill(:, ii);
        [temp_interp, ~] = complex_img_interpolate( temp, step );
        F_Vill_interp(:,ii) = temp_interp(F_pos_in_NA);
    end
end
