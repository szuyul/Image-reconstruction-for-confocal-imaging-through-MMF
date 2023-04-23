function [ confocal_image_stack_H, confocal_image_stack_V ]... 
   = MMF_volumetric_reconstruction( sample, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions, T_fft_H, T_fft_V )
%% assemble TM at z = 0
n_z_positions = numel(z_positions);
% conversion between numerical defocus coeff. and physical distance
refocus_coeff_X = 7.1*z_positions; 
refocus_coeff_Y = 7.1*(z_positions+20);

[~, IA] = intersect(pre_mode_ind, H.mode_ind);  % mode_ind must be in pre_mode_ind, pre_mode_ind(IA) == mode_ind

delta = size(H.SF, 1) - 1;
step = delta/(2*img_physical_size);
shift_new_FOV = ceil((1-step)*img_digital_size/2);
F_size = numel( 1:step:(delta+1) );

ts = tic;
% interpolate in F domain for larger physical FOV
if step ~= 1 
    [H.dis_modeimg_F_Hill, H.dis_modeimg_F_Vill] = interpolate_F(H, step, F_pos_in_NA.H);
    [V.dis_modeimg_F_Hill, V.dis_modeimg_F_Vill] = interpolate_F(V, step, F_pos_in_NA.V);
end
toc(ts)

T_fw = [                           H.dis_modeimg_F_Hill(:, IA),                            sample.HV_ill_phase_offset*H.dis_modeimg_F_Vill(:, IA)  ;
        sample.HV_det_phase_offset*V.dis_modeimg_F_Hill(:, IA), sample.HV_ill_phase_offset*sample.HV_det_phase_offset*V.dis_modeimg_F_Vill(:, IA) ];
Tik_par = 0.1;
T_fw_inv  = Tikinv(T_fw, Tik_par);
% T_fw_inv  = T_fw'; 
%sample.T_2X = (sample.T_2X + sample.T_2X.')/2; % to make T_2X perfectly symmetric??

n = img_digital_size;
n_spatial_chs = numel(xy_pos_in_FOV);
n_dis_position_groups = ceil(n_spatial_chs/n);
n_input_sampling = size(sample.T_2X, 2); % for partial T2X measurement

figure
confocal_image_stack_H = zeros(img_digital_size, img_digital_size, n_z_positions);
confocal_image_stack_V = zeros(img_digital_size, img_digital_size, n_z_positions);
for ii = 1:n_z_positions
%% calculate the T_fw at refocused plane
    Fresnel_F_H = zern_aberr(F_size, [2;0], refocus_coeff_X(ii), F_pos_in_NA.H, [], [], 'vector');     % 08/30 H travels longer?
    Fresnel_F_V = zern_aberr(F_size, [2;0], refocus_coeff_Y(ii), F_pos_in_NA.V, [], [], 'vector');     

%% calculate the corresponding confocal images
    % refocus "after" matrix inversion since transfer function H is invertible
    T_fw_refocus_inv = T_fw_inv .* ([Fresnel_F_H.', Fresnel_F_V.']);
    confocal_image_H = complex(zeros(img_digital_size));
    confocal_image_V = complex(zeros(img_digital_size));

    for jj = 1:n_dis_position_groups
        if jj == n_dis_position_groups
            idx = (n*(jj-1)+1):n_spatial_chs;
        else
            idx = (n*(jj-1)+1):(n*jj); % address the last pixels chunk
        end
        % H ch
        in_gating = T_fw_refocus_inv(:, 1:size(T_fft_H, 1))*T_fft_H(:,idx);
        out_gating = in_gating.';
        temp = out_gating*(sample.T_2X)*in_gating(1:n_input_sampling,:);
        confocal_image_H(xy_pos_in_FOV(idx)) = diag(temp);

        % V ch
        in_gating = T_fw_refocus_inv(:, (size(T_fft_H, 1)+1):end)*T_fft_V(:,idx);
        out_gating = in_gating.';
        temp = out_gating*(sample.T_2X)*in_gating(1:n_input_sampling,:);
        confocal_image_V(xy_pos_in_FOV(idx)) = diag(temp);
    end
    confocal_image_stack_H(:,:, ii) = circshift(confocal_image_H, [shift_new_FOV, shift_new_FOV]);
    confocal_image_stack_V(:,:, ii) = circshift(confocal_image_V, [shift_new_FOV, shift_new_FOV]);

%% show the results
    gamma = 1/2.2;
    temp_H = abs(confocal_image_stack_H(:,:, ii)).^2;
    subplot(2, n_z_positions, ii);              complex_imagesc(imadjust(temp_H,[0 max(temp_H(:))],[0 max(temp_H(:))],gamma));
    temp_V = abs(confocal_image_stack_V(:,:, ii)).^2;
    subplot(2, n_z_positions, ii+n_z_positions);complex_imagesc(imadjust(temp_V,[0 max(temp_V(:))],[0 max(temp_V(:))],gamma));

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
    F_Hill_interp = complex(zeros(numel(F_pos_in_NA), n_pre_mode_ind));
    F_Vill_interp = complex(zeros(numel(F_pos_in_NA), n_pre_mode_ind));
    for ii = 1:n_pre_mode_ind
        temp(DC_aperture == 1) = P.dis_modeimg_F_Hill(:, ii);
        [temp_interp, ~] = complex_img_interpolate( temp, step );
        F_Hill_interp(:,ii) = temp_interp(F_pos_in_NA);
        
        temp(DC_aperture == 1) = P.dis_modeimg_F_Vill(:, ii);
        [temp_interp, ~] = complex_img_interpolate( temp, step );
        F_Vill_interp(:,ii) = temp_interp(F_pos_in_NA);
    end
end
