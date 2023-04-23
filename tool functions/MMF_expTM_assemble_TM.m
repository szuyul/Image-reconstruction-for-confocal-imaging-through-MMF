function [T_HH, T_VH, T_HV, T_VV, corrco, phased] = MMF_expTM_assemble_TM(img_stacks, H, V, n_positions, pre_mode_ind)
% input:
% img_stacks is a (delta+1) by (delta+1) by 3*(# of scanned channels) by 2(H and V detection) 4D matrix of acquired data
% n_positions is the total # of focal lattice
% 
% output:
% T_HH{1} is a (# of F channels) by (# of scanned channels) matrix
% T_HH{2} is a delta+1 by delta+1 by n_positions matrix of fully sampled TM
% corrco is a (# of scanned channels) by 1 complex numbers of ref. mode correlation trace
% phased is a (# of scanned channels) by 1 doubles of ref. mode phase trace

H_img_stack = img_stacks(:,:,:, 1); % H detection
V_img_stack = img_stacks(:,:,:, 2); % V detection

%% H detection, H illumination
T_HH = cell(1,2);
[T_HH{1}, T_HH{2}, corrco, phased] = Hilbert_transform(H_img_stack, H, 1, n_positions, pre_mode_ind);
%% V detection, H illumination
T_VH = cell(1,2);
[T_VH{1}, T_VH{2}] = Hilbert_transform(V_img_stack, V, 1, n_positions, pre_mode_ind);
%% H detection, V illumination
T_HV = cell(1,2);
[T_HV{1}, T_HV{2}] = Hilbert_transform(H_img_stack, H, 2, n_positions, pre_mode_ind);
%% V detection, V illumination 
T_VV = cell(1,2);
[T_VV{1}, T_VV{2}] = Hilbert_transform(V_img_stack, V, 2, n_positions, pre_mode_ind);

%% monitor the phase drift
n_scanned_positions = 1:numel(pre_mode_ind);

close all
figure('Position', [100, 100, 500, 400])
plot(n_scanned_positions, abs(corrco)); hold on
yyaxis right 
plot(n_scanned_positions, phased)
legend({'correlation';'phase drift'})
xlim([0, numel(pre_mode_ind)])
grid on

fprintf('\nProgram paused. Press enter to continue.\n');
pause;
end

function [modeimg_F, modeimg, corrco, phased] = Hilbert_transform(img_stack, P, t_offset, n_positions, pre_mode_ind)
% input: 
% img_stack is a delta+1 by delta+1 by 3*numel(pre_mode_ind) 3D matrix of raw camera fringed images of MMF in response
%    to (H, V, H ref...) inputs in a single polarization state detection
% t_offset: 
%       1 = H illumination at mode position 
%       2 = V illumination at mode position
%       3 = H illumination at reference position
% n_positions are total # of input focus lattice
% 
% output:
% modeimg_F is a (# of F channels) by numel(pre_mode_ind) matrix of the chosen polarization detection in F domain
% modeimg is a delta+1 by delta+1 by n_positions matrix of full MMF complex output
% corrco is a numel(pre_mode_ind) by 1 complex correlation trace for monitoring phase drifting 

n_frames = size(img_stack,3);
num_spatial_channels = n_frames/3;

mode_stack = img_stack(:,:, t_offset:3:(3*(num_spatial_channels-1)+t_offset));
ref_stack = img_stack(:,:, 3:3:end);

delta = size(mode_stack,1)-1;

modeimg_F =   zeros( sum(sum( P.SF )), num_spatial_channels);                       % for full TM saving [# of output modes in k space, # of input modes]
modeimg = zeros(delta+1,delta+1, n_positions);

corrco = zeros(num_spatial_channels, 1);

Fccd = fftshift(fft2( ref_stack(:,:,1) ));
comref_F = Fccd( P.SF == 1 );
ref_F_m = norm( comref_F );
corrco(1) = 1;

for ii = 1:num_spatial_channels 
    % reference
    Fccd = fftshift(fft2( ref_stack(:,:,ii) ));
    moderef_F = Fccd( P.SF == 1 );
    imag_F_m = norm( moderef_F );
    corrco(ii) = sum(sum( comref_F.*conj(moderef_F) )) / ( ref_F_m*imag_F_m );
    % mode
    Fccd = fftshift(fft2( mode_stack(:,:,ii) ));
    modeimg_F(:, ii) = Fccd( P.SF == 1 );
    
    modeimg(:,:,pre_mode_ind(ii)) = ifft2(ifftshift(circshift((P.SF.*Fccd),[P.kxshift P.kyshift]))); 
end
phased = unwrap(angle(corrco));
%phasedint = phased;
%phasedint = spline(3:3:3*num_spatial_channels, phased, t_offset:3:(3*(num_spatial_channels-1)+t_offset)); % this will lead to error in the beginning extrapolation 

modeimg_F   = bsxfun(@times,modeimg_F,  reshape(exp(1i*phased(1:num_spatial_channels)),[1,num_spatial_channels]));
modeimg(:,:,pre_mode_ind) = bsxfun(@times,modeimg(:,:,pre_mode_ind),reshape(exp(1i*phased(1:num_spatial_channels)),[1,1,num_spatial_channels]));

end