function [ confocal_XH, confocal_YV, confocal_ret, confocal_oax, confocal_H_scat, confocal_V_scat ]... 
   = MMF_confocal_multicontrasts( sample, img_digital_size, img_physical_size, xy_pos_in_FOV, F_pos_in_NA, H, V, pre_mode_ind, z_positions, T_fft_X, T_fft_Y,...
                                  phase_offset_X, phase_offset_Y)
%% physical z to numerical coefficients
n_z_positions = numel(z_positions);
% conversion between numerical defocus coeff. and physical distance
refocus_coeff_X = 7.1*z_positions; 
refocus_coeff_Y = 7.1*(z_positions);

[~, IA] = intersect(pre_mode_ind, H.mode_ind);  % mode_ind must be in pre_mode_ind, pre_mode_ind(IA) == mode_ind

delta = size(H.SF, 1) - 1;
step = delta/(2*img_physical_size);
shift_new_FOV = ceil((1-step)*img_digital_size/2);
F_size = delta+1;

%% interpolate in F domain for larger physical FOV
if step ~= 1 
    [H.dis_modeimg_F_Hill, H.dis_modeimg_F_Vill] = interpolate_F(H, step, F_pos_in_NA.H);
    [V.dis_modeimg_F_Hill, V.dis_modeimg_F_Vill] = interpolate_F(V, step, F_pos_in_NA.V);
    [~, F_size] = complex_img_interpolate( zeros(delta+1), step );
end

T_fw = [                           H.dis_modeimg_F_Hill(:, IA),                            sample.HV_ill_phase_offset*H.dis_modeimg_F_Vill(:, IA)  ;
        sample.HV_det_phase_offset*V.dis_modeimg_F_Hill(:, IA), sample.HV_ill_phase_offset*sample.HV_det_phase_offset*V.dis_modeimg_F_Vill(:, IA) ];
Tik_par = 0.1;
T_fw_inv = Tikinv(T_fw, Tik_par); 

sample.T_2X = (sample.T_2X + sample.T_2X.')/2;

%% calculate images of versatile contrasts from reflection matrix at an OP
pos_offset = zeros(img_digital_size);
for ff = 1:img_digital_size
    pos_offset(:,ff) = abs( ((1-ff):(img_digital_size-ff)) );
end
confocal_XH = zeros(img_digital_size); % complex image
confocal_YV = zeros(size(confocal_XH));

confocal_ret = zeros(img_digital_size); % real-valued iamge
confocal_oax = zeros(img_digital_size); 

confocal_H_scat = zeros(size(confocal_XH)); % real-valued image
confocal_V_scat = zeros(size(confocal_XH));

%%%%%%%% correction to distal misalignment %%%%%%%%
if nargin == 11
    Fresnel_X = ones(img_digital_size^2,1);
    Fresnel_Y = ones(img_digital_size^2,1);
else
    Fresnel_X = zern_aberr(img_digital_size, [0,1,1,2;0,1,-1,0], phase_offset_X/2, 1:img_digital_size^2, [], [], 'vector').';
    Fresnel_Y = zern_aberr(img_digital_size, [0,1,1,2;0,1,-1,0], phase_offset_Y/2, 1:img_digital_size^2, [], [], 'vector').'; 
end
Fresnel_X_2 = Fresnel_X.^2;
Fresnel_Y_2 = Fresnel_Y.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for ii = 1:n_z_positions
%% calculate the T_fw at refocused plane
Fresnel_F_H = zern_aberr(F_size, [2;0], refocus_coeff_X(ii), F_pos_in_NA.H, [], [], 'vector');     % 08/30 H travels longer?
Fresnel_F_V = zern_aberr(F_size, [2;0], refocus_coeff_Y(ii), F_pos_in_NA.V, [], [], 'vector');     
            
%% calculate the corresponding confocal images
% refocus "after" matrix inversion since transfer function H is invertible
T_fw_refocus_inv = T_fw_inv .* ([Fresnel_F_H.', Fresnel_F_V.']);
confocal_image_X = complex(zeros(img_digital_size));
confocal_image_Y = complex(zeros(img_digital_size));

n = img_digital_size;
n_spatial_chs = numel(xy_pos_in_FOV);
n_dis_position_groups = ceil(n_spatial_chs/n);
n_input_sampling = size(sample.T_2X, 2); % for partial T2X measurement

for jj = 1:n_dis_position_groups
    if jj == n_dis_position_groups
        idx = (n*(jj-1)+1):n_spatial_chs;
    else
        idx = (n*(jj-1)+1):(n*jj); % address the last pixels chunk
    end
%% co-polarization (phase)
    % X ch
    in_gating_X = T_fw_refocus_inv(:, 1:size(T_fft_X, 1))*T_fft_X(:,idx);
    out_gating_X = in_gating_X.';
    temp_XX = out_gating_X*(sample.T_2X)*in_gating_X(1:n_input_sampling,:) *diag(Fresnel_X_2(idx)) ;
    confocal_image_X(xy_pos_in_FOV(idx)) = diag(temp_XX);
    % multiple-scattering
    confocal_H_scat(xy_pos_in_FOV(idx)) = sum(pos_offset.*(abs(temp_XX).^2),1)./sum(abs(temp_XX).^2,1);
    
    % Y ch
    in_gating_Y = T_fw_refocus_inv(:, (size(T_fft_X, 1)+1):end)*T_fft_Y(:,idx);
    out_gating_Y = in_gating_Y.';
    temp_YY = out_gating_Y*(sample.T_2X)*in_gating_Y(1:n_input_sampling,:) *diag(Fresnel_Y_2(idx)) ;
    confocal_image_Y(xy_pos_in_FOV(idx)) = diag(temp_YY);
    % multiple-scattering
    confocal_V_scat(xy_pos_in_FOV(idx)) = sum(pos_offset.*(abs(temp_YY).^2),1)./sum(abs(temp_YY).^2,1);
    
%% cross-polarization (birefringence)
    temp_YX = diag(Fresnel_Y(idx)) * (out_gating_Y*(sample.T_2X)*in_gating_X(1:n_input_sampling,:)) * diag(Fresnel_X(idx)) ;
    temp_XY = diag(Fresnel_X(idx)) * (out_gating_X*(sample.T_2X)*in_gating_Y(1:n_input_sampling,:)) * diag(Fresnel_Y(idx)) ;
    
    % original naive method
    % see Eqs. 14 & 16 in de Boer, et al., "Polarization sensitive optical coherence tomography - a review"
    %[confocal_bi(xy_pos_in_FOV(idx)), confocal_di(xy_pos_in_FOV(idx))] = retardance_R(diag(temp_HH), diag(temp_VH), diag(temp_HV), diag(temp_VV));
    
    Js = reshape([diag(temp_XX), diag(temp_YX), diag(temp_XY), diag(temp_YY)].', [4, n]);
    %{
    % 
    Js = reshape([diag(temp_XX), diag(temp_YX), diag(temp_XY), diag(temp_YY)].', [4, n]);
    [r,d] = JonesDecomp(Js, true);
    confocal_bi_1(xy_pos_in_FOV(idx)) = r(1,:);
    confocal_di_1(xy_pos_in_FOV(idx)) = d(1,:);
    confocal_bi_2(xy_pos_in_FOV(idx)) = r(2,:);
    confocal_di_2(xy_pos_in_FOV(idx)) = d(2,:);
    confocal_bi_3(xy_pos_in_FOV(idx)) = r(3,:);
    confocal_di_3(xy_pos_in_FOV(idx)) = d(3,:);
    %}
    [ret, oax] = oax_from_J(Js);
    confocal_ret(xy_pos_in_FOV(idx)) = ret(:);
    confocal_oax(xy_pos_in_FOV(idx)) = oax(:);
    
    jj
end
confocal_XH(:,:, ii) = circshift(confocal_image_X, [shift_new_FOV, shift_new_FOV]);
confocal_YV(:,:, ii) = circshift(confocal_image_Y, [shift_new_FOV, shift_new_FOV]);

%% show the results
gamma = 1/2.2;
temp_H = abs(confocal_XH(:,:, ii)).^2;
subplot(2, n_z_positions, ii);              complex_imagesc(imadjust(temp_H,[0 max(temp_H(:))],[0 max(temp_H(:))],gamma));
temp_V = abs(confocal_YV(:,:, ii)).^2;
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

function [retard, diatten] = retardance_R(R_HH, R_VH, R_HV, R_VV)
    % see Eqs. 14 & 16 in de Boer, et al., "Polarization sensitive optical coherence tomography - a review"
    % R_XX are n by 1 col. vector (diag. from reflection submatrix)
    n = numel(R_HH);
    Js = reshape([R_HH, R_VH, R_HV, R_VV].', [2, 2, n]); % Jones matrix of each spatial channel
    retard = zeros(n, 1);
    diatten = zeros(n,1);
    for ii = 1:n
        [~, vals] = eig(Js(:,:,ii), 'vector'); % 2 eigenvalues in a column vector
        retard(ii) = abs(angle(vals(1)/vals(2)))-pi/2;
        diatten(ii) = abs( abs(vals(1))^2 - abs(vals(2))^2 ) / (abs(vals(1))^2 + abs(vals(2))^2);
    end
end

function [ret, oax] = oax_from_J(Js)
% Js = Jr*Jd, polar decomp for each 2-by-2 Jones matrix 
%   Js = U S V'
%   Jr = U V'
%   Jd = V S V'
% 
% input:
% Js is a 4-by-n_matrices matrix, each column is a 2x2 J matrix
% 
% output:
% ret is a 3-by-n_matrices matrix, 
% oax is a 3-by-n_matrices matrix, 

    [~, Jd2] = JonesDecomp(MatrixMultiply(conj(Js([1;3;2;4],:)),Js)); % Jd2 = Js'*Js = V*S^2*V'
    Jd = makeJones(Jd2/2); 
    Jdinv = Jd([4;2;3;1],:).*[1;-1;-1;1];
    rr = JonesDecomp(makeJones(JonesDecomp(MatrixMultiply(Js,Jdinv))));
    %rr is a 'cleaned' retardation vector
    %rr_amp = sqrt(sum(rr.^2,1));
    %m = (rr_amp+pi - mod(rr_amp+pi,2*pi))/(2*pi);
    %rr = rr  - 2*pi*m.*rr./rr_amp;
    rr2 = rr  - 2*pi*rr./sqrt(sum(rr.^2,1));
    rr(repmat(sqrt(sum(rr.^2,1))>pi, [3,1])) = rr2(repmat(sqrt(sum(rr.^2,1))>pi, [3,1])); 
    % unwrap higher-order rr into base band, while keep rr already in base band untouched

    ret = sqrt(sum(rr.^2,1));
    oax = atan2(rr(1,:),rr(2,:));
end

function  Out = MatrixMultiply(A,B)
% computes matrix products, assuming matrices linearized along first dimension, in column order
% input: 
% A and B are N^2-by-n_matrices matrices 
% 
% output: 
% Out is a N^2-by-n_matrix, the matrix product of A and B in each column

N = sqrt(size(A,1));% number of elements along one side of the matrix

dimA = size(A);
dimB = size(B);

dimA = [dimA,ones(1,numel(dimB)-numel(dimA))];
dimB = [dimB,ones(1,numel(dimA)-numel(dimB))];

S1.type = '()';
S1.subs = repmat({':'},1,ndims(A));
S2.type = '()';
S2.subs = repmat({':'},1,ndims(B));
 
Out = zeros(max(dimA,dimB));
for ind3 = 1:N
    ind1 = repmat((ind3-1)*N+1:ind3*N,[1,N])';
    ind2 = repmat((ind3:N:N*N),[N,1]);
    S1.subs{1} = ind1;
    S2.subs{1} = ind2;
    Out = Out + subsref(A,S1).*subsref(B,S2);
end
end

function [r,d] = JonesDecomp(J,polar,pureDiatt)
% JonesDecomp computes the retardation and diattenuation from a Jones matrix, 
% provided as input argument (the matrix logarithm is part of the function, using the 'concurrent' decomposition.)
%
% input: 
% J is either 2x2, or 4 x whatever, where the 4 elements are [J11;J21;J12;J22]
% polar is optional argument, it performs the polar decomposition, decomposing the Jones matrix 
%   into a sequence of a diattenuation and a retardation matrix, J = Jr*Jd.
% 
% output: 
% r is a 3-by-whatever retardation matrix
% d is a 3-by-whatever diattenuation matrix

dim = size(J);

if nargin<2
    polar = false;
end
if nargin<3
    pureDiatt = false;
end

if ~polar
    if dim(1) == 2 && dim(2) == 2
        J = J/sqrt(det(J));
        q = cat(1,(J(1,1)-J(2,2)),(J(2,1)+J(1,2)),(-1i*J(2,1)+1i*J(1,2)))/2;
        tr = trace(J)/2;
        c = acosh(tr);
        csin = c/sinh(c);
        csin(c==0) = 1;
        f = 2*q*csin;
        r = -imag(f);
        d = real(f);

    elseif dim(1)==4
        detJ = sqrt(J(1,:).*J(4,:)-J(2,:).*J(3,:));
        J = bsxfun(@rdivide,J(:,:),detJ);
        q = cat(1,(J(1,:)-J(4,:)),(J(2,:)+J(3,:)),(-1i*J(2,:)+1i*J(3,:)))/2;
        tr = (J(1,:) + J(4,:))/2;
        c = acosh(tr);
        csin = c./sinh(c);
        csin(c==0) = 1;
        f = 2*bsxfun(@times,q,csin);
%        f = 2*bsxfun(@times,q,c./sinh(c));
        r = reshape(-imag(f),[3,dim(2:end)]);
        d = reshape(real(f),[3,dim(2:end)]);
    end
else% polar decomposition
    if dim(1) == 2 && dim(2) == 2
        J = J/sqrt(det(J));
        Jd = J'*J;
        [~,d] = JonesDecomp(Jd);
        d = d/2;
        r = JonesDecomp(J*makeJones(zeros(size(d)),-d));
    elseif dim(1)==4
        detJ = sqrt(J(1,:).*J(4,:)-J(2,:).*J(3,:));
        J = bsxfun(@rdivide,J(:,:),detJ);
        [~,d] = JonesDecomp(MatrixMultiply(conj(J([1;3;2;4],:)),J));
        d = d/2;
        r = JonesDecomp(MatrixMultiply(J,makeJones(zeros(size(d)),-d)));
        % makeJones(zeros(size(d)),-d) generate J matrix with 0 retardation and inverse diattenuation
        r = reshape(r,[3,dim(2:end)]);
        d = reshape(d,[3,dim(2:end)]);
    end
end
end

function Jout = makeJones(r,d)
%Jout = makeJones(r,d) generates the general Jones matrix, defined by the retaration 'r' and diattenuation 'd'
%
% input:
% r is a 3-by-whatever matrix (an array of any matching dimension in the additional dimensions)
% d is a 3-by-whatever matrix
% If r is only a three component vector, Jout will be reshape to 2x2 elements, 
%    otherwise Jout is 4 along the first dimension, and matches r for the additional dimensions.
%
% output:
% Jout is a 4-by-whatever matrix, with each column as a 2x2 matrix 

if nargin<2 || isempty(d)
    d = zeros(size(r));
end

dim = size(r);

f = (d - 1i*r)/2;
c = sqrt(sum(f.^2));

sinch = sinh(c)./c;
sinch(c==0) = 1;
Jout = bsxfun(@times,[1;0;0;1],cosh(c(:,:))) + sinch(:,:).*(bsxfun(@times,[1;0;0;-1],f(1,:)) + bsxfun(@times,[0;1;1;0],f(2,:)) + bsxfun(@times,[0;1i;-1i;0],f(3,:)));

if numel(r)==3
    Jout = reshape(Jout,[2,2]);
else
    Jout = reshape(Jout,[4,dim(2:end)]);
end
end
