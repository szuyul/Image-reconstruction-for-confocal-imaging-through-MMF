function MMF_confocal_multicontrasts_plotting(confocal_XH, confocal_YV, confocal_ret, confocal_oax, confocal_H_scat, confocal_V_scat, phase_offset_X, phase_offset_Y)
    
img_digital_size = size(confocal_XH, 1);

close all
figure('Position', [100, 100, 800, 1300])

% intensity imaging
gamma = 1/2.2;
temp = abs(confocal_XH).^2;
X_max = max(temp(:));
subplot(421); clean_image(log(1e-3+temp/X_max)); caxis([-5 0]); title('intensity in X'); colorbar
temp = abs(confocal_YV).^2;
Y_max = max(temp(:));
subplot(422); clean_image(log(1e-3+temp/Y_max)); caxis([-5 0]); title('intensity in Y'); colorbar

% phase imaging
[~, aberr] = zern_aberr(img_digital_size, [0,1,1,2;0,1,-1,0], phase_offset_X, []);
temp = confocal_XH.*aberr;
%temp = temp./abs(temp);
subplot(423); complex_imagesc( temp ); title('phase in X')
[~, aberr] = zern_aberr(img_digital_size, [0,1,1,2;0,1,-1,0], phase_offset_Y, []);
temp = confocal_YV.*aberr;
%temp = temp./abs(temp);
subplot(424); complex_imagesc( temp ); title('phase in Y')

% retardance + optic axis
subplot(425); clean_image( confocal_ret ); caxis([0 2*pi]); title('retardance'); colorbar
subplot(426); clean_image( confocal_oax ); caxis([-pi pi]); title('optic axis'); colorbar

% scattering
gamma = 1/2.2;
temp = confocal_H_scat;
subplot(427); clean_image( log(1e-3 + temp/max(temp(:))) ); caxis([-5 0]); title('scattering in X'); colorbar
temp = confocal_V_scat;
subplot(428); clean_image( log(1e-3 + temp/max(temp(:))) ); caxis([-5 0]); title('scattering in Y'); colorbar

end

function clean_image(img, clrmap, varargin)
    imagesc(img); axis image; set(gca, 'XTick', [], 'YTick', []);
    if nargin == 1
        colormap('gray');
    else
        colormap(clrmap);
    end
end