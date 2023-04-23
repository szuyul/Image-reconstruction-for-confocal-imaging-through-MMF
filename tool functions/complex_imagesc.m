%complex_imagesc plots the complex image on the current figure
%
% input: 
% img is an 2D image either flattened or not
% amp_range is the range of amplitude for colormap
%
% output: in the current figure, plot the complex field without ticks
%
%
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [  ] = complex_imagesc(img, amp_range, varargin)

if ismember(1, size(img))
    img_dim = sqrt(length( img ));
    oimg = reshape( img , [img_dim img_dim]);
else 
    oimg = img;
end

imgmax = max(max(abs(oimg)));
field = oimg/imgmax;
if nargin == 2    % if the amplitude range is given
    image(HueOverLum(angle(field), abs(field), colormap(gca, cmap('C6')), [-pi, pi], amp_range)); 
else              % otherwise normalize the amplitude
    image(HueOverLum(angle(field), abs(field), colormap(gca, cmap('C6')), [-pi, pi], [0 1])); 
end

set(gca,'XTick',[], 'YTick',[]);    axis image

%{
%% generate complex phase colormap
[cmap_x, cmap_y] = meshgrid(linspace(-1, 1, 101),linspace(-1, 1, 101));
circ_cmap = cmap_x + 1i*cmap_y;
circ_cmap(sqrt(cmap_x.^2 + cmap_y.^2) > 1) = 0;
axes('pos',[.14 .14 .1 .1])
image(HueOverLum(angle(circ_cmap), abs(circ_cmap), colormap(cmap('C6')), [-pi, pi], [0 1]));
set(gca,'XTick',[], 'YTick',[]);    axis image
%}
end

