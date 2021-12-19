function [Tx_RGB, Ty_RGB]= ImageRegistration(thresholdValue)
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
%
% FOR STUDENTS
%
% This function registers the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' and returns the shifts for each image and each layer
% Red, Green and Blue. The shifts are calculated relatively to the first
% image 'LR_Tiger_01.tif'. Each low-resolution image is 64 x64 pixels.
%
%
% OUTPUT:   Tx_RGB: horizontal shifts, a 40x3 matrix
%           Ty_RGB: vertical shifts, a 40x3 matrix
%


% NOTE: _Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%       _Tx_RGB(20,2) is the horizontal shift of the Green layer of the
%       20th image relatively to the Green layer of the firs image.
%
%
% OUTLINE OF THE ALGORITHM:
%
% 1.The first step is to compute the continuous moments m_00, m_01 and m_10
% of each low-resolution image using the .mat file called:
% PolynomialReproduction_coef.mat. This file contains three matrices
% 'Coef_0_0', 'Coef_1_0' and 'Coef_0_1' used to calculate the continuous
% moments.
%
% 2.The second step consists in calculating the barycenters of the Red,
% Green and Blue layers of the low-resolution images.
%
% 3.By computing the difference between the barycenters of corresponding
% layers between two images, the horizontal and vertical shifts can be
% retrieved for each layer.
%
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************


% Load the coefficients for polynomial reproduction
load('PolynomialReproduction_coef.mat','Coef_0_0','Coef_1_0','Coef_0_1');

% -------- include your code here -----------

% thresholdValue = 0.1;

Number_of_Imgs = 40;
baryCenterX = zeros(Number_of_Imgs,3);
baryCenterY = zeros(Number_of_Imgs,3);


% Import the low-resolution RGB image
for picture_index = 1:1:Number_of_Imgs
    picture_name = ['LR_Tiger_' sprintf('%02d', picture_index)  '.tif'];
    LR_img = imread(picture_name);
    % rescale the image from 0-255 to 0-1 to reduce the noise created by the
    % background of the image and improve the accuracy of moments
    LR_img = rescale(LR_img, 0, 1);
    % threshold the normalised image
    LR_img = LR_img .* (LR_img >= thresholdValue);
    
    % correspond to the R,G and B layers
    
    for layer = 1:1:3
        LR_img_layer = LR_img(:,:,layer);
        
        % calculate the moments
        m00 = sum(sum(Coef_0_0 .* LR_img_layer));
        m01 = sum(sum(Coef_0_1 .* LR_img_layer));
        m10 = sum(sum(Coef_1_0 .* LR_img_layer));
        
        % calculate the bary center position
        baryCenterX(picture_index, layer) = m10/m00;
        baryCenterY(picture_index, layer) = m01/m00;
    end
end

Tx_RGB = zeros(Number_of_Imgs,3);
Ty_RGB = zeros(Number_of_Imgs,3);
for layer = 1:1:3
    Tx_RGB(:, layer) = baryCenterX(:, layer) - baryCenterX(1, layer);
    Ty_RGB(:, layer) = baryCenterY(:, layer) - baryCenterY(1, layer);
end

end


