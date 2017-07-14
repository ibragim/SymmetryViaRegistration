Image = imread('gantrycrane.png');
Image = imrotate(Image,45,'crop');

[angles, midPoints, segLengths, IOut] = symmetryViaRegistration2D(Image,'RegMethod','dense');

subplot(1,2,1)
imshow(Image), title('Input')

subplot(1,2,2)
imshow(IOut), title('Output')

% Reference:
% Finding Mirror Symmetry via Registration
% Marcelo Cicconet, David G. C. Hildebrand, Hunter Elliott
% https://arxiv.org/abs/1611.05971