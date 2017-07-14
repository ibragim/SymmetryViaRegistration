Image = imread('gantrycrane.png');
Image = imrotate(Image,45,'crop');

[angles, midPoints, segLengths, IOut] = symmetryViaRegistration2D(Image,'RegMethod','nxc');

subplot(1,2,1)
imshow(Image), title('Input')

subplot(1,2,2)
imshow(IOut), title('Output')