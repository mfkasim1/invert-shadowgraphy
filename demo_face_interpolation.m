% This file is to demonstrate the face interpolation using forward and inverse
% algorithms
% Author: Muhammad Firmansyah Kasim (University of Oxford, 2018)

add_libs;

% load the face images
face1 = mean(double(imread('test-figures/face-01-gray.png')), 3);
face2 = mean(double(imread('test-figures/face-02-gray.png')), 3);

% get the face part of the image
mask1 = (face1 < 255);
mask2 = (face2 < 255);

% get the mean values of the faces
mean1 = mean(face1(mask1));
mean2 = mean(face2(mask2));

% set the background to the mean value
face1(~mask1) = mean1;
face2(~mask2) = mean2;

clear options;
options.alpha = 0.1;
options.minstep = 1e-3;
options.interp = 'linear';
phi = main_inverse(face1, face2, options);

for i = linspace(0, 1, 40)
  imagesc(main_forward(face1, phi*i));
  colormap gray;
  pbaspect([1,1,1]);
  xticks([]);
  yticks([]);
  pause(0.05);
end
