function [mask_texture] = loadMask()

disp('WELCOME to loadMask')

global w

% path to mask folder
mask_folder = fullfile(pwd,'stimuli','mask');

% get texture of the mask
mask_texture = getTexturesFromHD(mask_folder, w);

end