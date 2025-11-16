function [mask_texture] = loadMask()

disp('WELCOME to loadMask')

global w exp_num

% path to mask folder
mask_folder = fullfile(pwd,'stimuli',['mask_exp', num2str(exp_num)]);

% get texture of the mask
mask_texture = getTexturesFromHD(mask_folder, w);

end