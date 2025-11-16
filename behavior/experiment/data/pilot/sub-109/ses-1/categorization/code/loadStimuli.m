
% LOADSTIMULI loads all the stimuli in all types and orientation from HD. Wasteful and should be used before the experiment is running.
% output:
% -------
% loads all textures from HD, and puts them into the global textures
% structure (texture_struct) matching the folder structure of the stimuli
% folder

function [] = loadStimuli(file_list)

disp('WELCOME to loadStimuli')

global texture_struct w

% path to stimuli folder
PreFolderName = [pwd,filesep,'stimuli\'];
img_types = fieldnames(file_list);
cate = fieldnames(file_list.(img_types{1}));
texture_struct = file_list;

% loops through the folders an loads all stimuli
for j = 1:length(img_types)
    for jj = 1:length(cate)

        FolderName = fullfile(PreFolderName, img_types{j}, cate{jj});
        new_textures = getTexturesFromHD(FolderName, w);
        texture_struct.(img_types{j}).(cate{jj}) = new_textures;

        % check if texture struct matches file_list
        for check = 1:length(texture_struct.(img_types{j}).(cate{jj}))
            textue_struct_texture = texture_struct.(img_types{j}).(cate{jj})(check);
            file_list_texture = file_list.(img_types{j}).(cate{jj})(check).texture;

            % if textures are not identical, give error message
            if ~(textue_struct_texture == file_list_texture)
                showMessage('Textures are NOT identical!');
                WaitSecs(2);
            end

        end
    end
end
