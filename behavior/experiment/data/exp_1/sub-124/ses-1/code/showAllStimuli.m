function [] = showAllStimuli(file_list)

global gray_color w ScreenHeight originalHeight originalWidth ScreenWidth category spaceBar kb

Screen('FillRect', w, gray_color);
drawFrame();

% check which images to show
if strcmp(category, ' ')
    show_cates = ["kitchen", "bathroom"];
else
    show_cates = category;
end

% define image size
image_height = round(ScreenHeight/6 - ScreenHeight*0.01);
image_width = round((image_height/originalHeight) * originalWidth);

for show_cate = {show_cates}
    show_cate = char(show_cate);

    % get texture pointers
    texture_array = [file_list.raw.(show_cate).texture];

    % in how many slides should it be devided (each 5x% grid of images)
    slides = ceil(length(texture_array)/25);

    count = 0;
    for slide = 1:slides

        % make a 5x5 grid of images
        for row = 1:5
            for col = 1:5
                count = count + 1;

                % get center of that image
                center_of_image = [(ScreenWidth/6)*col; (ScreenHeight/6)*row];

                % make image dimensions
                x = transpose(center_of_image) - [image_width/2 image_height/2];
                y = transpose(center_of_image) + [image_width/2 image_height/2];

                % draw image to texture
                Screen('DrawTexture',w, texture_array(count),[],[x y]);

            end
        end

        % flip texture
        Screen('Flip', w, [], 1);

        WaitSecs(1)

        % wait for space bar
        [~, ~, wait_resp] = KbCheck(kb);
        while ~wait_resp(spaceBar)
            [~, ~, wait_resp] = KbCheck(kb);
        end

    end
end
end