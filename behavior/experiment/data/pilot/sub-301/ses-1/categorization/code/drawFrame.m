%DRAWFRAME
% draws a frame around the stimuli
function [] = drawFrame()

    global FRAME_COLOR
    global w gray_color FRAME_WIDTH x_pos_frame y_pos_frame

    Screen('FillRect', w, FRAME_COLOR, [x_pos_frame y_pos_frame] + [-FRAME_WIDTH -FRAME_WIDTH FRAME_WIDTH FRAME_WIDTH]);
    Screen('FillRect', w, gray_color, [x_pos_frame y_pos_frame]);

end
