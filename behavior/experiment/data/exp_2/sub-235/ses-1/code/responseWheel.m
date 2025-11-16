% adapted from Daniel Kaiser (https://doi.org/10.1162/jocn_a_01891)

function [resp, resp_time] = responseWheel(tr)

global w kb center stimSizeHeight white black green gray_color spaceBar ScreenHeight ScreenWidth abortKey task

%% set up

% all participants should have same random numbers
rng(tr)

stimpositionmid=center;
fixposition=stimpositionmid;

%show cursor
SetMouse(fixposition(1),fixposition(2));
ShowCursor('Hand');

% rectangle surounding the middle
picrectangle1=[stimpositionmid(1)-stimSizeHeight,stimpositionmid(2)-stimSizeHeight];
picrectangle2=[stimpositionmid(1)+stimSizeHeight,stimpositionmid(2)+stimSizeHeight];
middle_rec=[picrectangle1,picrectangle2];  %Rectangle for the Picture

% rectangle in the middle
sf=0.98; % scaling factor
picrectangle1=[stimpositionmid(1)-stimSizeHeight*sf,stimpositionmid(2)-stimSizeHeight*sf];
picrectangle2=[stimpositionmid(1)+stimSizeHeight*sf,stimpositionmid(2)+stimSizeHeight*sf];
middle_rec2=[picrectangle1,picrectangle2];  %Rectangle for the heel

% rectangle for rating wheel
sf=1.35;
picrectangle1=[stimpositionmid(1)-stimSizeHeight*sf,stimpositionmid(2)-stimSizeHeight*sf];
picrectangle2=[stimpositionmid(1)+stimSizeHeight*sf,stimpositionmid(2)+stimSizeHeight*sf];
wheelrectangle=[picrectangle1,picrectangle2];  %Rectangle for the heel

% large ractangle sourding rating wheel
sf=1.37;
picrectangle1=[stimpositionmid(1)-stimSizeHeight*sf,stimpositionmid(2)-stimSizeHeight*sf];
picrectangle2=[stimpositionmid(1)+stimSizeHeight*sf,stimpositionmid(2)+stimSizeHeight*sf];
large_rec=[picrectangle1,picrectangle2];  %Rectangle for the heel

a1=360/8;
start_angle=360*rand;  % A random start
cm=colormap(gray(9))*255;

%%
% large ractangle sourding rating wheel
Screen('FillOval', w, black, large_rec);

% response options
for i=1:7
    a0=start_angle+(i-1)*a1;
    Screen('FillArc',w,cm(i+2,:),wheelrectangle,a0,a1);
    r=1.2*(wheelrectangle(3)-wheelrectangle(1))/2;
    x_text=fixposition(1)+r*cos(deg2rad(start_angle+(i-1)*a1-1.5*a1));
    y_text=fixposition(2)+r*sin(deg2rad(start_angle+(i-1)*a1-1.5*a1));
    DrawFormattedText(w,num2str(i),x_text,y_text,black);
end

% make inner black line of response wheel
Screen('FillOval', w, black, middle_rec);

% leave one position open
a0=start_angle+(8-1)*a1;
Screen('FillArc',w,gray_color,large_rec,a0,a1);
r=1.2*(large_rec(3)-large_rec(1))/2;

% fill inner gray part of response wheel
Screen('FillOval', w, gray_color, middle_rec2);

% make text and flip
DrawFormattedText(w,task,'center','center',black);
Screen('Flip',w);

% keep waiting for them to select something
resp = 0;
while resp == 0 || keyCode(spaceBar) ~= 1

    [~,resp_time,keyCode] = KbCheck(kb);
    [x,y,buttons]=GetMouse;

    % add option to abort using esc
    if keyCode(abortKey)
        error('Experiment has been aborted');
    end

    % when mouse is clicked
    if buttons(1)~=0

        % check all boxes
        for i=1:7
            % make boxes
            r=max([ScreenHeight, ScreenWidth]);
            ax0=deg2rad(start_angle+(i-1)*a1-90);
            ax1=deg2rad(start_angle+i*a1-90);
            if ax0>ax1
                t=linspace(ax1,ax0);
            else
                t=linspace(ax0,ax1);
            end

            x_arc=[fixposition(1),fixposition(1)+r*cos(t),fixposition(1)];
            y_arc=[fixposition(2),fixposition(2)+r*sin(t),fixposition(2)];
            cursor_dist=dist([x,y;fixposition(1),fixposition(2)]');

            % if response is in this box store box
            if inpolygon(x,y,x_arc,y_arc)==1 && cursor_dist(2)>=(middle_rec(3)-middle_rec(1))/2 && cursor_dist(2)<=(wheelrectangle(3)-wheelrectangle(1))/2

                Screen('FillOval', w, black, large_rec);

                for j=1:7
                    a0=start_angle+(j-1)*a1;
                    % check if box is the selected box
                    if j==i
                        % make green if yes
                        Screen('FillArc',w,green,wheelrectangle,a0,a1);
                    else
                        Screen('FillArc',w,cm(j+2,:),wheelrectangle,a0,a1);
                    end
                    r=1.2*(wheelrectangle(3)-wheelrectangle(1))/2;
                    x_text=fixposition(1)+r*cos(deg2rad(start_angle+(j-1)*a1-1.5*a1));
                    y_text=fixposition(2)+r*sin(deg2rad(start_angle+(j-1)*a1-1.5*a1));
                    DrawFormattedText(w,num2str(j),x_text,y_text,black);
                end

                % make inner black line of response wheel
                Screen('FillOval', w, black, middle_rec);

                % leave one position open
                a0=start_angle+(8-1)*a1;
                Screen('FillArc',w,gray_color,large_rec,a0,a1);
                r=1.2*(large_rec(3)-large_rec(1))/2;

                % fill inner gray part of response wheel
                Screen('FillOval', w, gray_color, middle_rec2);

                % make text and flip
                DrawFormattedText(w,task,'center','center',black);
                Screen('Flip',w);

                % store response
                resp=i;
            end
        end
    end
end
HideCursor;

KbWait(-1);


end