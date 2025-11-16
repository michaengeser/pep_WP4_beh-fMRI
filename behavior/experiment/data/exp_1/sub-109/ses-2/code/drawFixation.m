%DRAWFIXATION
%Draws fixation on the screen buffer without "flipping"

function [ ] = drawFixation(innerCross)

global  ppd w center DIAMOUT_FIXATION DIAMIN_FIXATION white black

% if no color is specified its black
if nargin<1
    innerCross = white;
end

%Change size based on dva
cx = center(1,1);
cy = center(2,1);

Screen('FillOval', w, black, [cx-DIAMOUT_FIXATION/2*ppd, cy-DIAMOUT_FIXATION/2*ppd, cx+DIAMOUT_FIXATION/2*ppd, cy+DIAMOUT_FIXATION/2 * ppd], DIAMOUT_FIXATION*ppd);
Screen('FillOval', w, innerCross, [cx-DIAMIN_FIXATION/2*ppd, cy-DIAMIN_FIXATION/2*ppd, cx+DIAMIN_FIXATION/2*ppd, cy+DIAMIN_FIXATION/2 * ppd], DIAMOUT_FIXATION*ppd);

% comment out cross
%Screen('DrawLine', w, innerCross, cx-DIAMOUT_FIXATION/2*ppd, cy,cx+DIAMOUT_FIXATION/2*ppd, cy, DIAMIN_FIXATION*ppd);
%Screen('DrawLine', w, innerCross, cx, cy-DIAMOUT_FIXATION/2*ppd,cx, cy+DIAMOUT_FIXATION/2 * ppd, DIAMIN_FIXATION*ppd);
% 
% global w FIXATION_FONT_SIZE FIXATION screenScaler fontSize FIXATION_COLOR
% 
% Screen('TextSize', w, round(FIXATION_FONT_SIZE*screenScaler));
% 
% DrawFormattedText(w, textProcess(FIXATION), 'center', 'center', FIXATION_COLOR);
% 
% Screen('TextSize', w, round(fontSize*screenScaler));

end