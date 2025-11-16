% In this script, the experiment runtime parameters are being set. This
% is where the recording modes, hardware and physical parameters are being
% set. Before running the experiment, this file should be checked to be
% sure that the parameters are set correctly.
% This script has no input nor ouput. Every parameters being set are
% global variables. This means that they can be called by the different
% functions. 

function initRuntimeParameters
% Recording modalities:
global Behavior EYE_TRACKER 
% Photodiode parameters:
global PHOTODIODE DIOD_DURATION DIOD_SIZE DIOD_ON_COLOUR DIOD_OFF_COLOUR FRAME_ANTICIPATION
% Eyetracker parameters:
global DISTANCE_SCREEN_TRACKER HEAD_FIXED
% Hardware parameters:
global SCREEN_SIZE_CM REF_RATE_OPTIMAL viewDistanceBottomTop viewDistance 
% Debugging and code parameters:
global SHOW_PRACTICE DEBUG RESOLUTION_FORCE NO_FULLSCREEN WINDOW_RESOLUTION SHOW_INSTRUCTIONS 
% experimental parameters
global task sub_num
% keys
global CalibrationKey ValidationKey NO_KEY RESTART_KEY ABORT_KEY abortKey upKey downKey PauseKey RestartKey YesKey
global oneKey twoKey threeKey fourKey fiveKey sixKey sevenKey eightKey spaceBar 
global F_Key J_Key KITCHEN_KEY BATHROOM_KEY RightKey LeftKey valid_resp_keys
% text
global fontType fontSize fontColor
% colors
global black gray red green
% optics
global FRAME_WIDTH MAX_VISUAL_ANGEL VIEWING_DISTANCE FRAME_COLOR FIXATION_COLOR FIXATION_FONT_SIZE DIAMOUT_FIXATION DIAMIN_FIXATION FILE_POSTFIX
% other
global TRUE FALSE

%% DEBUG parameters
DEBUG = 0; % 0 = no debug | 1 = regular debug | 2 = fast debug
SHOW_INSTRUCTIONS = 1;
SHOW_PRACTICE = 1; % skip the practice run
RESOLUTION_FORCE = 0; % the program will complain if optimal refresh rate is not possible on this screen
NO_FULLSCREEN = 0; % enable windowed mode for dubugging
WINDOW_RESOLUTION = [10 10 1200 800];
FRAME_ANTICIPATION = 0.1; % optimizes timing in PTB

%% Recording modalities
EYE_TRACKER = 0; % Must be set to 1 if recording with Eyetracker
Behavior = 1; %Set to 1 if recording with Behavior only

%% Hardware and physical parameters:
REF_RATE_OPTIMAL = 60; % in Hz. Screen refresh rate.
viewDistance = 70; % Default viewing distance (if no viewDist argument sent with the function call)
SCREEN_SIZE_CM = [53.2 29.8]; % screen [width, height] in centimeters, change it to fit your setting
viewDistanceBottomTop = [72 72]; % IN CM!! Distance between the participant head and the top and bottom of the screen. Only needed if HEAD_FIXED on.

%% Eyetracker parameters:
DISTANCE_SCREEN_TRACKER = 90; % Distance between the eyetracker lense and the computer screen. Only needed for REMOTE MODE!
HEAD_FIXED = 0; % Head fixed must be set to 0 if remote mode

%% Photodiode parameters:
PHOTODIODE = 0; % Must be set to 1 for the photodiode to be presented
DIOD_ON_COLOUR = 255; % Color of the photodiode when turned on (255 white, 0 black)
DIOD_OFF_COLOUR = 1;  % Color of the photodiode when off (255 white, 0 black)
DIOD_SIZE = 100; % Size of the square where the photodiode is presented (in pixels)
DIOD_DURATION = 3; % Duration of the photodiode flash when turned on (in frames)

%% Constant parameters

% Add a frame for the bebug mode to not mix things up:
FRAME_WIDTH = 0; % 0 for delete
if DEBUG; FRAME_WIDTH = 1; end
FRAME_COLOR = [39,241,44];

% TEXT
fontType = 'David';
fontSize = 50; % general text size, not target
FIXATION_COLOR = [205 33 42];
FIXATION_FONT_SIZE = 20;
fontColor = 0; % black;

% color
black = [255, 255, 255];
gray = [125, 125, 125];
red = [255, 0, 0];
green = [0, 255, 0];

% image format
FILE_POSTFIX = '*tif';

VIEWING_DISTANCE = viewDistance; % in centimeters
MAX_VISUAL_ANGEL = [6,6]; % in degrees | "on a rectangular aperture at an average visual angle of 6? by 4?"

% Size of the fixation in DVA:
DIAMOUT_FIXATION = 0.6; % diameter of outer circle (degrees)
DIAMIN_FIXATION = 0.1; % diameter of inner circle (degrees)

%% Response params
KbName('UnifyKeyNames');
CalibrationKey= KbName('C');
upKey         =  KbName('UpArrow');
downKey       =  KbName('DownArrow');
RightKey      =  KbName('RightArrow');
LeftKey       =  KbName('LeftArrow');
PauseKey      =  KbName('Q');
RestartKey    =  KbName('R');
abortKey      =  KbName('ESCAPE'); % ESC aborts experiment
YesKey        =  KbName('Y');
F_Key         =  KbName('F');
J_Key         =  KbName('J');
spaceBar      =  KbName('SPACE');
oneKey        =  KbName('1!');
twoKey        =  KbName('2@');
threeKey      =  KbName('3#');
fourKey       =  KbName('4$');
fiveKey       =  KbName('5%');
sixKey        =  KbName('6^');
sevenKey      =  KbName('7&');
eightKey      =  KbName('8*');
ValidationKey = KbName('V');

if mod(sub_num, 2)
    % odd subject number
    KITCHEN_KEY = F_Key;
    BATHROOM_KEY = J_Key;
else
    % even subject number
    KITCHEN_KEY = J_Key;
    BATHROOM_KEY = F_Key;
end


% define which keys are valid response keys
if strcmp(task, 'categorization')

    % cagtegorization response keys
    valid_resp_keys = [KITCHEN_KEY, BATHROOM_KEY];

else

    % rating response key
    valid_resp_keys = [oneKey, twoKey,threeKey, fourKey, fiveKey, sixKey, sevenKey];

end

% program codes
ABORT_KEY = 14;
RESTART_KEY = 13;
NO_KEY = 0;

TRUE = 1;
FALSE = 0;

end