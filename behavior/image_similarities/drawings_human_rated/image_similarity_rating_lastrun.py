#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2023.2.3),
    on Sat May 18 13:32:11 2024
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

import psychopy
psychopy.useVersion('2023.2.3')


# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout
from psychopy.tools import environmenttools
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard

# Run 'Before Experiment' code from init_block_message
block_message = "Block 0"
# --- Setup global variables (available in all functions) ---
# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# Store info about the experiment session
psychopyVersion = '2023.2.3'
expName = 'image_similarity_rating'  # from the Builder filename that created this script
expInfo = {
    'participant': '200',
    'session': ['bathroom','kitchen'],
    'session_number': '',
    'date': data.getDateStr(),  # add a simple timestamp
    'expName': expName,
    'psychopyVersion': psychopyVersion,
}


def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # temporarily remove keys which the dialog doesn't need to show
    poppedKeys = {
        'date': expInfo.pop('date', data.getDateStr()),
        'expName': expInfo.pop('expName', expName),
        'psychopyVersion': expInfo.pop('psychopyVersion', psychopyVersion),
    }
    # show participant info dialog
    dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    # restore hidden keys
    expInfo.update(poppedKeys)
    # return expInfo
    return expInfo


def setupData(expInfo, dataDir=None):
    """
    Make an ExperimentHandler to handle trials and saving.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    dataDir : Path, str or None
        Folder to save the data to, leave as None to create a folder in the current directory.    
    Returns
    ==========
    psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'data/sub-%s/ses-%s/sub-%s_ses-%s_%s_events' % (expInfo['participant'], expInfo['session'], expInfo['participant'], expInfo['session'], expName)
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='/Users/kaiserlab/Documents/GitHub/pep_wp4/image_similarities/drawings_human_rated/image_similarity_rating_lastrun.py',
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='time'
    )
    thisExp.setPriority('thisRow.t', priority.CRITICAL)
    thisExp.setPriority('expName', priority.LOW)
    # return experiment handler
    return thisExp


def setupLogging(filename):
    """
    Setup a log file and tell it what level to log at.
    
    Parameters
    ==========
    filename : str or pathlib.Path
        Filename to save log file and data files as, doesn't need an extension.
    
    Returns
    ==========
    psychopy.logging.LogFile
        Text stream to receive inputs from the logging system.
    """
    # this outputs to the screen, not a file
    logging.console.setLevel(logging.EXP)
    # save a log file for detail verbose info
    logFile = logging.LogFile(filename+'.log', level=logging.EXP)
    
    return logFile


def setupWindow(expInfo=None, win=None):
    """
    Setup the Window
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    win : psychopy.visual.Window
        Window to setup - leave as None to create a new window.
    
    Returns
    ==========
    psychopy.visual.Window
        Window in which to run this experiment.
    """
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=[1280, 720], fullscr=True, screen=0,
            winType='pyglet', allowStencil=False,
            monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='height'
        )
        if expInfo is not None:
            # store frame rate of monitor if we can measure it
            expInfo['frameRate'] = win.getActualFrameRate()
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [0,0,0]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'height'
    win.mouseVisible = False
    win.hideMessage()
    return win


def setupInputs(expInfo, thisExp, win):
    """
    Setup whatever inputs are available (mouse, keyboard, eyetracker, etc.)
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window in which to run this experiment.
    Returns
    ==========
    dict
        Dictionary of input devices by name.
    """
    # --- Setup input devices ---
    inputs = {}
    ioConfig = {}
    ioSession = ioServer = eyetracker = None
    
    # create a default keyboard (e.g. to check for escape)
    defaultKeyboard = keyboard.Keyboard(backend='event')
    # return inputs dict
    return {
        'ioServer': ioServer,
        'defaultKeyboard': defaultKeyboard,
        'eyetracker': eyetracker,
    }

def pauseExperiment(thisExp, inputs=None, win=None, timers=[], playbackComponents=[]):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    inputs : dict
        Dictionary of input devices by name.
    win : psychopy.visual.Window
        Window for this experiment.
    timers : list, tuple
        List of timers to reset once pausing is finished.
    playbackComponents : list, tuple
        List of any components with a `pause` method which need to be paused.
    """
    # if we are not paused, do nothing
    if thisExp.status != PAUSED:
        return
    
    # pause any playback components
    for comp in playbackComponents:
        comp.pause()
    # prevent components from auto-drawing
    win.stashAutoDraw()
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # make sure we have a keyboard
        if inputs is None:
            inputs = {
                'defaultKeyboard': keyboard.Keyboard(backend='Pyglet')
            }
        # check for quit (typically the Esc key)
        if inputs['defaultKeyboard'].getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win, inputs=inputs)
        # flip the screen
        win.flip()
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, inputs=inputs, win=win)
    # resume any playback components
    for comp in playbackComponents:
        comp.play()
    # restore auto-drawn components
    win.retrieveAutoDraw()
    # reset any timers
    for timer in timers:
        timer.reset()


def run(expInfo, thisExp, win, inputs, globalClock=None, thisSession=None):
    """
    Run the experiment flow.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    psychopy.visual.Window
        Window in which to run this experiment.
    inputs : dict
        Dictionary of input devices by name.
    globalClock : psychopy.core.clock.Clock or None
        Clock to get global time from - supply None to make a new one.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    # mark experiment as started
    thisExp.status = STARTED
    # make sure variables created by exec are available globally
    exec = environmenttools.setExecEnvironment(globals())
    # get device handles from dict of input devices
    ioServer = inputs['ioServer']
    defaultKeyboard = inputs['defaultKeyboard']
    eyetracker = inputs['eyetracker']
    # make sure we're running in the directory for this experiment
    os.chdir(_thisDir)
    # get filename from ExperimentHandler for convenience
    filename = thisExp.dataFileName
    frameTolerance = 0.001  # how close to onset before 'same' frame
    endExpNow = False  # flag for 'escape' or other condition => quit the exp
    # get frame duration from frame rate in expInfo
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0  # could not measure, so guess
    
    # Start Code - component code to be run after the window creation
    
    # --- Initialize components for Routine "preparation" ---
    # Run 'Begin Experiment' code from scale_set_up
    oldt=0
    x_size=8.560
    y_size=5.398
    screen_height=0
    
    if win.units=='norm':
        x_scale=.05
        y_scale=.1
        dbase = .0001
        unittext=' norm units'
        vsize=2
    elif win.units=='pix':
        x_scale=60
        y_scale=40
        dbase = .1
        unittext=' pixels'
        vsize=win.size[1]
    else:
        x_scale=.05
        y_scale=.05
        dbase = .0001
        unittext=' height units'
        vsize=1
    
    # --- Initialize components for Routine "instructions" ---
    startResponse = keyboard.Keyboard()
    slide_img = visual.ImageStim(
        win=win,
        name='slide_img', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=(0, 0), size=(1.7, 1.7*0.56),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    
    # --- Initialize components for Routine "trial" ---
    # Run 'Begin Experiment' code from positions
    img_x_pos = 7
    img1pos = -img_x_pos
    img2pos = img_x_pos
    
    polygon = visual.Rect(
        win=win, name='polygon',
        width=(8.2*x_scale, 8.2*x_scale*0.75)[0], height=(8.2*x_scale, 8.2*x_scale*0.75)[1],
        ori=0.0, pos=(img1pos*x_scale, 4*y_scale), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor='black', fillColor='black',
        opacity=None, depth=-1.0, interpolate=True)
    polygon_2 = visual.Rect(
        win=win, name='polygon_2',
        width=(8.2*x_scale, 8.2*x_scale*0.75)[0], height=(8.2*x_scale, 8.2*x_scale*0.75)[1],
        ori=0.0, pos=(img2pos*x_scale, 4*y_scale), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor='black', fillColor='black',
        opacity=None, depth=-2.0, interpolate=True)
    image1 = visual.ImageStim(
        win=win,
        name='image1', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-3.0)
    image2 = visual.ImageStim(
        win=win,
        name='image2', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=512.0, interpolate=True, depth=-4.0)
    rating_1 = visual.Slider(win=win, name='rating_1',
        startValue=None, size=1.0, pos=[0,0], units=win.units,
        labels=("1 - not similar at all", 2, 3, 4, 5, 6, "7 - very similar"), ticks=(1, 2, 3, 4, 5, 6, 7), granularity=1.0,
        style='rating', styleTweaks=('triangleMarker',), opacity=None,
        labelColor='white', markerColor='black', lineColor='White', colorSpace='rgb',
        font='Verdana', labelHeight=0.35*y_scale,
        flip=False, ori=0.0, depth=-5, readOnly=False)
    question = visual.TextStim(win=win, name='question',
        text='Content',
        font='Verdana',
        pos=[0,0], height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-6.0);
    rating_2 = visual.Slider(win=win, name='rating_2',
        startValue=None, size=1.0, pos=[0,0], units=win.units,
        labels=("1 - not similar at all", 2, 3, 4, 5, 6, "7 - very similar"), ticks=(1, 2, 3, 4, 5, 6, 7), granularity=1.0,
        style='rating', styleTweaks=('triangleMarker',), opacity=None,
        labelColor='white', markerColor='black', lineColor='White', colorSpace='rgb',
        font='Verdana', labelHeight=0.35*y_scale,
        flip=False, ori=0.0, depth=-7, readOnly=False)
    question_2 = visual.TextStim(win=win, name='question_2',
        text='Drawing style/ability',
        font='Verdana',
        pos=[0,0], height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-8.0);
    block_text = visual.TextStim(win=win, name='block_text',
        text='',
        font='Open Sans',
        pos=(16*x_scale, 9.8*y_scale), height=0.02, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-10.0);
    
    # --- Initialize components for Routine "ITI" ---
    text = visual.TextStim(win=win, name='text',
        text=None,
        font='Open Sans',
        pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    
    # --- Initialize components for Routine "practice_over" ---
    text_2 = visual.TextStim(win=win, name='text_2',
        text='Practice is over.\nWell done!\n\nPress space to proceed to experiment.',
        font='Open Sans',
        pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    space_practice_over = keyboard.Keyboard()
    
    # --- Initialize components for Routine "breakBlocks" ---
    # Run 'Begin Experiment' code from break_2
    blockN = 1
    maxBlock = 10
    
    startBlock = keyboard.Keyboard()
    
    # --- Initialize components for Routine "ITI" ---
    text = visual.TextStim(win=win, name='text',
        text=None,
        font='Open Sans',
        pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    
    # --- Initialize components for Routine "trial" ---
    # Run 'Begin Experiment' code from positions
    img_x_pos = 7
    img1pos = -img_x_pos
    img2pos = img_x_pos
    
    polygon = visual.Rect(
        win=win, name='polygon',
        width=(8.2*x_scale, 8.2*x_scale*0.75)[0], height=(8.2*x_scale, 8.2*x_scale*0.75)[1],
        ori=0.0, pos=(img1pos*x_scale, 4*y_scale), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor='black', fillColor='black',
        opacity=None, depth=-1.0, interpolate=True)
    polygon_2 = visual.Rect(
        win=win, name='polygon_2',
        width=(8.2*x_scale, 8.2*x_scale*0.75)[0], height=(8.2*x_scale, 8.2*x_scale*0.75)[1],
        ori=0.0, pos=(img2pos*x_scale, 4*y_scale), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor='black', fillColor='black',
        opacity=None, depth=-2.0, interpolate=True)
    image1 = visual.ImageStim(
        win=win,
        name='image1', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-3.0)
    image2 = visual.ImageStim(
        win=win,
        name='image2', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=512.0, interpolate=True, depth=-4.0)
    rating_1 = visual.Slider(win=win, name='rating_1',
        startValue=None, size=1.0, pos=[0,0], units=win.units,
        labels=("1 - not similar at all", 2, 3, 4, 5, 6, "7 - very similar"), ticks=(1, 2, 3, 4, 5, 6, 7), granularity=1.0,
        style='rating', styleTweaks=('triangleMarker',), opacity=None,
        labelColor='white', markerColor='black', lineColor='White', colorSpace='rgb',
        font='Verdana', labelHeight=0.35*y_scale,
        flip=False, ori=0.0, depth=-5, readOnly=False)
    question = visual.TextStim(win=win, name='question',
        text='Content',
        font='Verdana',
        pos=[0,0], height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-6.0);
    rating_2 = visual.Slider(win=win, name='rating_2',
        startValue=None, size=1.0, pos=[0,0], units=win.units,
        labels=("1 - not similar at all", 2, 3, 4, 5, 6, "7 - very similar"), ticks=(1, 2, 3, 4, 5, 6, 7), granularity=1.0,
        style='rating', styleTweaks=('triangleMarker',), opacity=None,
        labelColor='white', markerColor='black', lineColor='White', colorSpace='rgb',
        font='Verdana', labelHeight=0.35*y_scale,
        flip=False, ori=0.0, depth=-7, readOnly=False)
    question_2 = visual.TextStim(win=win, name='question_2',
        text='Drawing style/ability',
        font='Verdana',
        pos=[0,0], height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-8.0);
    block_text = visual.TextStim(win=win, name='block_text',
        text='',
        font='Open Sans',
        pos=(16*x_scale, 9.8*y_scale), height=0.02, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-10.0);
    
    # --- Initialize components for Routine "end" ---
    message_end = visual.TextStim(win=win, name='message_end',
        text='',
        font='Verdana',
        pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    key_resp = keyboard.Keyboard()
    
    # create some handy timers
    if globalClock is None:
        globalClock = core.Clock()  # to track the time since experiment started
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)
    
    # --- Prepare to start Routine "preparation" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('preparation.started', globalClock.getTime())
    # keep track of which components have finished
    preparationComponents = []
    for thisComponent in preparationComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "preparation" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in preparationComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "preparation" ---
    for thisComponent in preparationComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('preparation.stopped', globalClock.getTime())
    # Run 'End Routine' code from condition_matrix
    # Specify the file path based on the 'session'
    trial_condition_path = "conditions/conditions_" + expInfo['session'] + ".csv"
    practice_file_path = "conditions/practice_cond_" + expInfo['session'] + ".csv"
    instructions_file_path = "instructions/instructions_" + expInfo['session'] + ".xlsx"
    
    # Set the random seed based on the participant number
    seed_num = int(expInfo['participant'])
    # the Routine "preparation" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    instructionsLoop = data.TrialHandler(nReps=1.0, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(instructions_file_path),
        seed=None, name='instructionsLoop')
    thisExp.addLoop(instructionsLoop)  # add the loop to the experiment
    thisInstructionsLoop = instructionsLoop.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisInstructionsLoop.rgb)
    if thisInstructionsLoop != None:
        for paramName in thisInstructionsLoop:
            globals()[paramName] = thisInstructionsLoop[paramName]
    
    for thisInstructionsLoop in instructionsLoop:
        currentLoop = instructionsLoop
        thisExp.timestampOnFlip(win, 'thisRow.t')
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                inputs=inputs, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
        )
        # abbreviate parameter names if possible (e.g. rgb = thisInstructionsLoop.rgb)
        if thisInstructionsLoop != None:
            for paramName in thisInstructionsLoop:
                globals()[paramName] = thisInstructionsLoop[paramName]
        
        # --- Prepare to start Routine "instructions" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('instructions.started', globalClock.getTime())
        startResponse.keys = []
        startResponse.rt = []
        _startResponse_allKeys = []
        slide_img.setImage(slide)
        # keep track of which components have finished
        instructionsComponents = [startResponse, slide_img]
        for thisComponent in instructionsComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "instructions" ---
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *startResponse* updates
            waitOnFlip = False
            
            # if startResponse is starting this frame...
            if startResponse.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                startResponse.frameNStart = frameN  # exact frame index
                startResponse.tStart = t  # local t and not account for scr refresh
                startResponse.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(startResponse, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'startResponse.started')
                # update status
                startResponse.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(startResponse.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(startResponse.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if startResponse.status == STARTED and not waitOnFlip:
                theseKeys = startResponse.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                _startResponse_allKeys.extend(theseKeys)
                if len(_startResponse_allKeys):
                    startResponse.keys = _startResponse_allKeys[-1].name  # just the last key pressed
                    startResponse.rt = _startResponse_allKeys[-1].rt
                    startResponse.duration = _startResponse_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # *slide_img* updates
            
            # if slide_img is starting this frame...
            if slide_img.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                slide_img.frameNStart = frameN  # exact frame index
                slide_img.tStart = t  # local t and not account for scr refresh
                slide_img.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(slide_img, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'slide_img.started')
                # update status
                slide_img.status = STARTED
                slide_img.setAutoDraw(True)
            
            # if slide_img is active this frame...
            if slide_img.status == STARTED:
                # update params
                pass
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in instructionsComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "instructions" ---
        for thisComponent in instructionsComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('instructions.stopped', globalClock.getTime())
        # check responses
        if startResponse.keys in ['', [], None]:  # No response was made
            startResponse.keys = None
        instructionsLoop.addData('startResponse.keys',startResponse.keys)
        if startResponse.keys != None:  # we had a response
            instructionsLoop.addData('startResponse.rt', startResponse.rt)
            instructionsLoop.addData('startResponse.duration', startResponse.duration)
        # the Routine "instructions" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
    # completed 1.0 repeats of 'instructionsLoop'
    
    
    # set up handler to look after randomisation of conditions etc
    practice = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(practice_file_path),
        seed=seed_num, name='practice')
    thisExp.addLoop(practice)  # add the loop to the experiment
    thisPractice = practice.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisPractice.rgb)
    if thisPractice != None:
        for paramName in thisPractice:
            globals()[paramName] = thisPractice[paramName]
    
    for thisPractice in practice:
        currentLoop = practice
        thisExp.timestampOnFlip(win, 'thisRow.t')
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                inputs=inputs, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
        )
        # abbreviate parameter names if possible (e.g. rgb = thisPractice.rgb)
        if thisPractice != None:
            for paramName in thisPractice:
                globals()[paramName] = thisPractice[paramName]
        
        # --- Prepare to start Routine "trial" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('trial.started', globalClock.getTime())
        # Run 'Begin Routine' code from positions
        if random() > 0.5:
            img1pos = -img_x_pos
            img2pos = img_x_pos
        else:
            img1pos = img_x_pos
            img2pos = -img_x_pos
        image1.setPos((img1pos*x_scale, 4*y_scale))
        image1.setSize((8*x_scale, 8*x_scale*0.75))
        image1.setImage(stim1)
        image2.setPos((img2pos*x_scale, 4*y_scale))
        image2.setSize((8*x_scale, 8*x_scale*0.75))
        image2.setImage(stim2)
        rating_1.reset()
        rating_1.setPos((0, -3.5*y_scale))
        rating_1.setSize((12*x_scale, 0.75*y_scale))
        question.setPos((0, -2*y_scale))
        question.setHeight(0.7*y_scale)
        rating_2.reset()
        rating_2.setPos((0, -7*y_scale))
        rating_2.setSize((12*x_scale, 0.75*y_scale))
        question_2.setPos((0, -5.5*y_scale))
        question_2.setHeight(0.7*y_scale)
        # Run 'Begin Routine' code from make_block_text
        block_message = 'Block ' + str(blockN-1) + ' of ' + str(maxBlock)
        block_text.setText(block_message)
        # keep track of which components have finished
        trialComponents = [polygon, polygon_2, image1, image2, rating_1, question, rating_2, question_2, block_text]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "trial" ---
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *polygon* updates
            
            # if polygon is starting this frame...
            if polygon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon.frameNStart = frameN  # exact frame index
                polygon.tStart = t  # local t and not account for scr refresh
                polygon.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon.started')
                # update status
                polygon.status = STARTED
                polygon.setAutoDraw(True)
            
            # if polygon is active this frame...
            if polygon.status == STARTED:
                # update params
                pass
            
            # *polygon_2* updates
            
            # if polygon_2 is starting this frame...
            if polygon_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon_2.frameNStart = frameN  # exact frame index
                polygon_2.tStart = t  # local t and not account for scr refresh
                polygon_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_2.started')
                # update status
                polygon_2.status = STARTED
                polygon_2.setAutoDraw(True)
            
            # if polygon_2 is active this frame...
            if polygon_2.status == STARTED:
                # update params
                pass
            
            # *image1* updates
            
            # if image1 is starting this frame...
            if image1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image1.frameNStart = frameN  # exact frame index
                image1.tStart = t  # local t and not account for scr refresh
                image1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'image1.started')
                # update status
                image1.status = STARTED
                image1.setAutoDraw(True)
            
            # if image1 is active this frame...
            if image1.status == STARTED:
                # update params
                pass
            
            # *image2* updates
            
            # if image2 is starting this frame...
            if image2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image2.frameNStart = frameN  # exact frame index
                image2.tStart = t  # local t and not account for scr refresh
                image2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'image2.started')
                # update status
                image2.status = STARTED
                image2.setAutoDraw(True)
            
            # if image2 is active this frame...
            if image2.status == STARTED:
                # update params
                pass
            
            # *rating_1* updates
            
            # if rating_1 is starting this frame...
            if rating_1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                rating_1.frameNStart = frameN  # exact frame index
                rating_1.tStart = t  # local t and not account for scr refresh
                rating_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(rating_1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'rating_1.started')
                # update status
                rating_1.status = STARTED
                rating_1.setAutoDraw(True)
            
            # if rating_1 is active this frame...
            if rating_1.status == STARTED:
                # update params
                pass
            
            # *question* updates
            
            # if question is starting this frame...
            if question.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                question.frameNStart = frameN  # exact frame index
                question.tStart = t  # local t and not account for scr refresh
                question.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(question, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'question.started')
                # update status
                question.status = STARTED
                question.setAutoDraw(True)
            
            # if question is active this frame...
            if question.status == STARTED:
                # update params
                pass
            
            # *rating_2* updates
            
            # if rating_2 is starting this frame...
            if rating_2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                rating_2.frameNStart = frameN  # exact frame index
                rating_2.tStart = t  # local t and not account for scr refresh
                rating_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(rating_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'rating_2.started')
                # update status
                rating_2.status = STARTED
                rating_2.setAutoDraw(True)
            
            # if rating_2 is active this frame...
            if rating_2.status == STARTED:
                # update params
                pass
            
            # *question_2* updates
            
            # if question_2 is starting this frame...
            if question_2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                question_2.frameNStart = frameN  # exact frame index
                question_2.tStart = t  # local t and not account for scr refresh
                question_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(question_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'question_2.started')
                # update status
                question_2.status = STARTED
                question_2.setAutoDraw(True)
            
            # if question_2 is active this frame...
            if question_2.status == STARTED:
                # update params
                pass
            
            # *block_text* updates
            
            # if block_text is starting this frame...
            if block_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_text.frameNStart = frameN  # exact frame index
                block_text.tStart = t  # local t and not account for scr refresh
                block_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'block_text.started')
                # update status
                block_text.status = STARTED
                block_text.setAutoDraw(True)
            
            # if block_text is active this frame...
            if block_text.status == STARTED:
                # update params
                pass
            # Run 'Each Frame' code from wait_for_both_sliders
            keys = event.getKeys()
            
            # need to manually check for escape, as our key checking will interfere with
            # Builder's escape check:
            if 'escape' in keys: 
                core.quit()
            
            if 'space' in keys:
                # check how many ratings have been completed:
                completed_ratings = 0
                for scale in [rating_1, rating_2]:
                    if scale.getRating() is not None:
                        completed_ratings = completed_ratings + 1
            
                if completed_ratings == 2:
                    continueRoutine = False # end now
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "trial" ---
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('trial.stopped', globalClock.getTime())
        # Run 'End Routine' code from positions
        thisExp.addData('stim1pos', img1pos)
        thisExp.addData('stim2pos', img2pos)
        thisExp.addData('block', blockN-1)
        
        practice.addData('rating_1.response', rating_1.getRating())
        practice.addData('rating_1.rt', rating_1.getRT())
        practice.addData('rating_2.response', rating_2.getRating())
        practice.addData('rating_2.rt', rating_2.getRT())
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "ITI" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('ITI.started', globalClock.getTime())
        # keep track of which components have finished
        ITIComponents = [text]
        for thisComponent in ITIComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "ITI" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 0.3:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *text* updates
            
            # if text is starting this frame...
            if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text.frameNStart = frameN  # exact frame index
                text.tStart = t  # local t and not account for scr refresh
                text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'text.started')
                # update status
                text.status = STARTED
                text.setAutoDraw(True)
            
            # if text is active this frame...
            if text.status == STARTED:
                # update params
                pass
            
            # if text is stopping this frame...
            if text.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > text.tStartRefresh + 0.3-frameTolerance:
                    # keep track of stop time/frame for later
                    text.tStop = t  # not accounting for scr refresh
                    text.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'text.stopped')
                    # update status
                    text.status = FINISHED
                    text.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in ITIComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "ITI" ---
        for thisComponent in ITIComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('ITI.stopped', globalClock.getTime())
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.300000)
        thisExp.nextEntry()
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
    # completed 1.0 repeats of 'practice'
    
    
    # --- Prepare to start Routine "practice_over" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('practice_over.started', globalClock.getTime())
    space_practice_over.keys = []
    space_practice_over.rt = []
    _space_practice_over_allKeys = []
    # keep track of which components have finished
    practice_overComponents = [text_2, space_practice_over]
    for thisComponent in practice_overComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "practice_over" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text_2* updates
        
        # if text_2 is starting this frame...
        if text_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_2.frameNStart = frameN  # exact frame index
            text_2.tStart = t  # local t and not account for scr refresh
            text_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_2.started')
            # update status
            text_2.status = STARTED
            text_2.setAutoDraw(True)
        
        # if text_2 is active this frame...
        if text_2.status == STARTED:
            # update params
            pass
        
        # *space_practice_over* updates
        waitOnFlip = False
        
        # if space_practice_over is starting this frame...
        if space_practice_over.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            space_practice_over.frameNStart = frameN  # exact frame index
            space_practice_over.tStart = t  # local t and not account for scr refresh
            space_practice_over.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(space_practice_over, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'space_practice_over.started')
            # update status
            space_practice_over.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(space_practice_over.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(space_practice_over.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if space_practice_over.status == STARTED and not waitOnFlip:
            theseKeys = space_practice_over.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
            _space_practice_over_allKeys.extend(theseKeys)
            if len(_space_practice_over_allKeys):
                space_practice_over.keys = _space_practice_over_allKeys[-1].name  # just the last key pressed
                space_practice_over.rt = _space_practice_over_allKeys[-1].rt
                space_practice_over.duration = _space_practice_over_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in practice_overComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "practice_over" ---
    for thisComponent in practice_overComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('practice_over.stopped', globalClock.getTime())
    # check responses
    if space_practice_over.keys in ['', [], None]:  # No response was made
        space_practice_over.keys = None
    thisExp.addData('space_practice_over.keys',space_practice_over.keys)
    if space_practice_over.keys != None:  # we had a response
        thisExp.addData('space_practice_over.rt', space_practice_over.rt)
        thisExp.addData('space_practice_over.duration', space_practice_over.duration)
    thisExp.nextEntry()
    # the Routine "practice_over" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(trial_condition_path),
        seed=seed_num, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            globals()[paramName] = thisTrial[paramName]
    
    for thisTrial in trials:
        currentLoop = trials
        thisExp.timestampOnFlip(win, 'thisRow.t')
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                inputs=inputs, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
        )
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                globals()[paramName] = thisTrial[paramName]
        
        # --- Prepare to start Routine "breakBlocks" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('breakBlocks.started', globalClock.getTime())
        # Run 'Begin Routine' code from break_2
        break_instructions = visual.TextStim(
            win = win,
            text = f"Great job! You can take a break now. \n\n Try to stay focused! \n\n Press space to start block {blockN} of {maxBlock}.",
            height = 0.05,
            color = "black")
        
        # do not show if block is not over yet
        if trials.thisN > 0 and trials.thisN % 63 != 0:
            continueRoutine = False
        else:
            # do not show before first block
            if blockN == 1:
                blockN = blockN + 1
                continueRoutine = False
            else:
                break_instructions.setAutoDraw(True)
                blockN = blockN + 1
        
        startBlock.keys = []
        startBlock.rt = []
        _startBlock_allKeys = []
        # keep track of which components have finished
        breakBlocksComponents = [startBlock]
        for thisComponent in breakBlocksComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "breakBlocks" ---
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *startBlock* updates
            waitOnFlip = False
            
            # if startBlock is starting this frame...
            if startBlock.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                startBlock.frameNStart = frameN  # exact frame index
                startBlock.tStart = t  # local t and not account for scr refresh
                startBlock.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(startBlock, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'startBlock.started')
                # update status
                startBlock.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(startBlock.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(startBlock.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if startBlock.status == STARTED and not waitOnFlip:
                theseKeys = startBlock.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                _startBlock_allKeys.extend(theseKeys)
                if len(_startBlock_allKeys):
                    startBlock.keys = _startBlock_allKeys[-1].name  # just the last key pressed
                    startBlock.rt = _startBlock_allKeys[-1].rt
                    startBlock.duration = _startBlock_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in breakBlocksComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "breakBlocks" ---
        for thisComponent in breakBlocksComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('breakBlocks.stopped', globalClock.getTime())
        # Run 'End Routine' code from break_2
        break_instructions.setAutoDraw(False)
        # check responses
        if startBlock.keys in ['', [], None]:  # No response was made
            startBlock.keys = None
        trials.addData('startBlock.keys',startBlock.keys)
        if startBlock.keys != None:  # we had a response
            trials.addData('startBlock.rt', startBlock.rt)
            trials.addData('startBlock.duration', startBlock.duration)
        # the Routine "breakBlocks" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "ITI" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('ITI.started', globalClock.getTime())
        # keep track of which components have finished
        ITIComponents = [text]
        for thisComponent in ITIComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "ITI" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 0.3:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *text* updates
            
            # if text is starting this frame...
            if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text.frameNStart = frameN  # exact frame index
                text.tStart = t  # local t and not account for scr refresh
                text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'text.started')
                # update status
                text.status = STARTED
                text.setAutoDraw(True)
            
            # if text is active this frame...
            if text.status == STARTED:
                # update params
                pass
            
            # if text is stopping this frame...
            if text.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > text.tStartRefresh + 0.3-frameTolerance:
                    # keep track of stop time/frame for later
                    text.tStop = t  # not accounting for scr refresh
                    text.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'text.stopped')
                    # update status
                    text.status = FINISHED
                    text.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in ITIComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "ITI" ---
        for thisComponent in ITIComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('ITI.stopped', globalClock.getTime())
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.300000)
        
        # --- Prepare to start Routine "trial" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('trial.started', globalClock.getTime())
        # Run 'Begin Routine' code from positions
        if random() > 0.5:
            img1pos = -img_x_pos
            img2pos = img_x_pos
        else:
            img1pos = img_x_pos
            img2pos = -img_x_pos
        image1.setPos((img1pos*x_scale, 4*y_scale))
        image1.setSize((8*x_scale, 8*x_scale*0.75))
        image1.setImage(stim1)
        image2.setPos((img2pos*x_scale, 4*y_scale))
        image2.setSize((8*x_scale, 8*x_scale*0.75))
        image2.setImage(stim2)
        rating_1.reset()
        rating_1.setPos((0, -3.5*y_scale))
        rating_1.setSize((12*x_scale, 0.75*y_scale))
        question.setPos((0, -2*y_scale))
        question.setHeight(0.7*y_scale)
        rating_2.reset()
        rating_2.setPos((0, -7*y_scale))
        rating_2.setSize((12*x_scale, 0.75*y_scale))
        question_2.setPos((0, -5.5*y_scale))
        question_2.setHeight(0.7*y_scale)
        # Run 'Begin Routine' code from make_block_text
        block_message = 'Block ' + str(blockN-1) + ' of ' + str(maxBlock)
        block_text.setText(block_message)
        # keep track of which components have finished
        trialComponents = [polygon, polygon_2, image1, image2, rating_1, question, rating_2, question_2, block_text]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "trial" ---
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *polygon* updates
            
            # if polygon is starting this frame...
            if polygon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon.frameNStart = frameN  # exact frame index
                polygon.tStart = t  # local t and not account for scr refresh
                polygon.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon.started')
                # update status
                polygon.status = STARTED
                polygon.setAutoDraw(True)
            
            # if polygon is active this frame...
            if polygon.status == STARTED:
                # update params
                pass
            
            # *polygon_2* updates
            
            # if polygon_2 is starting this frame...
            if polygon_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon_2.frameNStart = frameN  # exact frame index
                polygon_2.tStart = t  # local t and not account for scr refresh
                polygon_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_2.started')
                # update status
                polygon_2.status = STARTED
                polygon_2.setAutoDraw(True)
            
            # if polygon_2 is active this frame...
            if polygon_2.status == STARTED:
                # update params
                pass
            
            # *image1* updates
            
            # if image1 is starting this frame...
            if image1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image1.frameNStart = frameN  # exact frame index
                image1.tStart = t  # local t and not account for scr refresh
                image1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'image1.started')
                # update status
                image1.status = STARTED
                image1.setAutoDraw(True)
            
            # if image1 is active this frame...
            if image1.status == STARTED:
                # update params
                pass
            
            # *image2* updates
            
            # if image2 is starting this frame...
            if image2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image2.frameNStart = frameN  # exact frame index
                image2.tStart = t  # local t and not account for scr refresh
                image2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'image2.started')
                # update status
                image2.status = STARTED
                image2.setAutoDraw(True)
            
            # if image2 is active this frame...
            if image2.status == STARTED:
                # update params
                pass
            
            # *rating_1* updates
            
            # if rating_1 is starting this frame...
            if rating_1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                rating_1.frameNStart = frameN  # exact frame index
                rating_1.tStart = t  # local t and not account for scr refresh
                rating_1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(rating_1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'rating_1.started')
                # update status
                rating_1.status = STARTED
                rating_1.setAutoDraw(True)
            
            # if rating_1 is active this frame...
            if rating_1.status == STARTED:
                # update params
                pass
            
            # *question* updates
            
            # if question is starting this frame...
            if question.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                question.frameNStart = frameN  # exact frame index
                question.tStart = t  # local t and not account for scr refresh
                question.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(question, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'question.started')
                # update status
                question.status = STARTED
                question.setAutoDraw(True)
            
            # if question is active this frame...
            if question.status == STARTED:
                # update params
                pass
            
            # *rating_2* updates
            
            # if rating_2 is starting this frame...
            if rating_2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                rating_2.frameNStart = frameN  # exact frame index
                rating_2.tStart = t  # local t and not account for scr refresh
                rating_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(rating_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'rating_2.started')
                # update status
                rating_2.status = STARTED
                rating_2.setAutoDraw(True)
            
            # if rating_2 is active this frame...
            if rating_2.status == STARTED:
                # update params
                pass
            
            # *question_2* updates
            
            # if question_2 is starting this frame...
            if question_2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                question_2.frameNStart = frameN  # exact frame index
                question_2.tStart = t  # local t and not account for scr refresh
                question_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(question_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'question_2.started')
                # update status
                question_2.status = STARTED
                question_2.setAutoDraw(True)
            
            # if question_2 is active this frame...
            if question_2.status == STARTED:
                # update params
                pass
            
            # *block_text* updates
            
            # if block_text is starting this frame...
            if block_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_text.frameNStart = frameN  # exact frame index
                block_text.tStart = t  # local t and not account for scr refresh
                block_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'block_text.started')
                # update status
                block_text.status = STARTED
                block_text.setAutoDraw(True)
            
            # if block_text is active this frame...
            if block_text.status == STARTED:
                # update params
                pass
            # Run 'Each Frame' code from wait_for_both_sliders
            keys = event.getKeys()
            
            # need to manually check for escape, as our key checking will interfere with
            # Builder's escape check:
            if 'escape' in keys: 
                core.quit()
            
            if 'space' in keys:
                # check how many ratings have been completed:
                completed_ratings = 0
                for scale in [rating_1, rating_2]:
                    if scale.getRating() is not None:
                        completed_ratings = completed_ratings + 1
            
                if completed_ratings == 2:
                    continueRoutine = False # end now
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "trial" ---
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('trial.stopped', globalClock.getTime())
        # Run 'End Routine' code from positions
        thisExp.addData('stim1pos', img1pos)
        thisExp.addData('stim2pos', img2pos)
        thisExp.addData('block', blockN-1)
        
        trials.addData('rating_1.response', rating_1.getRating())
        trials.addData('rating_1.rt', rating_1.getRT())
        trials.addData('rating_2.response', rating_2.getRating())
        trials.addData('rating_2.rt', rating_2.getRT())
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
    # completed 1.0 repeats of 'trials'
    
    
    # --- Prepare to start Routine "end" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('end.started', globalClock.getTime())
    message_end.setText('The end!\n\nThank you for participating')
    key_resp.keys = []
    key_resp.rt = []
    _key_resp_allKeys = []
    # keep track of which components have finished
    endComponents = [message_end, key_resp]
    for thisComponent in endComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "end" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *message_end* updates
        
        # if message_end is starting this frame...
        if message_end.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            message_end.frameNStart = frameN  # exact frame index
            message_end.tStart = t  # local t and not account for scr refresh
            message_end.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(message_end, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'message_end.started')
            # update status
            message_end.status = STARTED
            message_end.setAutoDraw(True)
        
        # if message_end is active this frame...
        if message_end.status == STARTED:
            # update params
            pass
        
        # *key_resp* updates
        waitOnFlip = False
        
        # if key_resp is starting this frame...
        if key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            key_resp.frameNStart = frameN  # exact frame index
            key_resp.tStart = t  # local t and not account for scr refresh
            key_resp.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(key_resp, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'key_resp.started')
            # update status
            key_resp.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(key_resp.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(key_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if key_resp.status == STARTED and not waitOnFlip:
            theseKeys = key_resp.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
            _key_resp_allKeys.extend(theseKeys)
            if len(_key_resp_allKeys):
                key_resp.keys = _key_resp_allKeys[-1].name  # just the last key pressed
                key_resp.rt = _key_resp_allKeys[-1].rt
                key_resp.duration = _key_resp_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in endComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "end" ---
    for thisComponent in endComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('end.stopped', globalClock.getTime())
    # check responses
    if key_resp.keys in ['', [], None]:  # No response was made
        key_resp.keys = None
    thisExp.addData('key_resp.keys',key_resp.keys)
    if key_resp.keys != None:  # we had a response
        thisExp.addData('key_resp.rt', key_resp.rt)
        thisExp.addData('key_resp.duration', key_resp.duration)
    thisExp.nextEntry()
    # the Routine "end" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # mark experiment as finished
    endExperiment(thisExp, win=win, inputs=inputs)


def saveData(thisExp):
    """
    Save data from this experiment
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    filename = thisExp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)


def endExperiment(thisExp, inputs=None, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    inputs : dict
        Dictionary of input devices by name.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # mark experiment handler as finished
    thisExp.status = FINISHED
    # shut down eyetracker, if there is one
    if inputs is not None:
        if 'eyetracker' in inputs and inputs['eyetracker'] is not None:
            inputs['eyetracker'].setConnectionState(False)
    logging.flush()


def quit(thisExp, win=None, inputs=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
    inputs : dict
        Dictionary of input devices by name.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    thisExp.abort()  # or data files will save again on exit
    # make sure everything is closed down
    if win is not None:
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed before quitting
        win.flip()
        win.close()
    if inputs is not None:
        if 'eyetracker' in inputs and inputs['eyetracker'] is not None:
            inputs['eyetracker'].setConnectionState(False)
    logging.flush()
    if thisSession is not None:
        thisSession.stop()
    # terminate Python process
    core.quit()


# if running this experiment as a script...
if __name__ == '__main__':
    # call all functions in order
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    inputs = setupInputs(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win, 
        inputs=inputs
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win, inputs=inputs)
