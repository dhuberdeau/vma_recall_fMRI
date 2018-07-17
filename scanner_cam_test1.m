%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
AssertOpenGL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

screens=Screen('Screens');
screenNumber=min(screens);
[win, rect] = Screen('OpenWindow', screenNumber, [0 0 800 450]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
dev_list = Screen('VideoCaptureDevices');
grabber = Screen('OpenVideoCapture', win, dev_list(5).DeviceIndex);
Screen('StartVideoCapture', grabber, 60, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
keyIsDown = 0;

while(~keyIsDown)
    [ keyIsDown, keyTime, keyCode ] = KbCheck; 
    [tex, pts, nrdropped, imtext] = Screen('GetCapturedImage', win, grabber, 1, [], 2);
    Screen('DrawTexture', win, tex);
    Screen('Flip', win)
    Screen('Close', tex);
    pause(.017);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
if exist('grabber')
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE