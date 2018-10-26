file_base_name = 'board_';
try
    screens=Screen('Screens');
    screenNumber=max(screens);
    [win, rect] = Screen('OpenWindow', screenNumber);%, [0 0 800 450]);

    grab = Screen('OpenVideoCapture', win, 0);
    Screen('StartVideoCapture', grab, 60, 1);
    for i = 1:4
        pause(10)
        [tex, pts, nrdropped, imtext]=Screen('GetCapturedImage', win, grab, 1, [], 2);
        img = permute(imtext([3,2,1], :,:), [3,2,1]);
        imwrite(img, [file_base_name, num2str(i), '.jpg']);
        Screen('DrawTexture', win, tex);
        Screen('Flip', win)
        Screen('Close', tex);
        pause(1)
    end
    Screen('StopVideoCapture', grab);
    Screen('CloseVideoCapture', grab);
catch err
   sca
   rethrow(err)
   Screen('StopVideoCapture', grab);
   Screen('CloseVideoCapture', grab);
end