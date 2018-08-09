 %% setup keyboard or fMRI trigger detection:
    trigSet = 15;
    triggering = struct('kbNum', []);
    
    DEVICE_NAME = 'Dell KB216 Wired Keyboard';
    [index, devName] = GetKeyboardIndices;
    for device = 1:length(index)
        if strcmp(devName(device),DEVICE_NAME)
            triggering.kbNum = index(device);
        end
    end
    triggering.kbNum = 7; %test to see why dev isn't working.
    %% wait for TR burn-in to begin new block
    KbQueueCreate(triggering.kbNum);
    KbQueueStart(triggering.kbNum);
%     Screen('DrawText', win, 'Please press any key to begin the next block.', round(screen_dim1/2), round(screen_dim2/2));
%     Screen('Flip', win);
%     %pause;

    burnInCount = 0; 
    press_sets = cell(3, 1000);
    k_it = 1;
    while burnInCount < 5
        [keyIsDown,keyCode] = KbQueueCheck(triggering.kbNum);
         keyPressed = find(keyCode);
        press_sets{1, k_it} = keyIsDown;
        press_sets{2, k_it} = keyCode;
        press_sets{3, k_it} = keyPressed;
        k_it = k_it + 1;
        
            if keyIsDown == 1 && ismember(keyPressed,trigSet)
                burnInCount = burnInCount + 1;
            end

         clear keyIsDown; clear keyCode; clear keyPressed;
         pause(.1);
    end

    fprintf('Burn in time complete. Starting experiment. \n')