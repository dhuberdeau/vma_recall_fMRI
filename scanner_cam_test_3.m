%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
AssertOpenGL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

test_location = 'SSS'; %options: SSS, BIC, MRRC

switch test_location
    case 'SSS'
        screen_dims = [1920, 1080];
        res1 = 1920;
        res2 = 1080;
        DISC_SIZE = 40;
    case 'BIC'
        screen_dims = [1600 900]; % for MRRC camera
        res1 = 1280;%1920;
        res2 = 1024;%1080;
        
        DISC_SIZE = 32;%40 %40 for SSS camera, 32 for BIC scanner camera

    case 'MRRC'
        
    otherwise
        warning('test location unrecognized: using BIC configuration')
        screen_dims = [1600 900]; % for MRRC camera
        res1 = 1280;%1920;
        res2 = 1024;%1080;
        
        DISC_SIZE = 32;%40 %40 for SSS camera, 32 for BIC scanner camera
end
ind1 = repmat((1:res2)', 1, res1);
ind2 = repmat((1:res1), res2, 1);

ind1_d = repmat((1:DISC_SIZE:res2)', 1, res1/DISC_SIZE);
ind2_d = repmat((1:DISC_SIZE:res1), res2/DISC_SIZE, 1);

x = nan(1, 10000);
y = nan(1, 10000);
tim = nan(1, 10000);
cursor_dims = [-10 -10 10 10]';

screens=Screen('Screens');
screenNumber=max(screens);
[win, rect] = Screen('OpenWindow', screenNumber, 0); %, [0 0 800 450]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
dev_list = Screen('VideoCaptureDevices');
grabber = Screen('OpenVideoCapture', win, dev_list(5).DeviceIndex);
Screen('StartVideoCapture', grabber, 60, 1);
RMIN = 0;
RMAX = .025;
SUBWIN_SIZE = 75;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
timr = tic;
keyIsDown = 0;
k_samp = 1;
while(~keyIsDown)
    [ keyIsDown, keyTime, keyCode ] = KbCheck; 
    [tex, pts, nrdropped, imtext] = Screen('GetCapturedImage', win, grabber, 1, [], 2);
    
%     img_ = imtext;
    img_ = imtext(:, 1:DISC_SIZE:end, 1:DISC_SIZE:end);
    img = permute(img_([3,2,1], :,:), [3,2,1]);
    b = rgb2hsv(img);

    switch test_location
        case 'MRRC'
            im_r = b(:,:,3) > REFL_TH; %For MRRC
        case 'BIC'
            im_r = inRange(b, [RMAX 1 1], [RMIN 0.5 0.5]); % for BIC or SSS cameras
        case 'SSS'
            im_r = inRange(b, [RMAX 1 1], [RMIN 0.5 0.5]); % for BIC or SSS cameras
        otherwise
            im_r = inRange(b, [RMAX 1 1], [RMIN 0.5 0.5]); % for BIC or SSS cameras
    end

    trk_y_rd = (median(ind1_d(im_r)));
    trk_x_rd = (median(ind2_d(im_r)));
    if ~isempty(trk_y_rd) && ~isempty(trk_x_rd)
        img_ = imtext(:, max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]), max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]));
        img = permute(img_([3 2 1], :, :), [3 2 1]);
        c_r = rgb2hsv(img);
        
        switch test_location
            case 'SSS'
                im_r = inRange(c_r, [RMAX 1 1], [RMIN 0.5 0.5]);
            case 'BIC'
                im_r = inRange(c_r, [RMAX 1 1], [RMIN 0.5 0.5]);
            case 'MRRC'
                im_r = c_r(:,:,3) > REFL_TH;
            otherwise
                im_r = inRange(c_r, [RMAX 1 1], [RMIN 0.5 0.5]);
        end

        rel_ind2 = ind2(max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]),max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]));
        rel_ind1 = ind1(max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]),max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]));
        trk_y_r = median(rel_ind1(im_r))*screen_dims(1)/res1;
        trk_x_r = (res1 - median(rel_ind2(im_r)))*screen_dims(2)/res2;

        Screen('FillOval', win, [200;0;0], [trk_x_r trk_y_r trk_x_r trk_y_r]' + cursor_dims);
        Screen('Flip', win);
        Screen('Close', tex);

        x(k_samp) = trk_x_r;
        y(k_samp) = trk_y_r;
        tim(k_samp) = toc(timr);
        k_samp = k_samp + 1;
    else
        Screen('Flip', win);
        Screen('Close', tex);

        x(k_samp) = nan;
        y(k_samp) = nan;
        tim(k_samp) = toc(timr);
        k_samp = k_samp + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
if exist('grabber')
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
