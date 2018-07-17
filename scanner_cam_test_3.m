%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
AssertOpenGL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
screen_dims = [1920, 1080];
res1 = 1280;%1920;
res2 = 1024;%1080;
ind1 = repmat((1:res2)', 1, res1);
ind2 = repmat((1:res1), res2, 1);
DISC_SIZE = 1;

ind1_d = repmat((1:DISC_SIZE:res2)', 1, res1/DISC_SIZE);
ind2_d = repmat((1:DISC_SIZE:res1), res2/DISC_SIZE, 1);
% SUBWIN_SIZE = 75;
% ind1 = repmat((1:res2)', 1, res1);
% ind2 = repmat((1:res1), res2, 1);

x = nan(1, 10000);
y = nan(1, 10000);
tim = nan(1, 10000);
cursor_dims = [-10 -10 10 10]';

screens=Screen('Screens');
screenNumber=min(screens);
[win, rect] = Screen('OpenWindow', screenNumber, 0); %, [0 0 800 450]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
dev_list = Screen('VideoCaptureDevices');
grabber = Screen('OpenVideoCapture', win, dev_list(5).DeviceIndex);
Screen('StartVideoCapture', grabber, 60, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
timr = tic;
keyIsDown = 0;
k_samp = 1;
while(~keyIsDown)
    [ keyIsDown, keyTime, keyCode ] = KbCheck; 
    [tex, pts, nrdropped, imtext] = Screen('GetCapturedImage', win, grabber, 1, [], 2);
    
    img_ = imtext;
    img = permute(img_([3,2,1], :,:), [3,2,1]);
    b = rgb2hsv(img);

    im_r = b(:,:,3) > .75;

    trk_y_rd = (median(ind1_d(im_r)));
    trk_x_rd = (median(ind2_d(im_r)));
    
%     img_ = imtext(:, max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]), max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]));
%     img = permute(img_([3 2 1], :, :), [3 2 1]);
%     %             c_r = rgb2hsv(img(max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]), max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]), :));
%     c_r = rgb2hsv(img);
%     %                 im_r = inRange(c_r, [.02 1 1], [0 0.5 0.5]);
%     im_r = c_r(:,:,3) > .75;
%     rel_ind2 = ind2(max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]),max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]));
%     rel_ind1 = ind1(max([(trk_y_rd - SUBWIN_SIZE),1]):min([(trk_y_rd + SUBWIN_SIZE),res2]),max([(trk_x_rd - SUBWIN_SIZE),1]):min([(trk_x_rd + SUBWIN_SIZE), res1]));
%     trk_y_r = median(rel_ind1(im_r))*screen_dims(1)/res1;
%     trk_x_r = median(rel_ind2(im_r))*screen_dims(2)/res2;

    trk_x_r = (res1 - trk_x_rd)*screen_dims(1)/res1;
    trk_y_r = trk_y_rd*screen_dims(2)/res2;
    

%     Screen('DrawTexture', win, tex);
    Screen('FillOval', win, [200;0;0], [trk_x_r trk_y_r trk_x_r trk_y_r]' + cursor_dims);
    Screen('Flip', win);
    Screen('Close', tex);
    
    x(k_samp) = trk_x_r;
    y(k_samp) = trk_y_r;
    tim(k_samp) = toc(timr);
    k_samp = k_samp + 1;
%     pause(.005);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
if exist('grabber')
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
