% Use this script at the beginning of each experiment to calibrate any
% potential rotation of the camera lens relative to the participant's
% workspace.

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

cursor_dims = [-10 -10 10 10]';

load('camera_params.mat');
load('mm_per_pix.mat');

screens=Screen('Screens');
screenNumber=max(screens);
[win, rect] = Screen('OpenWindow', screenNumber, 0, [0 0 800 450]);% 0); %, [0 0 800 450]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
dev_list = Screen('VideoCaptureDevices');
grabber = Screen('OpenVideoCapture', win, dev_list(5).DeviceIndex);
Screen('StartVideoCapture', grabber, 60, 1);

RMIN = [-.05, .9, .4];
RMAX = [.05, 1.1, .75];
GMIN = [.3, .9, .2]; %[.3 .95 .2]
GMAX = [.4, 1.1, .4]; %[.8 1.1 .5]
BMIN = [.4 .9 .5];
BMAX = [.6 1.1 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA KAPTURE
N_PICs = 10;
x_mean = nan(2, N_PICs);
y_mean = nan(2, N_PICs);
windw = 75;
try
    for i_pic = 1:size(x_mean,2)
        % grab a snapshot from camera
        [tex, pts, nrdropped, imtext] = Screen('GetCapturedImage', win, grabber, 1, [], 2);

%         disp('checkpoint: A');
        
        % convert to HSV
        img = permute(imtext([3,2,1], :,:), [3,2,1]);
        img = undistortImage(img, camera_params);
        b = rgb2hsv(img);
%         disp('checkpoint: B');
        
        % find pixels with specified color
        im_g_ = inRange(b, GMAX, GMIN);
        im_b_ = inRange(b, BMAX, BMIN);
%         disp('checkpoint: C');
        
        % find location of colored markers and zero-out all other pixels
        trk_y_b_ = (median(ind1(im_b_)));
        trk_x_b_ = (median(ind2(im_b_)));
        im_b_(ind1(:,1) > (trk_y_b_ + windw) | ind1(:,1) < (trk_y_b_ - windw),...
            ind2(1,:) > (trk_x_b_ + windw) | ind2(1,:) < (trk_x_b_ - windw)) = 0;

        trk_y_g_ = (median(ind1(im_g_)));
        trk_x_g_ = (median(ind2(im_g_)));
        im_g_(ind1(:,1) > (trk_y_g_ + windw) | ind1(:,1) < (trk_y_g_ - windw),...
            ind2(1,:) > (trk_x_g_ + windw) | ind2(1,:) < (trk_x_g_ - windw)) = 0;
%         disp('checkpoint: D');
        
        % get geometric properties of region defined by colored marker
        stat_b = regionprops(im_b_, 'BoundingBox');
        stat_g = regionprops(im_g_, 'BoundingBox');
%         disp('checkpoint: E');
        
        % Get bounding box of marker regions
        boxs_b = nan(4,length(stat_b));
        box_size_b = nan(1, length(stat_b));
        for i_r = 1:length(stat_b)
            boxs_b(:, i_r) = stat_b(i_r).BoundingBox';
            box_size_b(i_r) = norm([boxs_b(3,i_r); boxs_b(4, i_r)]);
        end

        boxs_g = nan(4,length(stat_g));
        box_size_g = nan(1, length(stat_g));disp('checkpoint: A');
        for i_r = 1:length(stat_g)
            boxs_g(:, i_r) = stat_g(i_r).BoundingBox';
            box_size_g(i_r) = norm([boxs_g(3,i_r); boxs_g(4, i_r)]);
        end
%         disp('checkpoint: F');
        
        [~, k_box_b] = max(box_size_b);
        [~, k_box_g] = max(box_size_g);

        % get coordinates of colored markers:
        trk_x_b = boxs_b(1, k_box_b) + .5*boxs_b(3, k_box_b);
        trk_y_b = boxs_b(2, k_box_b) + .5*boxs_b(4, k_box_b);

        trk_x_g = boxs_g(1, k_box_g) + .5*boxs_g(3, k_box_g);
        trk_y_g = boxs_g(2, k_box_g) + .5*boxs_g(4, k_box_g);

        % add location of this snapshot's colored marking to list if a mean is
        % desired
        x_mean(1, i_pic) = trk_x_b; 
        x_mean(2, i_pic) = trk_x_g;
        y_mean(1, i_pic) = trk_y_b;
        y_mean(2, i_pic) = trk_y_g;
%         disp('checkpoint: G');
    end
catch err_
   sca 
end
%%
% show image of workspace
tex_ = Screen('MakeTexture', win, img);
Screen('DrawTexture', win, tex_);

% show location of detected colored patches
Screen('FillOval', win, [200;200;0], mean([x_mean(1,:)', y_mean(1,:)', x_mean(1,:)', y_mean(1,:)'], 1)' + cursor_dims);
Screen('FillOval', win, [200;0;0], mean([x_mean(2,:)', y_mean(2,:)', x_mean(2,:)', y_mean(2,:)'], 1)' + cursor_dims);
Screen('Flip', win);
Screen('Close', tex);

Screen('StopVideoCapture', grabber);
Screen('CloseVideoCapture', grabber);
%%
keyIsDown = 0;
while(~keyIsDown)
    [keyIsDown, keyTime, keyCode ] = KbCheck; 
    pause(.1);
end

% get mean locations of colored patches
v_b = mean([x_mean(1,:)', y_mean(1,:)'],1)';
v_g = mean([x_mean(2,:)', y_mean(2,:)'],1)';

sca

%% 
x = v_b - v_g;

confirm_mm_conversion_fact = 20/x(1); %mm per pixel

try
    if abs(confirm_mm_conversion_fact - mm_pix) > .05
        warning('pixel-mm conversion may be inaccurate');
    end
catch
   warning('Unable to confirm pixel-to-mm conversion'); 
end

figure;
imshow(img); hold on;
plot(v_g(1) + [0 x(1)], v_g(2) + [0 x(2)], 'k-')
for i_ = 1:length(stat_b)
    plot(stat_b(i_).BoundingBox(1)+[0, stat_b(i_).BoundingBox(3)], stat_b(i_).BoundingBox(2) + [0, stat_b(i_).BoundingBox(4)], 'w-'); 
end
for i_ = 1:length(stat_g)
    plot(stat_g(i_).BoundingBox(1)+[0, stat_g(i_).BoundingBox(3)], stat_g(i_).BoundingBox(2) + [0, stat_g(i_).BoundingBox(4)], 'w-'); 
end

angle_error = atan2d(x(2), x(1));
mov_field_left_corner = v_g;
mov_field_right_corner = v_b;
desired_field_dims = [10, 20]./confirm_mm_conversion_fact;

save tracker_field_of_view_calibration angle_error mov_field_left_corner mov_field_right_corner desired_field_dims confirm_mm_conversion_fact

% calib_angle_error = atan2d(norm(cross([y(1), y(2), 0]', [x(1), x(2), 0]')), dot([y(1), y(2)]', [x(1) x(2)]'));

% save camera_angle_calibration angle_error calib_angle_error