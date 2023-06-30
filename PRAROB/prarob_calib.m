%% assignment 1 - get images
camera_info = imaqhwinfo('linuxvideo');
% use camera_info = imaqhwinfo('winvideo') on Windows machine

video_obj = videoinput('linuxvideo', camera_info.DeviceInfo.DeviceID, camera_info.DeviceInfo.DefaultFormat);

video_obj.ReturnedColorSpace = 'rgb';
f = figure('Visible', 'off'); 
vidRes = video_obj.VideoResolution;
imageRes = fliplr(vidRes);
hImage = imshow(zeros(imageRes));
num_images = 10;
for i=1:num_images
    % preview video object
    preview(video_obj, hImage);
    waitforbuttonpress;
    stoppreview(video_obj);
    % start to enable grab
    start(video_obj);
    img = getdata(video_obj);
    img = img(:, :, (1:3));
    imwrite(img, "calib_img"+num2str(i)+".png")
    stop(video_obj);
end

clear video_obj

%% Assignment 2
% After running the Camera Calibrator app, display the matrix of intrinsic
% parameters

cameraParams.Intrinsics.IntrinsicMatrix'

%% Assignment 3

% get one image
video_obj = videoinput('linuxvideo', camera_info.DeviceInfo.DeviceID, camera_info.DeviceInfo.DefaultFormat);
video_obj.ReturnedColorSpace = 'rgb';
f = figure('Visible', 'off'); 
vidRes = video_obj.VideoResolution;
imageRes = fliplr(vidRes);
hImage = imshow(zeros(imageRes));
% preview video object
preview(video_obj, hImage);
waitforbuttonpress;
stoppreview(video_obj);
% start to enable grab
start(video_obj);
img = getdata(video_obj);
img = img(:, :, (1:3));
stop(video_obj);
clear video_obj

%% process checkerboard

% grab an image 
video_obj = videoinput('linuxvideo', camera_info.DeviceInfo.DeviceID, camera_info.DeviceInfo.DefaultFormat);
video_obj.ReturnedColorSpace = 'rgb';
f = figure('Visible', 'off'); 
vidRes = video_obj.VideoResolution;
imageRes = fliplr(vidRes);
hImage = imshow(zeros(imageRes));
% preview video object
preview(video_obj, hImage);
waitforbuttonpress;
stoppreview(video_obj);
% start to enable grab
start(video_obj);
img = getdata(video_obj);
img = img(:, :, (1:3));
stop(video_obj);
clear video_obj

% using the instrinsic parameters we can undistort the image if necessary 
[img_undist,newOrigin] = undistortImage(img,cameraParams,'OutputView','full');

% calculate corners in the image
[imagePoints,boardSize] = detectCheckerboardPoints(img_undist);

% set size of a square in milimeters and calculate locations in the checkerboard frame
squareSize = 8;
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

% calculate extrinsic parameters
[R, t] = extrinsics(imagePoints, worldPoints, cameraParams)

[orientation_from_pattern, location_from_pattern] = extrinsicsToCameraPose(R, t)

%% Assignment 4
squareSize = 8;
checkerboard_points = generateCheckerboardPoints(boardSize, squareSize);

% create homogeneous coordinates
checkerboard_points_h = [checkerboard_points zeros(size(checkerboard_points,1),1) ones(size(checkerboard_points,1),1)];

% measure the pose of the pattern in the frame {R}
T = []%

worldPoints = transpose(T*transpose(checkerboard_points_h));
worldPoints = worldPoints(:,1:2);

[R, t] = extrinsics(imagePoints, worldPoints, cameraParams);

[orientation_from_world, location_from_world] = extrinsicsToCameraPose(R, t)

%% Assignment 5

% get one image where camera axis is vertical with respect to the pattern
video_obj = videoinput('linuxvideo', camera_info.DeviceInfo.DeviceID, camera_info.DeviceInfo.DefaultFormat);
video_obj.ReturnedColorSpace = 'rgb';
f = figure('Visible', 'off'); 
vidRes = video_obj.VideoResolution;
imageRes = fliplr(vidRes);
hImage = imshow(zeros(imageRes));
% preview video object
preview(video_obj, hImage);
waitforbuttonpress;
stoppreview(video_obj);
% start to enable grab
start(video_obj);
img = getdata(video_obj);
img = img(:, :, (1:3));
stop(video_obj);
clear video_obj

% insert processing to extract the colored circle from the image

% extract the centroid of the colored circle 

% using intrinsic camera matrix and distance of the pattern from the camera
% calculate the world coordinates of the colored circle

% z_pattern = % distance of the pattern from the camera
% u = % centroid x value in pixels
% v = % centroid y value in pixels
% cx, cy, fx, fy are instrinsic camera parameters, extract them from the
% matrix
% insert calculation
x_camera = [];
y_camera = [];
z_camera = [];

% convert to homogeneous coordinates
centroid_camera_h = [x_camera; y_camera; z_camera; 1];

% calculate centroid in world coordinate frame
T_world_camera = [orientation_from_world, location_from_world'; 0 0 0 1];

centroid_world_h = T_world_camera * centroid_camera_h;

centroid_world = centroid_world_h(1:3,1)
