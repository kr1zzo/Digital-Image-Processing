close all;

imgs = imageDatastore("D:\Documents\Fer\Senzori, percepcija i aktuacija u robotici\Laboratorijske vježbe\lab3\Slike"); %Podesite svoj path do kalibracijskih slika. U tom folderu moraju biti samo slike, bez drugih tipova datoteka!
[imagePoints,boardSize,~] = detectCheckerboardPoints(imgs.Files);

for i = 1:4
  I = imread(imgs.Files{i});
  subplot(2, 2, i);
  imshow(I);
  hold on;
  plot(imagePoints(:,1,i),imagePoints(:,2,i),'ro');
end

squareSize = 0.022;     %Ovdje postavite izmjerenu širinu/visinu kvadrata u metrima.
[worldPoints] = generateCheckerboardPoints(boardSize, squareSize);

[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints);

figure; 
imshow(imgs.Files{1}); 
hold on;
plot(imagePoints(:,1,1), imagePoints(:,2,1),'go');
plot(cameraParams.ReprojectedPoints(:,1,1),cameraParams.ReprojectedPoints(:,2,1),'r+');
legend('Detected Points','ReprojectedPoints');
hold off;



figure;
showExtrinsics(cameraParams,'patternCentric');

K = cameraParams.IntrinsicMatrix.'
