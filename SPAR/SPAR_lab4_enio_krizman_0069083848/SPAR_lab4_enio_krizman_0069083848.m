% Senzori, percepcija i aktuacija u robotici

% 4. laboratorijska vježba 2022./2023.

%% read images

image = imread('checkboard.png');
im1 = imread("traffic1.png");
im2 = imread("traffic2.png");

%% Zadatak 1

% Prije obrade pretvorite vrijednosti intenziteta u realne brojeve.
image = im2double(image);

% Za težinsku funkciju w(i, j) koristite Gaussovu masku 3x3.
% Gaussovu masku izraˇcunajte funkcijom fspecial( gaussian ).
w = fspecial('gaussian',[3 3]);

% Za računanje vertikalnih i horizontalnih usmjerenih derivacija koristite konvoluciju Sobelovim operatorom,
% fspecial("sobel") Za računanje konvolucije koristite funkciju imfilter. Usmjerene derivacije izračunajte za svaki piksel.
sobel = fspecial('sobel');
I_u = imfilter(image, sobel);
I_v = imfilter(image, sobel');


figure(2)
imshow(I_v)
title('vertikalni gradijent')

figure(3)
imshow(I_u)
title('Horizontalni gradijent')

k = 0.04;
n = size(I_u);
image_ = zeros(n(1) - 3, n(2) - 3);
C_h = zeros(n(1) - 3, n(2) - 3);

for u = 1:(480-3)
    for v = 1:(600-3)

        A = zeros(2);
        A_ = zeros(2);

        for i=1:3
            for j=1:3
                 A_ = w(i, j) * [(I_u(u+i,v+j))^2 I_u(u+i,v+j)*I_v(u+i,v+j);
                     I_u(u+i,v+j)*I_v(u+i,v+j) (I_v(u+i,v+j))^2];
                 A = A + A_;
            end
         end
         
         %Harrisova mjera
         C_h(u,v) = det(A) - k*trace(A)^2;

         if C_h(u,v) > 0.01
             %Za kuteve uzmite piksele ako je vrijednost Harrisove mjere veća od praga
             image_(u, v) = image(u, v);
         end           
    end
end

figure(4)
imshow(image_)
title('Odziv funkcije C_h')

figure(5)
imshow(image)
hold on
[x, y] = find(image_);
% Mark detected corners on the image
plot(y, x, 'ro');
title('Detektirani kutovi')
hold off;

%% Zadatak 2

im1 = imresize(img1, 0.5);
im2 = imresize(img2, 0.5);

ww = 45;
w = round(ww/2);

Ix_m = conv2(im1, [-1 1; -1 1], 'valid');
Iy_m = conv2(im1, [-1 -1; 1 1], 'valid');
It_m = conv2(im1, ones(2), 'valid') + conv2(im2, -ones(2), 'valid');

u = zeros(size(im1));
v = zeros(size(im2));

for i = w+1:size(Ix_m,1)-w
    for j = w+1:size(Ix_m,2)-w
        Ix = Ix_m(i-w:i+w, j-w:j+w);
        Iy = Iy_m(i-w:i+w, j-w:j+w);
        It = It_m(i-w:i+w, j-w:j+w);
        
        Ix = Ix(:);
        Iy = Iy(:);
        b = -It(:);
        
        A = [Ix Iy];
        nu = pinv(A)*b;
        
        u(i,j) = nu(1);
        v(i,j) = nu(2);
    end
end

% downsize u and v
u_deci = u(1:10:end, 1:10:end);
v_deci = v(1:10:end, 1:10:end);
% get coordinate for u and v in the original frame
[m, n] = size(img1);
[X,Y] = meshgrid(1:n, 1:m);
X_deci = X(1:20:end, 1:20:end);
Y_deci = Y(1:20:end, 1:20:end);

% Optički tok računajte samo za značajke koje detektirate već 
% implementiranim Harrisovim detektorom
corners = detectHarrisFeatures(img1);
figure(6);
imshow(img1);
hold on;
%Prag postavite tako da imate oko 30 detektiranih značajki
points = points.selectStrongest(30);
% Ne koristite značajke koje se nalaze
% za manje od 50 piksela od ruba slika
border = 50;
validPoints = points.Location(:, 1) > border & points.Location(:, 1) < size(image1, 2) - border & ...
              points.Location(:, 2) > border & points.Location(:, 2) < size(image1, 1) - border;
points = points(validPoints);

plot(points);
title("Lucas-Kanade optički tok")
quiver(X_deci, Y_deci, u_deci, v_deci, 'y');
hold off
