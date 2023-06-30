% Senzori perecpija i aktuacija u robotici 
% 3. laboratorisjka vježba

%% 2. zadatak
clear all;
close all;
%% 2. a) + b)

%---z=1(a),z=2(b)-------------
z = 2;
%----------------------------

if z == 1
    load("measurementsNoisy4.mat");
    n = 4;
elseif z==2
    load("measurementsNoisy100.mat");
    n = 100;
end
x_i = X(1,:)
y_i = X(2,:)
x1
u1 = x1(1,:)
v1  = x1(2,:)

M1 = []
M2 = []
M3 = []

for i = 1:n
    M11 = [-X(1,i) -X(2,i) -1 0 0 0 X(1,i)*x1(1,i) X(2,i)*x1(1,i) x1(1,i);
        0 0 0 -X(1,i) -X(2,i) -1 X(1,i)*x1(2,i) X(2,i)*x1(2,i) x1(2,i);];
    M1 = [M1; M11];

    M22 = [-X(1,i) -X(2,i) -1 0 0 0 X(1,i)*x2(1,i) X(2,i)*x2(1,i) x2(1,i);
        0 0 0 -X(1,i) -X(2,i) -1 X(1,i)*x2(2,i) X(2,i)*x2(2,i) x2(2,i);];
    M2 = [M2; M22];
    
    M33 = [-X(1,i) -X(2,i) -1 0 0 0 X(1,i)*x3(1,i) X(2,i)*x3(1,i) x3(1,i);
        0 0 0 -X(1,i) -X(2,i) -1 X(1,i)*x3(2,i) X(2,i)*x3(2,i) x3(2,i);];
    M3 = [M3; M33];

end
M1
M2
M3

% U1,S1,V1 represent the left singular vectors, singular values, and right singular vectors respectively.
[U1,S1,V1] = svd(M1);
[U2,S2,V2] = svd(M2);
[U3,S3,V3] = svd(M3);

H1 = [reshape(V1(:,9),3,3)']
H2 = [reshape(V2(:,9),3,3)']
H3 = [reshape(V3(:,9),3,3)']


N1 = [H1(1,1)*H1(1,2),        H1(1,2)*H1(2,1) + H1(1,1)*H1(2,2),      H1(1,2)*H1(3,1) + H1(1,1)*H1(3,2),      H1(2,1)*H1(2,2),         H1(2,2)*H1(3,1) + H1(2,1)*H1(3,2),       H1(3,1)*H1(3,2);
      H1(1,1)^2 - H1(1,2)^2,  2*(H1(1,1)*H1(2,1) - H1(1,2)*H1(2,2)),  2*(H1(1,1)*H1(3,1) - H1(1,2)*H1(3,2)),  H1(2,1)^2 - H1(2,2)^2,   2*(H1(2,1)*H1(3,1) - H1(2,2)*H1(3,2)),   H1(3,1)^2 - H1(3,2)^2];
N2 = [H2(1,1)*H2(1,2),        H2(1,2)*H2(2,1) + H2(1,1)*H2(2,2),      H2(1,2)*H2(3,1) + H2(1,1)*H2(3,2),      H2(2,1)*H2(2,2),         H2(2,2)*H2(3,1) + H2(2,1)*H2(3,2),       H2(3,1)*H2(3,2);
      H2(1,1)^2 - H2(1,2)^2,  2*(H2(1,1)*H2(2,1) - H2(1,2)*H2(2,2)),  2*(H2(1,1)*H2(3,1) - H2(1,2)*H2(3,2)),  H2(2,1)^2 - H2(2,2)^2,   2*(H2(2,1)*H2(3,1) - H2(2,2)*H2(3,2)),   H2(3,1)^2 - H2(3,2)^2];
N3 = [H3(1,1)*H3(1,2),        H3(1,2)*H3(2,1) + H3(1,1)*H3(2,2),      H3(1,2)*H3(3,1) + H3(1,1)*H3(3,2),      H3(2,1)*H3(2,2),         H3(2,2)*H3(3,1) + H3(2,1)*H3(3,2),       H3(3,1)*H3(3,2);
      H3(1,1)^2 - H3(1,2)^2,  2*(H3(1,1)*H3(2,1) - H3(1,2)*H3(2,2)),  2*(H3(1,1)*H3(3,1) - H3(1,2)*H3(3,2)),  H3(2,1)^2 - H3(2,2)^2,   2*(H3(2,1)*H3(3,1) - H3(2,2)*H3(3,2)),   H3(3,1)^2 - H3(3,2)^2];

N = [N1;N2;N3]

[b1, B1, A1, K1] = get_K(N)

K1

K_inv = inv(K1)

Rt1 = get_Rt(K_inv, H1)
Rt2 = get_Rt(K_inv, H2)
Rt3 = get_Rt(K_inv, H3)

error = 0;
for i=1:size(X, 2)
   pixel = x1(:, i);
   pixel_est = K1*Rt1*X(:, i);
   pixel_est = pixel_est / pixel_est(3);
   
   error = error + norm((pixel - pixel_est), 2);
end
error
%% functions

function [b, B, A, K] = get_K(N)
    [U,S,V] = svd(N);
    b = V(:, end);
    
    B = [b(1) b(2) b(3); 
        b(2) b(4) b(5); b(3) b(5) b(6)];
    B = B / B(3, 3);
    
    A = chol(B, 'lower');
    K = (inv(A))';
    
    K = K / K(3, 3)
end

function Rt = get_Rt(K_inv,H)
    n = 1 / norm(K_inv   * [H(1,1) H(2,1) H(3,1)]', 2);

    r1 = K_inv  * [H(1,1) H(2,1) H(3,1)]';
    r2 = K_inv  * [H(1,2) H(2,2) H(3,2)]';
    r3 = cross(r1, r2);
    
    t = n * (K_inv  * [H(1,3) H(2,3) H(3,3)]');
    
    
    Rt = [r1, r2, r3, t];

end