clear; clc; close all

%% ��������
c = 3e8; % ����,�������ȵ�λ��Ϊm��ʱ�䵥λ��Ϊs���Ա������ƥ��
eps0 = 8.85e-12; % ��糣��
mu0 = pi*4e-7; % ����ϵ��
le = 5; % z��Χ
tf = 1e-7;  % �������ʱ��
NDT = 10000;  % ʱ�䳤��
NDZ = 100; % ���������򳤶�
% ����FDTD��������dtС�ڵ���dz/v;
dt = tf/NDT;
dz = le/NDZ;
t = 0:dt:tf;
z = 0:dz:le;
f = 1e9;
wnige = 2*pi*f;  % Ƶ��
Z0 = 377;

%% FDTD���
E = zeros(NDT+1,NDZ+1);
H = zeros(NDT+1,NDZ);
for i = 1:NDT
    E(i+1,1) = cos(wnige*t(i+1));
    E(i+1,NDZ+1) = 0;
    H(i+1,1) = cos(wnige*t(i+1))/Z0;
    H(i+1,NDZ) = 0;
    for j = 1:NDZ-1
        E(i+1,j+1) = E(i,j+1) - dt/dz * 1/eps0 * (H(i,j+1) - H(i,j));
        H(i+1,j) = H(i,j) - dt/dz * 1/mu0 * (E(i+1,j+1) - E(i+1,j));
    end
    if mod(i,10) == 0
        figure(1)
        plot(z,E(i,:),'r-','Linewidth',2)
        figure(2)
        plot(z(2:end),H(i,:),'b-','Linewidth',2)
        pause(0.05)
    end
end