clear; clc; close all

%% 参数设置
c = 3e8; % 光速,后续长度单位都为m、时间单位都为s，以便与光速匹配
eps0 = 8.85e-12; % 介电常数
mu0 = pi*4e-7; % 导磁系数
le = 5; % z范围
tf = 1e-7;  % 最终求解时间
NDT = 10000;  % 时间长度
NDZ = 100; % 波传播方向长度
% 满足FDTD条件，即dt小于等于dz/v;
dt = tf/NDT;
dz = le/NDZ;
t = 0:dt:tf;
z = 0:dz:le;
f = 1e9;
wnige = 2*pi*f;  % 频率
Z0 = 377;

%% FDTD求解
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