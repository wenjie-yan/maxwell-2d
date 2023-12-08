clc;%清理命令行窗口
clear all;%清理工作区
%读取excel数据，同目录下
tx=xlsread('B1.xlsx');

Angle_1=tx(1,:);%第1行数据为角度
Value_1=tx(2,:);%第2行数据为极化值
Angle_2=tx(3,:);%第3行数据为角度
Value_2=tx(4,:);%第4行数据为极化值

% 创建极坐标图，不使用 deg2rad
polarplot(Angle_1 * (pi/180), Value_1, 'r', 'DisplayName', '转动上方罗盘');
hold on;
polarplot(Angle_2 * (pi/180), Value_2, 'b', 'DisplayName', '转动下方罗盘');
legend;



% 可以使用标题来添加标题
title('极坐标图');

% 可以使用下面的语句来调整标题的位置
% title('极坐标图', 'Position', [0, 1.4]);