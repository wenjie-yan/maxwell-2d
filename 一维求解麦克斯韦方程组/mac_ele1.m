% 初始化仿真参数和初始条件
Nz = 200;   % 网格点数
T = 500;    % 总时间步数
frequency = 0.1; % 频率
z0 = 1;     % 特性阻抗
Dx = 10;    % Dx 是初始化条件的宽度参数

Ey = cos(2 * pi * frequency * (0:T) / T); % 初始条件 Ex(0, t) = cos(wt)
Hy = zeros(1, Nz);
Esrc = Ey(1); % 源的初始条件

% 边界条件
Es = 0;  % 边界条件 Es
E3 = 0;  % 边界条件 E3

% 初始化可视化
SHOW_PLOT = true;  % 是否显示图形
h3 = figure;
h1 = plot(1:Nz, Ey(1:Nz), 'b', 1:Nz, Hy, 'r');
xlabel('z direction');
ylabel('Ey and Hy Fields');
axis([1 Nz -3 3]);
legend({'Ey','Hy'},'Location','SouthEast');

% 迭代循环
for i = 1:T
    % 更新 Hy
    for nz = 1:(Nz - 1)
        if nz == 1
            Hy(nz) = Hy(nz) + (1 / z0) * cos(2 * pi * frequency * i / T) - 0.5 * (Ey(nz + 1) - Esrc);
        else
            Hy(nz) = Hy(nz) + 0.5 * (Ey(nz + 1) - Ey(nz));
        end
    end
    Hy(Nz) = Hy(Nz) + 0.5 * (0 - Ey(Nz));
    
    % 更新 Ex
    for nz = 1:Nz
        if nz == 1
            Ey(nz) = Ey(nz) + 0.5 * (Hy(nz) - 0) - 0.5 * Esrc;
        else
            Ey(nz) = Ey(nz) + 0.5 * (Hy(nz) - Hy(nz - 1));
        end
    end
    
    % 在每个时间步中更新可视化
    if SHOW_PLOT
        remain = mod(i, 10);
        if remain == 0
            set(h1(1), 'YData', Ey(1:Nz));
            set(h1(2), 'YData', Hy);
            title(['Step Num: ' num2str(i) '/' num2str(T)]);
            drawnow;
        end
    end
end