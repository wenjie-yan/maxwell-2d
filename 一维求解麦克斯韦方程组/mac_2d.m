%一维波动方程，三层显格式求解法
h=0.1;
tau=0.1*h;
r=tau/h;N=1/h;
M=2/tau;
x=0:h:1;
t=0:tau:2;
u=sin(pi*x);%计算 t= 时刻的u值
u(1,11)=0;
for j=2:N
u(2,j)=0.5*r^2*u(1,j+1)+(1-r^2)*u(1,j)+0.5*r^2*u(1,j-1);
end
%定义x=0边界上的数值
for k=1:M+1
    u(k,1)=cos(pi/3*k);
end
%定义x=1边界上的数值
for k=1:M+1
    u(k,N+1)=0;
end
%迭代计算开始,差分格式
for k=2:M
    for j=2:N
        u(k+1,j)=r^2*u(k,j+1)+2*(1-r^2)*u(k,j)+r^2*u(k,j-1)-u(k-1,j);
        

    end
    % if mod(k, 10) == 0
    %     clf;
    %     hold on;
    %     plot(u, 'LineWidth', 2);
    %     ylim([-1.5, 1.5]);
    %     xlim([0, 30]);
    %     title(['Step ', num2str(T)])
    %     drawnow;
        
    % end
end
u(201,:)=zeros(1,11);
%计算 k=201 行的数值解
u2(201,11)=0;
for j=2:N
    u2(201,j)=r^2*u(200,j+1)+2*(1-r^2)*u(200,j)+r^2*u(200,j-1)-u(199,j);
end
u=u+u2;
u=rot90(u,2);%将矩阵u 旋转 180度赋值于u
%作出图像
[x,t]=meshgrid(0:0.1:1,0:0.01:2);%划分网格%作出数值解的函数图像
subplot(2,2,1) ;
mesh(x, t,u) ;
title('u(x,t)数值解的函数图像');
xlabel('x变量');
ylabel('t 变量');
zlabel('u值');
%作出精确解的函数图像
subplot (2, 2, 2) ;
ul=cos(pi*t).*sin(pi*x);
mesh(x, t,ul);
title('u(x,t)精确解的函数图像');
xlabel('x变量');
ylabel('t 变量');
zlabel('u 值');
%作出 t=0.5，1.0，1.5，2.0 时刻的绝对误差图像
subplot(2,2, 3) ;
wucha=abs(u-ul);
x=0:h:1;
plot(x, wucha(51,:),'g*-');
hold on
grid on
hold on
grid on
plot(x, wucha(101,:),'ro-');
hold on
plot(x, wucha(151,:),'ks-');
hold on
plot(x, wucha(201,:),'mp-');title('t=0.5，1.0，1.5，2.0 时刻的绝对误差函数图像');xlabel('x变量');ylabel( '绝对误差值;');legend(' t=0.5','t=1.0','t=1.5','t=2.0');%作出 t=0.5，1.0，1.5，2.0 时刻的数值解函数图像subplot (2, 2, 4) ;
x=0:h:1;
plot(x,u(51,:),'g*-');
hold on
grid on
plot(x,u(101,:),'ro-');
hold on
plot(x,u(151,:),'ks-');
hold on
plot(x,u(201,:),'mp-');title('t=0.5，1.0，1.5，2.0 时刻的数值解函数图像');xlabel('x变量'); ylabel( 'u值');legend('t=0.5','t=l.0','t=1.5','t=2.0');%当然也可以作出u(x,t)绝对误差的函数图像%mesh(x, t, wucha);
%title(’u(x,t)绝对误差的函数图像’);
%xlabel(’x变量);
%ylabel(’t 变量');
%zlabel(绝对误差值’);