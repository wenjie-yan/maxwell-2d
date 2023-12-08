clear
dx=0.01;
Lx=1;
x=-Lx:dx:Lx;
a=1;
A = (-2*eye(length(x))+diag(ones(1,length(x)-1),1)+diag(ones(1,length(x)-1),-1));
A(1,end)=1;
A(end,1)=1;
dt=0.0001;
t=0:dt:4;
u=zeros(length(x),length(t));
v=zeros(length(x),length(t));
%初始条件
Dx=0.1;
u0=exp(-x.^2/Dx^2);%高斯波包
v0=2*a*exp(-x.^2/Dx^2).*x/Dx^2;
u(:,1)=u0';
v(:,1)=v0';
for n=1:length(t)-1
    u(:,n+1)=u(:,n)+v(:,n)*dt;
    v(:,n+1)=v(:,n)+(a^2*A*u(:,n)/dx^2)*dt;
    if mod(n,100)==1
    plot(x,u(:,n+1))
    axis([x(1) x(end) 0 1.2])
    getframe
    end
end


