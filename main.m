% 导入90分钟外界边界条件
temp_90 = csvread('bianjietj.csv');
tic
u1 = MC_PDE_mex(temp_90,100);
toc

x =1:152;
t = 1:5401;
[X,T] = meshgrid(x,t);
surf(X,T,u1')

xlabel('Thickness')
ylabel('Time')
zlabel('Temperature')
shading flat
colorbar