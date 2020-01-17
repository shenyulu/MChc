function u = MC_PDE(temp_90,M)

  
    %layer:层数
    %temp_90：与皮肤密切接触的层次（x=L）随时间的变化
    %M：投点次数

%%%%%%%%%%%%%%第4层计算部分

    
%%%
%
% temp_90 ：90分钟外界边界条件
% k：热传导系数
% edge_first：边值条件
% n1：起始时间角标
% n2：终止时间角标
% initial_value：初值条件
% l_end：长度维终止点
%
%%%

%热传导方程系数
k = 0.028;
%内层边界条件（皮肤）
edge_first = repmat(37,1,5401);
%长度维起始行
edge_location = 1;
edge_location = int32(edge_location);
% 长度维终止点
l_end = 50;
l_end = int32(l_end);

%% x方向
% x方向划分间隔
h = 0.1;
% x方向：值域，x ∈ ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x方向迭代次数
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t方向
% t方向划分间隔
tau =1;
% 时间值域：t ∈ ( 1, t_max)
t_max = 5401;
% 时间迭代次数
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%投点次数
%M = 10;
temp_total = 0;

%系数求解
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);

%%初始化
u = zeros(x_times+1,t_tims) * nan;
% 初值条件
u(:,1) = 37;

% 边值条件
u(edge_location,:) = edge_first;
u(x_times+1, :) = temp_90(:,2)';

 

for j_value = (edge_location + 1) : l_end
    disp(j_value)
    for n = 2:5401
        % 投掷点
        temp = zeros(1,M);
        for point_num = 1 : M
            t = int32(n);
            l = int32(j_value);
            while(1)
                xi = rand();
                %count_num = 1;
                % u(j_value + 1, t)
                if xi > 0 && xi < coef * (1/2 * r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    
                    %检验函数
                    if l == 152
                        basicvalue = u(152,t);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    elseif   l == 1
                        basicvalue = u(1,t);
                        temp(point_num) = basicvalue;
                        signal = 1; 
                    elseif t == 1
                        basicvalue = u(l,1);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end
                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);

                    if signal == 1
                        break
                    end

                elseif xi >= coef * (1/2 * r) && xi < coef * ( r)
                    %count_num = count_num * 1/6;
                    l = l - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value - 1, t);
                elseif xi >= coef * ( r) && xi < coef
                    %count_num = count_num * 1/3;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                elseif xi >= coef && xi < coef * (1 + (1/2)*r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value + 1, t-1);
                elseif xi >= coef * (1 + (1/2)*r) && xi < 1
                    %count_num = count_num * 1/6;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                end
            end    
        end
        temp_total = sum(temp);
        %u(j_value, n) = coef * temp_total / M;
        u(j_value, n) =  temp_total / M;
         %temp = zeros(1,M);
         %temp_total = 0;
    end
end

%%%%%%%%%%%%%%第4层计算部分end

%%%%%%%%%%%%%%%第3层计算部分

        
%
% temp_90：90分钟外界边界条件
% k：热传导系数
% edge_first：边值条件
% n1：起始时间角标
% n2：终止时间角标
% initial_value：初值条件
% l_end：长度维终止点
%
%%%

% 导入90分钟外界边界条件
% temp_90 = csvread('bianjietj.csv');
%热传导方程系数
k = 0.045;
%长度维起始行
edge_location = 50;
% 长度维终止点
l_end = 86;
l_end = int32(l_end);

%% x方向
% x方向划分间隔
h = 0.1;
% x方向：值域，x ∈ ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x方向迭代次数
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t方向
% t方向划分间隔
tau =1;
% 时间值域：t ∈ ( 0, t_max)
t_max = 5401;
% 时间迭代次数
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%投点次数
%M = 10;
temp_total = 0;

%系数求解
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);

for j_value = (edge_location+1) : l_end
    disp(j_value)
    for n = 2:5401
        % 投掷点
        temp = zeros(1,M);
        for point_num = 1 : M
            t = int32(n);
            l = int32(j_value);
            while(1)
                xi = rand();
                %count_num = 1;
                % u(j_value + 1, t)
                if xi > 0 && xi < coef * (1/2 * r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    
                    %检验函数
                    if l == 152
                        basicvalue = u(152,t);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    elseif   l == 1
                        basicvalue = u(1,t);
                        temp(point_num) = basicvalue;
                        signal = 1; 
                    elseif t == 1
                        basicvalue = u(l,1);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end
                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);

                    if signal == 1
                        break
                    end

                elseif xi >= coef * (1/2 * r) && xi < coef * ( r)
                    %count_num = count_num * 1/6;
                    l = l - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value - 1, t);
                elseif xi >= coef * ( r) && xi < coef
                    %count_num = count_num * 1/3;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                elseif xi >= coef && xi < coef * (1 + (1/2)*r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value + 1, t-1);
                elseif xi >= coef * (1 + (1/2)*r) && xi < 1
                    %count_num = count_num * 1/6;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                end
            end    
        end
        temp_total = sum(temp);
        %u(j_value, n) = coef * temp_total / M;
        u(j_value, n) =  temp_total / M;
         %temp = zeros(1,M);
         %temp_total = 0;
    end
end

%%%%%%%%%%%%%%%第3层计算部分end

%%%%%%%%%%%%%%%%第2层计算部分

        

%
% temp_90：90分钟外界边界条件
% k：热传导系数
% edge_first：边值条件
% n1：起始时间角标
% n2：终止时间角标
% initial_value：初值条件
% l_end：长度维终止点
%
%%%


%热传导方程系数
k = 0.37;
%内层边界条件（皮肤）
% edge_first = repmat(37,1,5401);
%长度维起始行
edge_location = 86;
edge_location = int32(edge_location);
% 长度维终止点
l_end = 146;
l_end = int32(l_end);

%% x方向
% x方向划分间隔
h = 0.1;
% x方向：值域，x ∈ ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x方向迭代次数
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t方向
% t方向划分间隔
tau =1;
% 时间值域：t ∈ ( 0, t_max)
t_max = 5401;
% 时间迭代次数
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%投点次数
%M = 10;
temp_total = 0;

%系数求解
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);


for j_value = (edge_location+1) : l_end
    disp(j_value)
    for n = 2:5401
        % 投掷点
        temp = zeros(1,M);
        for point_num = 1 : M
            t = int32(n);
            l = int32(j_value);
            while(1)
                xi = rand();
                %count_num = 1;
                % u(j_value + 1, t)
                if xi > 0 && xi < coef * (1/2 * r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    
                    %检验函数
                    if l == 152
                        basicvalue = u(152,t);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    elseif   l == 1
                        basicvalue = u(1,t);
                        temp(point_num) = basicvalue;
                        signal = 1; 
                    elseif t == 1
                        basicvalue = u(l,1);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end
                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);

                    if signal == 1
                        break
                    end

                elseif xi >= coef * (1/2 * r) && xi < coef * ( r)
                    %count_num = count_num * 1/6;
                    l = l - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value - 1, t);
                elseif xi >= coef * ( r) && xi < coef
                    %count_num = count_num * 1/3;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                elseif xi >= coef && xi < coef * (1 + (1/2)*r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value + 1, t-1);
                elseif xi >= coef * (1 + (1/2)*r) && xi < 1
                    %count_num = count_num * 1/6;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                end
            end    
        end
        temp_total = sum(temp);
        %u(j_value, n) = coef * temp_total / M;
        u(j_value, n) =  temp_total / M;
         %temp = zeros(1,M);
         %temp_total = 0;
    end
end

%%%%%%%%%%%%%%%%第2层计算部分end

%%%%%%%%%%%%%%%%%第1层计算部分

        

%
% temp_90：90分钟外界边界条件
% k：热传导系数
% edge_first：边值条件
% n1：起始时间角标
% n2：终止时间角标
% initial_value：初值条件
% l_end：长度维终止点
%
%%%

% 导入90分钟外界边界条件
% temp_90 = csvread('bianjietj.csv');
%热传导方程系数
k = 0.082;
%内层边界条件（皮肤）
% edge_first = repmat(37,1,5401);
%长度维起始行
edge_location = 146;
edge_location = int32(edge_location);
% 长度维终止点
l_end = 152;
l_end = int32(l_end);

%% x方向
% x方向划分间隔
h = 0.1;
% x方向：值域，x ∈ ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x方向迭代次数
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t方向
% t方向划分间隔
tau =1;
% 时间值域：t ∈ ( 0, t_max)
t_max = 5401;
% 时间迭代次数
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%投点次数
%M = 10;
temp_total = 0;

%系数求解
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);

%注意l_lend边界无需遍历
for j_value = (edge_location+1) : (l_end - 1)
    disp(j_value)
    for n = 2:5401
        % 投掷点
        temp = zeros(1,M);
        for point_num = 1 : M
            t = int32(n);
            l = int32(j_value);
            while(1)
                xi = rand();
                %count_num = 1;
                % u(j_value + 1, t)
                if xi > 0 && xi < coef * (1/2 * r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    
                    %检验函数
                    if l == 152
                        basicvalue = u(152,t);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    elseif   l == 1
                        basicvalue = u(1,t);
                        temp(point_num) = basicvalue;
                        signal = 1; 
                    elseif t == 1
                        basicvalue = u(l,1);
                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end
                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);

                    if signal == 1
                        break
                    end

                elseif xi >= coef * (1/2 * r) && xi < coef * ( r)
                    %count_num = count_num * 1/6;
                    l = l - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value - 1, t);
                elseif xi >= coef * ( r) && xi < coef
                    %count_num = count_num * 1/3;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                elseif xi >= coef && xi < coef * (1 + (1/2)*r)
                    %count_num = count_num * 1/6;
                    l = l + 1;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value + 1, t-1);
                elseif xi >= coef * (1 + (1/2)*r) && xi < 1
                    %count_num = count_num * 1/6;
                    t = t - 1;

                    %检验函数
                    if l == 152 || l == 1 || t == 1 
                        basicvalue = u(l,t);

                        temp(point_num) = basicvalue;
                        signal = 1;
                    else
                        temp(point_num) = nan;
                        signal = 0;
                    end

                    %[temp(point_num),signal] = edge_check(u,l,t,count_num);
                    if signal == 1
                        break
                    end
                    %temp(point_num) = u(j_value, t-1);
                end
            end    
        end
        temp_total = sum(temp);
        %u(j_value, n) = coef * temp_total / M;
        u(j_value, n) =  temp_total / M;
         %temp = zeros(1,M);
         %temp_total = 0;
    end
end

%%%%%%%%%%%%%%%%%第1层计算部分
    
end