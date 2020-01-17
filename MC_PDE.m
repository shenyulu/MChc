function u = MC_PDE(temp_90,M)

  
    %layer:����
    %temp_90����Ƥ�����нӴ��Ĳ�Σ�x=L����ʱ��ı仯
    %M��Ͷ�����

%%%%%%%%%%%%%%��4����㲿��

    
%%%
%
% temp_90 ��90�������߽�����
% k���ȴ���ϵ��
% edge_first����ֵ����
% n1����ʼʱ��Ǳ�
% n2����ֹʱ��Ǳ�
% initial_value����ֵ����
% l_end������ά��ֹ��
%
%%%

%�ȴ�������ϵ��
k = 0.028;
%�ڲ�߽�������Ƥ����
edge_first = repmat(37,1,5401);
%����ά��ʼ��
edge_location = 1;
edge_location = int32(edge_location);
% ����ά��ֹ��
l_end = 50;
l_end = int32(l_end);

%% x����
% x���򻮷ּ��
h = 0.1;
% x����ֵ��x �� ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x�����������
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t����
% t���򻮷ּ��
tau =1;
% ʱ��ֵ��t �� ( 1, t_max)
t_max = 5401;
% ʱ���������
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%Ͷ�����
%M = 10;
temp_total = 0;

%ϵ�����
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);

%%��ʼ��
u = zeros(x_times+1,t_tims) * nan;
% ��ֵ����
u(:,1) = 37;

% ��ֵ����
u(edge_location,:) = edge_first;
u(x_times+1, :) = temp_90(:,2)';

 

for j_value = (edge_location + 1) : l_end
    disp(j_value)
    for n = 2:5401
        % Ͷ����
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
                    
                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

%%%%%%%%%%%%%%��4����㲿��end

%%%%%%%%%%%%%%%��3����㲿��

        
%
% temp_90��90�������߽�����
% k���ȴ���ϵ��
% edge_first����ֵ����
% n1����ʼʱ��Ǳ�
% n2����ֹʱ��Ǳ�
% initial_value����ֵ����
% l_end������ά��ֹ��
%
%%%

% ����90�������߽�����
% temp_90 = csvread('bianjietj.csv');
%�ȴ�������ϵ��
k = 0.045;
%����ά��ʼ��
edge_location = 50;
% ����ά��ֹ��
l_end = 86;
l_end = int32(l_end);

%% x����
% x���򻮷ּ��
h = 0.1;
% x����ֵ��x �� ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x�����������
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t����
% t���򻮷ּ��
tau =1;
% ʱ��ֵ��t �� ( 0, t_max)
t_max = 5401;
% ʱ���������
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%Ͷ�����
%M = 10;
temp_total = 0;

%ϵ�����
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);

for j_value = (edge_location+1) : l_end
    disp(j_value)
    for n = 2:5401
        % Ͷ����
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
                    
                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

%%%%%%%%%%%%%%%��3����㲿��end

%%%%%%%%%%%%%%%%��2����㲿��

        

%
% temp_90��90�������߽�����
% k���ȴ���ϵ��
% edge_first����ֵ����
% n1����ʼʱ��Ǳ�
% n2����ֹʱ��Ǳ�
% initial_value����ֵ����
% l_end������ά��ֹ��
%
%%%


%�ȴ�������ϵ��
k = 0.37;
%�ڲ�߽�������Ƥ����
% edge_first = repmat(37,1,5401);
%����ά��ʼ��
edge_location = 86;
edge_location = int32(edge_location);
% ����ά��ֹ��
l_end = 146;
l_end = int32(l_end);

%% x����
% x���򻮷ּ��
h = 0.1;
% x����ֵ��x �� ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x�����������
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t����
% t���򻮷ּ��
tau =1;
% ʱ��ֵ��t �� ( 0, t_max)
t_max = 5401;
% ʱ���������
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%Ͷ�����
%M = 10;
temp_total = 0;

%ϵ�����
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);


for j_value = (edge_location+1) : l_end
    disp(j_value)
    for n = 2:5401
        % Ͷ����
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
                    
                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

%%%%%%%%%%%%%%%%��2����㲿��end

%%%%%%%%%%%%%%%%%��1����㲿��

        

%
% temp_90��90�������߽�����
% k���ȴ���ϵ��
% edge_first����ֵ����
% n1����ʼʱ��Ǳ�
% n2����ֹʱ��Ǳ�
% initial_value����ֵ����
% l_end������ά��ֹ��
%
%%%

% ����90�������߽�����
% temp_90 = csvread('bianjietj.csv');
%�ȴ�������ϵ��
k = 0.082;
%�ڲ�߽�������Ƥ����
% edge_first = repmat(37,1,5401);
%����ά��ʼ��
edge_location = 146;
edge_location = int32(edge_location);
% ����ά��ֹ��
l_end = 152;
l_end = int32(l_end);

%% x����
% x���򻮷ּ��
h = 0.1;
% x����ֵ��x �� ( 0, x_max)
x_max = 15.2;
%x_initial = 0 : h :1;
% x�����������
x_times = x_max / h - 1;
x_times = int32(x_times);

%% t����
% t���򻮷ּ��
tau =1;
% ʱ��ֵ��t �� ( 0, t_max)
t_max = 5401;
% ʱ���������
t_tims = t_max / tau;
t_tims = int32(t_tims);

%%Ͷ�����
%M = 10;
temp_total = 0;

%ϵ�����
r = k * ( tau / (h^2) );
coef = 1 / (1 + r);

%ע��l_lend�߽��������
for j_value = (edge_location+1) : (l_end - 1)
    disp(j_value)
    for n = 2:5401
        % Ͷ����
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
                    
                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

                    %���麯��
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

%%%%%%%%%%%%%%%%%��1����㲿��
    
end