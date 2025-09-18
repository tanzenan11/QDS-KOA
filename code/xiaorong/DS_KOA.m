function [Foutput,time] = QDS_KOA(P,DA,CA,HA,Td,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,time,C)
% 四目标：站的数量，平滑程度，利润，能耗
% 规定参数
tic %开始计时
np = 100;  %行星数
Tmax = 2; %最大迭代次数
dim = 3; %目标数
[n,A] = size(P); % n代表零件数
favg = 0;   %存储平均适应度值
Q = zeros(20,10);
pareto_frontier = [10 -3000 40];  %计算HV值真实点
ub = 1;
lb = 0;
IGD_one_sum = 0;    %第一代个体的IGD值求和
IGD_one_max = 0;    %第一代个体的IGD值最大值
IGD_one_S2 = 0;
Gamma = 0.6;%Q-learning远见系数
Alpha = 0.4;%Q-learning学习系数
%Q-learning 领域搜索比例-固定比例
Proportion=[
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25
    0.25, 0.25, 0.25, 0.25];
Sun_Pos = zeros(1,n*2); % 一个矢量，包含迄今为止最好的解决方案，代表太阳
Sun_Score = inf(1,dim); % 一个标量变量，包含迄今为止最好的分数

%---------------------控制参数------------------------%
%
Tc=3;   %loop control constant
M0=0.1;  %the initial value
lambda=15; %the constant
% 轨道偏心率(e)
orbital=rand(1,np); %% Eq.(4)   生成与种群数个数相同的偏心率（0-1）
%%轨道周期(T)
T=abs(randn(1,np)); %% Eq.(5)   根据正态分布生成与种群数个数相同的偏心率
% 种群初始化
EC=[];  % 帕累托解
[Chrom,A1] = Shengcheng(np,n,P,DA,HA);
positions=rand(1,n); % 初始化行星的位置
% 对Positions排序
[A,B] = sort(positions);
for j = 1:np
    Chromone = Chrom(j,1:n);
    for k = 1:n
        Positions(j,Chromone(k)) = A(1,k);
    end
end

t = 0;   %函数求值计数器
for j = 1:np
    % %计算目标值
    [Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值
    if calculate_igd(pareto_frontier,fit(j,:)) <= calculate_igd(pareto_frontier,Sun_Score)
        Sun_Score=fit(j,:); % 更新到目前为止最好的分数
        Sun_Pos=Positions(j,:); % 更新目前为止最好的解决方案
        Sun_Pos(1,n+1) = fit(j,1);
        Sun_Pos(1,n+2) = fit(j,2);
        Sun_Pos(1,n+3) = fit(j,3);
    end
end
%把第一代个体的IGD值求和
for j = 1:np
    IGD_one_sum = IGD_one_sum + calculate_igd(pareto_frontier,fit(j,:));
    if calculate_igd(pareto_frontier,fit(j,:)) > IGD_one_max
        IGD_one_max = calculate_igd(pareto_frontier,fit(j,:));
    end
end

% 计算第一代个体的种群多样性
for j = 1:np
    IGD_one_S2 = IGD_one_S2 + abs(calculate_igd(pareto_frontier,fit(j,:)) - IGD_one_sum/np);
end

while t < Tmax
    [Order] = sort(PL_Fit);
    % 函数求值t时的最差适应度值
    worstFitness = Order(np); %% Eq.(11)
    M = M0*(exp(-lambda*(t/Tmax))); %% Eq. (12) %M是一个随时间（t）呈指数下降的函数，用于控制搜索精度
    % 计算机R表示目前最佳解和第i个解之间的欧氏距离
    for j = 1:np
        R(j) = 0;
        for k = 1:dim
            R(j)=R(j)+(Sun_Pos(n+k)-fit(j,k))^2; %% Eq.(7)    %欧式距离应用到拆卸线得想一下修改
        end
        R(j)=sqrt(R(j));
    end
    for j = 1:np
        Sum = 0;
        for k = 1:np
            Sum=Sum+(PL_Fit(k)-worstFitness);
        end
        MS(j) = rand*(sum(Sun_Score)/4-worstFitness)/(Sum); % Eq.(8)
        m(j)=(PL_Fit(j)-worstFitness)/(Sum); %% Eq.(9)
    end
    %Step 2: 定义引力(F)
    %根据万有引力定律计算太阳和第i颗行星之间的引力
    for j = 1:np
        Rnorm(j)=(R(j)-min(R))/(max(R)-min(R)); %% The normalized R (Eq.(24))XS和Xi之间的欧几里德距离归一化
        MSnorm(j)=(MS(j)-min(MS))/(max(MS)-min(MS)); %% The normalized MS MS质量归一化
        Mnorm(j)=(m(j)-min(m))/(max(m)-min(m)); %% The normalized m      m质量归一化
        Fg(j)=orbital(j)*M*((MSnorm(j)*Mnorm(j))/(Rnorm(j)*Rnorm(j)+eps))+(rand); %% Eq.(6)  计算引力
    end
    % A1表示物体I在时刻t的椭圆轨道的半长轴，
    for j = 1:np
        a1(j)=rand*(T(j)^2*(M*(MS(j)+m(j))/(4*pi*pi)))^(1/3); %% Eq.(23)
    end
    for j =1:np
        % a2是一个从-1到-2逐渐减小的循环控制参数
        a2=-1+-1*(rem(t,Tmax/Tc)/(Tmax/Tc)); %% Eq.(29)
        % 而n1是从1到−2的线性递减因子
        n1=(a2-1)*rand+1; %% Eq.(28)
        a=randi(np);  % 取出一个随机解
        b=randi(np); % 取出一个随机解
        rd=rand(1,n); %根据正态分布生成的向量 r5
        r=rand; %r是0到1的一个随机数 r4
        %  一个随机分配的二进制向量
        U1=rd<r; %% Eq.(21)    %好像反了需要最后效果求证
        O_P=Chrom(j,:); %% 存储第i个解
        O_Positions = Positions(j,:);%% 存储原本位置
        % Step 6: 更新行星
        if rand<rand
            % h是用于控制在时间t时太阳与当前行星之间的距离的自适应因子
            h=(1/(exp(n.*randn))); %% Eq.(27)
            %基于三个解决方案的平均向量:当前解决方案、迄今最佳解决方案和随机选择的解决方案
            Xm=(Positions(a,:)+Sun_Pos(1,1:n)+Positions(j,:))/3.0;
            Positions(j,:)=Positions(j,:).*U1+(Xm+h.*(Xm-Positions(b,:))).*(1-U1); %% Eq.(26)
        else
            % Step 3: 计算物体速度
            % 与当前行星的搜索方向相反或偏离的标志
            if rand<0.5 %% Eq.(18)
                f=1;
            else
                f=-1;
            end
            L=(M*(MS(j)+m(j))*abs((2/(R(j)+eps))-(1/(a1(j)+eps))))^(0.5); %% Eq.(15)  eps为一个极小值数
            U=rd>rand(1,n); %% 二值向量Eq.(17)
            if Rnorm(j)<0.5 %% Eq.(13)
                M=(rand.*(1-r)+r); %% Eq.(16)
                l=L*M*U; %% Eq.(14)
                Mv=(rand*(1-rd)+rd); %% Eq.(20)
                l1=L.*Mv.*(1-U);%% Eq.(19)
                V(j,:)=l.*(2*rand*Positions(j,:)-Positions(a,:))+l1.*(Positions(b,:)-Positions(a,:))+(1-Rnorm(j))*f*U1.*rand(1,n).*(ub-lb); %% Eq.(13a)
            else
                U2=rand>rand; %% Eq. (22)
                V(j,:)=rand.*L.*(Positions(a,:)-Positions(j,:))+(1-Rnorm(j))*f*U2*rand(1,n).*(rand*ub-lb);  %% Eq.(13b)
            end %% End IF

            % Step 4: 逃离局部最优
            % 更新标志f到相反的方向或离开当前行星的搜索方向
            if rand<0.5 %% Eq.(18)
                f=1;
            else
                f=-1;
            end
            % Step 5 跟新位置
            Positions(j,:)=((Positions(j,:)+V(j,:).*f)+(Fg(j)+abs(randn))*U.*(Sun_Pos(1,1:n)-Positions(j,:))); %% Eq.(25)
        end
        %根据更新后的位置生成新种群,并检查是否合理,并且重新生成是否拆卸反序列
        Positionsone = Positions(j,:);
        [A,Chromone] = sort(Positionsone);
        % 对更新个体合理性进行判断
        for k = 1:n
            a = Chromone(1,k);
            A1 = find(P(:, a) == -1);  %前置或集合
            A2 = find(P(:, a) == 1);   %前置与集合
            if ~isempty(A1)
                if ~any(ismember(Chromone(1,1:k-1),A1))   % 或集合不满足
                    Chromone(1,:) = O_P(1,1:n);              % 以优秀父代个体作为子代
                    Positions(j,:) = O_Positions;
                    break;
                end
            end
            if ~isempty(A2)
                if ~all(ismember(A2, Chromone(1,1:k-1)))  % 与集合不满足
                    Chromone(1,:) = O_P(1,1:n);              % 以优秀父代个体作为子代
                    Positions(j,:) = O_Positions;
                    break;
                end
            end
        end
        %存放个体
        Chrom(j,1:n) = Chromone(1,:);
        %确定是否需要重新选择是否拆卸
        if ~isequal(Chromone, O_P(1:numel(Chromone)))
            Chrom(j,n+1:n*2) = 0;
            %重新选择
            P0 = P;
            indices_ca = find(DA == 1);  % 需求操作
            indices_ha = find(HA == 1);  % 危险操作
            indices = [indices_ca,indices_ha];    %所有有需求和危险的操作
            len = length(indices);
            A = indices;
            dis = [];    %拆卸集合
            while ~isempty(A2)
                len = length(indices);
                for h = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                    temp = indices(h);   %操作
                    A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
                    A1 = A1';
                    A = [A,A1];
                end
                A = unique(A);   %必须拆卸的步骤
                dis = [dis,A];
                %获得前置步骤的前置步骤
                A2 = A(~ismember(A, indices));
                indices = A2;
                A = indices;
            end
            dis = unique(dis, 'rows', 'stable');   %必须拆卸的步骤
            notdis = setdiff(1:n, dis);               %其余步骤
            len = length(notdis);
            if len > 0
                random_number = randi([1, len]); % 生成 0 到列表长度的随机数
                indices = randperm(len, random_number); % 生成 random_number 个不重复的随机索引
                random_dis = notdis(indices); % 随机取出的拆卸操作
                B = random_dis;               %重新存储一次随机取出的拆卸操作
                random_dis_fm = random_dis;
                A2 = [0];    %保证进入循环
                while ~isempty(A2)
                    len = length(random_dis);
                    for h = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                        temp = random_dis(h);   %操作
                        A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
                        A1 = A1';
                        A = [A,A1];
                    end
                    A = unique(A);   %必须拆卸的步骤
                    dis = [dis,A];
                    %获得前置步骤的前置步骤
                    A2 = A(~ismember(A, random_dis));
                    random_dis = A2;
                    A = random_dis;
                end
                random_dis_fm = unique(random_dis_fm, 'rows', 'stable');   %必须拆卸的步骤
                [~, idx1] = ismember(dis, Chromone);     %必须拆卸步骤
                Chrom(j,n+idx1) = 1;
                [~, idx2] = ismember(random_dis_fm, Chromone);  %随机的拆卸步骤
                Chrom(j,n+idx2) = 1;
            else
                [~, idx1] = ismember(dis, Chromone);     %必须拆卸步骤
                Chrom(j,n+idx1) = 1;
            end
        else
            Chrom(j,n+1:n*2) = O_P(1,n+1:n*2);
        end
        %计算函数值
        [Chrom_nouse,fit_one,PL_Fit1] = Cal(Chrom(j,:),n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值;
        %Step 7: 精英保留策略, Eq.(30)
        if calculate_igd(pareto_frontier,fit_one) <= calculate_igd(pareto_frontier,fit(j,:)) % 保留优解
            fit(j,:)=fit_one; % 保留优解
            % 更新全局最优解
            if calculate_igd(pareto_frontier,fit(j,:))<calculate_igd(pareto_frontier,Sun_Score) % 目前解比全局最优解好
                Sun_Score=fit(j,:); % 更新全局最优解适应度值
                Sun_Pos=Positions(j,:); % 更新全局最优解方案
                Sun_Pos(1,n+1) = fit(j,1);
                Sun_Pos(1,n+2) = fit(j,2);
                Sun_Pos(1,n+3) = fit(j,3);
            end
        else
            Chrom(j,:)=O_P;
        end %% End IF
    end

    % %计算目标值
    [Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值
    %计算这一代的IGD和以及IGD最大值
    IGD_t_sum = 0;
    IGD_t_max = 0;
    %把第t代个体的IGD值求和
    for j = 1:np
        IGD_t_sum = IGD_t_sum + calculate_igd(pareto_frontier,fit(j,:));
        if calculate_igd(pareto_frontier,fit(j,:)) > IGD_t_max
            IGD_t_max = calculate_igd(pareto_frontier,fit(j,:));   %求出第t代个体的IGD最大值
        end
    end
    s1 = IGD_t_sum/IGD_one_sum;
    % 将种群多样性进行归一化
    IGD_t_S2 = 0;
    for j = 1:np
        IGD_t_S2 = IGD_t_S2 + abs(calculate_igd(pareto_frontier,fit(j,:)) - IGD_t_sum/np);
    end
    s2 = IGD_t_S2/IGD_one_S2;     %种群多样性进行归一化结果
    %将最佳适应度归一化
    s3 = IGD_t_max / IGD_one_max;
    s = 0.35*s1 + 0.35*s2 + 0.3*s3;
    % Q-learning领域搜索
    [new_Chrom,Q] = Q_learning(Chrom,s,n,P,Q,Gamma,Alpha,t,Tmax,Proportion,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,pareto_frontier);
    Chrom = Chrom(:,1:n*2);
    Chrom = [new_Chrom;Chrom];
    [Chrom,fit,PL_Fit1] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
    Chrom = [fit,Chrom];
    frontlabel=quickrank(Chrom,fit);
    %去除重复值
    %frontlabel=xipai(frontlabel,5);
    if size(frontlabel,1)>2*np
        frontlabel=frontlabel(1:2*np,:);
    else
        while ~(size(frontlabel,1)==2*np)  %费时3秒
            if size(frontlabel,1)<2*np
                [buchongpop,A1]=Shengcheng(2*np-size(frontlabel,1),n,P,DA,HA);
                [~,fitvalue,~]=Cal(buchongpop,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
                buchongpop=[fitvalue,buchongpop];
                lastpopution=[frontlabel(:,2:end);buchongpop];  %合成种群
                frontvalue=lastpopution(:,1:3);  %%%三个目标值
                frontlabel=quickrank(lastpopution,frontvalue);
                %frontlabel=xipai(frontlabel,3);
            end
        end
    end
    on_inferior=feizhipei(frontlabel(:,2:end),frontlabel(:,2:4));
    min(on_inferior(:,1));
    if (length(find(frontlabel(:,1)==1)))<np  %----------------单目标极值并没有被保存下来
        fnum=0;                                                             %当前前沿面
        while numel(frontlabel(:,1),frontlabel(:,1)<=fnum+1)<=np         %判断前多少个面的个体能完全放入下一代种群
            fnum=fnum+1;
        end   %至少需要前fnum层的个体

        newnum=numel(frontlabel(:,1),frontlabel(:,1)<=fnum);             %前fnum个面的个体数
        Chrom=zeros(np,n*2);
        Chrom(1:newnum,:)=frontlabel(1:newnum,5:end);                   %将前fnum个面的个体复制入下一代 -------
        a=find(frontlabel(:,1)==(fnum+1));
        a_population=frontlabel(a,:);
        a_population1=crowded(a_population(:,2:end),np-newnum); %------对第fnum+1的个体进行排序，选择个体
        Chrom(newnum+1:np,:)=a_population1(:,1:end); %%已经筛选出了popsize个个体 。
    else
        Chrom=crowded(frontlabel(:,2:end),np);
    end
    t=t+1; %% 迭代次数加1
    if t>Tmax %% 检查终止条件
        break;
    end %% End IF
end
time = toc;
save('Q.mat', 'Q');
[Chrom,fit,PL_Fit1] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);
frontvalue = Newranking(fit);
Chrom(:,2*n+1) = frontvalue;
Chrom(:,2*n+2) = fit(:,1);
Chrom(:,2*n+3) = fit(:,2);
Chrom(:,2*n+4) = fit(:,3);
Chrom=select(np,frontvalue,Chrom,fit);

%%%帕累托解
L=1;
para=[];  % 帕累托解
for i=1:np
    if  Chrom(i,2*n+1)==1
        para(L,:)=Chrom(i,:);
        L=L+1;
    end
end
% 去除重复个体
[unique_para, idx, ~] = unique(para, 'rows', 'stable');
% disp('para解集');
% disp(unique_para);
para = unique_para(:,n*2+2:n*2+4);
Foutput = para;
% 存储运行时间到变量 time

end

% 按照不同目标分离数据
% % 假设数据保存在名为 para 的变量中
% smoothness = para(:, 2);
% profit = para(:, 3);
% energy_consumption = para(:, 4);
% total_sites = para(:, 1);
%
% % 使用总站数创建一个颜色向量，颜色可以根据总站数的不同值来变化
% colors = total_sites; % 假设使用总站数作为颜色值
%
% % 绘制三维散点图，并根据颜色向量设置颜色
% scatter3(smoothness, profit, energy_consumption, [], colors, 'filled');
% xlabel('平滑度');
% ylabel('利润');
% zlabel('能耗');
% title('三维散点图表示平滑度、利润和能耗关系（根据总站数着色）');
%
% % 添加颜色条
% colorbar; % 显示颜色条，显示不同颜色对应的总站数范围


%%删除重复值
function  [frontlabel]=xipai(frontlabel,cengnum)
fl=[];
lab=[];
labelpop=[];
a=[];
for k=1:cengnum
    weizhi=find(frontlabel(:,1)==k);
    a1=find(frontlabel(:,1)==k);
    a=[a;a1];
    labelpop=frontlabel(weizhi,:);
    if length(weizhi)==1
        fl=[fl; labelpop];
        labelpop=[];
        continue
    end
    for k1=1:(length(weizhi)-1)
        for k2=(k1+1):length(weizhi)
            chazhi1=(labelpop(k1,2:4)==labelpop(k2,2:4));
            if sum(chazhi1)==(length(chazhi1))
                lab=[lab,k2];
            end
        end
        labelpop(lab,:)=[];
        weizhi=find( labelpop(:,1)==k);
        lab=[];
    end
    fl=[fl; labelpop];
    labelpop=[];
end
frontlabel=[fl;frontlabel((length(a)+1):end,:)];
end

%% 生成种群
function [Chrom,A1] = Shengcheng(np,n,P,DA,HA)
Chrom = zeros(np,n*2);
% 生成操作顺序
for i = 1:np
    P0 = P;
    bused = [];   %存放已选择的任务
    for j = 1:n
        A1 = [];
        A2 = [];
        % 生成or集合
        a = 1;
        for k = 1:n   %遍历行
            if ismember(-1,P0(k,:)) == 1 %存在or约束,记录行数（前驱零件）
                A1(1,a) = k;
                a = a + 1;
            end
        end
        for k = 1:n   %遍历列
            if ismember(-1,P0(:,k)) == 1 %存在or约束,记录列数（前驱零件）
                A1(1,a) = k;
                a = a + 1;
            end
        end
        % 遍历列生成优先任务集合
        c = 1;
        for k = 1:n
            if any(P0(:,k)) == 0     %k零件没有任何前驱约束
                A2(1,c) = k;
                c = c + 1;
            end
        end
        %保证选出来的任务没有被安排过
        A2 = setdiff(A2,bused);
        [~,n1] = size(A2);
        val = A2(1,randi(n1)); %在可以加工的任务中随机选取选取一个任务
        p1 = [];
        if ismember(val,A1) == 1 %若所选任务存在or关系
            c = find(P0(val,:) == -1); % 找到val or约束的位置
            [~,n2] = size(c);         % or约束对象个数
            %消除or关系影响
            for k = 1:n2
                p1 = P0(:,c(1,k));
                p1(p1 == -1) = 0;
                P0(:,c(1,k)) = p1;
            end
            %消除and关系影响
            p1 = P0(val,:);  %将val对其后继零件影响消除
            p1 = 0;
            P0(val,:) = p1;
        else
            %消除and关系影响
            p1 = P0(val,:);  %将val对其后继零件影响消除
            p1 = 0;
            P0(val,:) = p1;
        end
        Chrom(i,j) = val;
        % 判断选出的工件属性是

        bused(1,j) = val;   %将已选择的任务进行保存，保证下次不会选中
    end
end
P0 = P;
indices_ca = find(DA == 1);  % 需求操作
indices_ha = find(HA == 1);  % 危险操作
indices = [indices_ca,indices_ha];    %所有有需求和危险的操作
len = length(indices);
A = indices;
dis = [];    %拆卸集合
while ~isempty(A2)
    len = length(indices);
    for i = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
        temp = indices(i);   %操作
        A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
        A1 = A1';
        A = [A,A1];
    end
    A = unique(A);   %必须拆卸的步骤
    dis = [dis,A];
    %获得前置步骤的前置步骤
    A2 = A(~ismember(A, indices));
    indices = A2;
    A = indices;
end
dis = unique(dis, 'rows', 'stable');   %必须拆卸的步骤
for j = 1:np
    notdis = setdiff(1:n, dis);               %其余步骤
    len = length(notdis);
    if len>0
        random_number = randi([1, len]); % 生成 0 到列表长度的随机数
        indices = randperm(len, random_number); % 生成 random_number 个不重复的随机索引
        random_dis = notdis(indices); % 随机取出的拆卸操作
        B = random_dis;               %重新存储一次随机取出的拆卸操作
        random_dis_fm = random_dis;
        A2 = [0];    %保证进入循环
        while ~isempty(A2)
            len = length(random_dis);
            for i = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                temp = random_dis(i);   %操作
                A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
                A1 = A1';
                A = [A,A1];
            end
            A = unique(A);   %必须拆卸的步骤
            random_dis_fm = [random_dis_fm,A];
            %获得前置步骤的前置步骤
            A2 = A(~ismember(A, random_dis));
            random_dis = A2;
            A = random_dis;
        end
        random_dis_fm = unique(random_dis_fm, 'rows', 'stable');   %必须拆卸的步骤
        ChromOne = Chrom(j,:);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
        [~, idx2] = ismember(random_dis_fm, ChromOne);  %随机的拆卸步骤
        Chrom(j,n+idx2) = 1;
    else
        ChromOne = Chrom(j,:);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
    end
end
end

%% 计算目标值
function [Chrom,fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh)
[np,~] = size(Chrom);
fit = zeros(np,3);          %目标值

for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    Temp = Chromone(indices);   %需要拆卸的步骤
    len = length(Temp);
    Profit =  sum(Pr(Temp));   %拆卸得到的利润
    %对需要拆卸的步骤进行解码
    g = 1;           %工厂号
    Mt = zeros(100,1);  %工厂时间
    Mr = 0;   %人工站数目
    Mm = 0;   %机器站数目
    cost = 0; %成本
    profit = 0;   %利润
    ec = 0;       %能耗
    MM=zeros(15,15);
    MT=zeros(15,15);
    MP=zeros(15,1);
    for j = 1:len
        a = Temp(j);   %拆卸步骤
        da = DA(a);    %需求属性
        ca = CA(a);    %复杂属性
        ha = HA(a);    %危害属性
        t = Td(a);      %拆卸时间
        if ca~= 1 && ha ~= 1     
            if j ==1
                p = 1;
                h=1;
                Mm = Mm + 1;
            end
            newP = 1;
            if newP ~= p                 
                g = g + 1;  
                h=1;
                Mm = Mm + 1;  
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t;  
                ec = ec + Ert*t;      
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t;
                ec = ec + Ert*t;       
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Crt*t;
                    ec = ec - Ert*t;
                    g = g + 1;
                    h=1;
                    Mm = Mm + 1;  
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Crt*t;
                    ec = ec + Ert*t;       %计算能耗
                end
            end
        elseif ca == 1 && ha ~= 1    
            if j == 1
                p = -1;
                h=1;
                Mr = Mr + 1;
            end
            newP = -1;      
            if newP ~= p                
                g = g + 1;  
                h=1;
                Mr = Mr + 1;
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Cpt*t;
                ec = ec + Ept*t;       %计算能耗
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Cpt*t;
                ec = ec + Ept*t;       %计算能耗
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Cpt*t;
                    ec = ec - Ept*t;      
                    g = g + 1;
                    h=1;
                    Mr = Mr + 1;  
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Cpt*t;
                    ec = ec + Ept*t;       %计算能耗
                end
            end
        else             
            if j == 1
                p = 1;
                h=1;
                Mm = Mm + 1;  
            end
            newP = 1;      
            if newP ~= p                   %和前一个工厂不同
                g = g + 1; 
                h=1;
                Mm = Mm + 1;  
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
            else
                Mt(g,1) = Mt(g,1) + t;
                cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
                if Mt(g,1) > C
                    Mt(g,1) = Mt(g,1) - t;
                    cost = cost - Crt*t - Ch*t;  %注意加上处理危险任务的额外成本
                    ec = ec - Ept*t - Eh*t;       %计算加上处理危险任务的额外能耗
                    g = g + 1;
                    h=1;
                    Mm = Mm + 1; 
                    Mt(g,1) = Mt(g,1) + t;
                    cost = cost + Crt*t + Ch*t;  %注意加上处理危险任务的额外成本
                    ec = ec + Ept*t + Eh*t;       %计算加上处理危险任务的额外能耗
                end
            end
        end
        p = newP;          
        %Chrom(i,2*n+j) = p;
        MM(g,h)=Temp(j);
        MT(g,h)=Td(Temp(j));
        MP(g,1)=p;
        h=h+1;
    end
    Mt = Mt(1:g,:);
    % 计算负载均衡
    F = 0;
    Ct = max(Mt);
    EC=0;
    for j = 1:g
        F = F + (Ct - Mt(j,1))^2;
    end
    %计算利润
    profit = Profit - cost - g * Cw - Mm * Cr;
    fit(i,1) = sqrt(F)* 3 / 2;
    fit(i,2) = (-profit)* 2 / 3;
    fit(i,3) = (ec+EC)* 3 / 2;
    PL_Fit(i,1) = sum(fit(i,:))/3;   %加权平均值
end
end

%% 计算IGD指标
function igd_value = calculate_igd(pareto_frontier, pareto_set)

% 计算真实 Pareto 前沿集合的大小
pareto_frontier_size = size(pareto_frontier, 1);

% 计算算法生成的 Pareto 解集的大小
pareto_set_size = size(pareto_set, 1);

% 初始化 IGD 值为 0
igd_value = 0.0;

% 遍历算法生成的 Pareto 解集中的每个解
for i = 1:pareto_set_size
    min_distance = Inf; % 初始化最小距离为正无穷大
    % 计算当前解与真实 Pareto 前沿集合中每个解之间的距离
    for j = 1:pareto_frontier_size
        % 计算 Euclidean 距离
        distance = norm(pareto_set(i, :) - pareto_frontier(j, :));
        % 更新最小距离
        if distance < min_distance
            min_distance = distance;
        end
    end
    % 将最小距离加入 IGD 值
    igd_value = igd_value + min_distance;
end

% 计算平均最小距离
igd_value = igd_value / pareto_set_size;
end


%% Q-learning 领域搜索
function [Chrom,Q] = Q_learning(Chrom,s,n,P,Q,Gamma,Alpha,t,Tmax,Proportion,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,pareto_frontier)
Epsilon = 0.9 - 0.85*t/(Tmax*10);
[state,action] = size(Q);

if s == 0
    i = 1;   %表示此状态
    if Epsilon < rand
        % 找到第一行的最大值及其索引
        [maxValue, maxIndex] = max(Q(1, :));
        % 找到所有最大值的索引
        maxIndices = find(Q(1, :) == maxValue);
        % 如果有多个最大值，随机选择一个索引
        randomIndex = maxIndices(randi(length(maxIndices)));   %动作指引
        if randomIndex == 1
            %领域搜索一
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(1,:));
        elseif randomIndex == 2
            %领域搜索二
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(2,:));
        elseif randomIndex == 3
            %领域搜索三
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(3,:));
        elseif randomIndex == 4
            %领域搜索四
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(4,:));
        elseif randomIndex == 5
            %领域搜索五
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(5,:));
        elseif randomIndex == 6
            %领域搜索六
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(6,:));
        elseif randomIndex == 7
            %领域搜索七
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(7,:));
        elseif randomIndex == 8
            %领域搜索八
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(8,:));
        elseif randomIndex == 9
            %领域搜索九
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(9,:));
        else
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(10,:));
        end
        %更新Q值表
        Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max);

    else
        randomIndex = unidrnd(action);    %动作指引
        if randomIndex == 1
            %领域搜索一
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(1,:));
        elseif randomIndex == 2
            %领域搜索二
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(2,:));
        elseif randomIndex == 3
            %领域搜索三
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(3,:));
        elseif randomIndex == 4
            %领域搜索四
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(4,:));
        elseif randomIndex == 5
            %领域搜索五
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(5,:));
        elseif randomIndex == 6
            %领域搜索六
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(6,:));
        elseif randomIndex == 7
            %领域搜索七
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(7,:));
        elseif randomIndex == 8
            %领域搜索八
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(8,:));
        elseif randomIndex == 9
            %领域搜索九
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(9,:));
        else
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(10,:));
        end
        %更新Q值表
        Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max);
    end
elseif s > 1
    i = 20;
    if Epsilon < rand
        % 找到第一行的最大值及其索引
        [maxValue, maxIndex] = max(Q(1, :));
        % 找到所有最大值的索引
        maxIndices = find(Q(1, :) == maxValue);
        % 如果有多个最大值，随机选择一个索引
        randomIndex = maxIndices(randi(length(maxIndices)));   %动作指引
        if randomIndex == 1
            %领域搜索一
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(1,:));
        elseif randomIndex == 2
            %领域搜索二
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(2,:));
        elseif randomIndex == 3
            %领域搜索三
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(3,:));
        elseif randomIndex == 4
            %领域搜索四
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(4,:));
        elseif randomIndex == 5
            %领域搜索五
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(5,:));
        elseif randomIndex == 6
            %领域搜索六
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(6,:));
        elseif randomIndex == 7
            %领域搜索七
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(7,:));
        elseif randomIndex == 8
            %领域搜索八
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(8,:));
        elseif randomIndex == 9
            %领域搜索九
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(9,:));
        else
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(10,:));
        end
        %更新Q值表
        Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max);

    else
        randomIndex = unidrnd(action);    %动作指引
        if randomIndex == 1
            %领域搜索一
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(1,:));
        elseif randomIndex == 2
            %领域搜索二
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(2,:));
        elseif randomIndex == 3
            %领域搜索三
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(3,:));
        elseif randomIndex == 4
            %领域搜索四
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(4,:));
        elseif randomIndex == 5
            %领域搜索五
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(5,:));
        elseif randomIndex == 6
            %领域搜索六
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(6,:));
        elseif randomIndex == 7
            %领域搜索七
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(7,:));
        elseif randomIndex == 8
            %领域搜索八
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(8,:));
        elseif randomIndex == 9
            %领域搜索九
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(9,:));
        else
            Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(10,:));
        end
        %更新Q值表
        Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max);
    end
else
    for i = 1:20
        if s > 1.5/(20)*(i-1) && s <= 1.5/(20)*(i)
            if Epsilon < rand
                % 找到第一行的最大值及其索引
                [maxValue, maxIndex] = max(Q(i, :));
                % 找到所有最大值的索引
                maxIndices = find(Q(i, :) == maxValue);
                % 如果有多个最大值，随机选择一个索引
                randomIndex = maxIndices(randi(length(maxIndices)));   %动作指引
                if randomIndex == 1
                    %领域搜索一
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(1,:));
                elseif randomIndex == 2
                    %领域搜索二
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(2,:));
                elseif randomIndex == 3
                    %领域搜索三
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(3,:));
                elseif randomIndex == 4
                    %领域搜索四
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(4,:));
                elseif randomIndex == 5
                    %领域搜索五
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(5,:));
                elseif randomIndex == 6
                    %领域搜索六
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(6,:));
                elseif randomIndex == 7
                    %领域搜索七
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(7,:));
                elseif randomIndex == 8
                    %领域搜索八
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(8,:));
                elseif randomIndex == 9
                    %领域搜索九
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(9,:));
                else
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(10,:));
                end
                %更新Q值表
                Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max);
            else
                randomIndex = unidrnd(action);      %动作指引
                if randomIndex == 1
                    %领域搜索一
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(1,:));
                elseif randomIndex == 2
                    %领域搜索二
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(2,:));
                elseif randomIndex == 3
                    %领域搜索三
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(3,:));
                elseif randomIndex == 4
                    %领域搜索四
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(4,:));
                elseif randomIndex == 5
                    %领域搜索五
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(5,:));
                elseif randomIndex == 6
                    %领域搜索六
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(6,:));
                elseif randomIndex == 7
                    %领域搜索七
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(7,:));
                elseif randomIndex == 8
                    %领域搜索八
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(8,:));
                elseif randomIndex == 9
                    %领域搜索九
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(9,:));
                else
                    Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion(10,:));
                end
                %更新Q值表
                Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max);
            end
        end
    end
end
end

%% 领域搜索
function Chrom = Domain_search(Chrom,n,P,DA,HA,Proportion_one)
[np,~] = size(Chrom);
%按照Proportion_one中的比例分开Chrom
start_idx = [0 cumsum(np * Proportion_one)];
Chrom_split = cell(1,4);
for i = 1:length(Proportion_one)
    if i == length(Proportion_one)
        Chrom_split{i} = Chrom(start_idx(i) + 1:end, :);
    else
        Chrom_split{i} = Chrom(start_idx(i) + 1:start_idx(i + 1), :);
    end
end
Chrom1 = Chrom_split{1};
Chrom2 = Chrom_split{2};
Chrom3 = Chrom_split{3};
Chrom4 = Chrom_split{4};
%按照不同的领域搜索策略进行搜索
Chrom1 = Serach1(Chrom1,n);
%按照不同的领域搜索策略进行搜索
Chrom2 = Serach2(Chrom2,n,P);
%按照不同的领域搜索策略进行搜索
Chrom3 = Serach3(Chrom3,n,P);
%按照不同的领域搜索策略进行搜索
Chrom4 = Serach4(Chrom4,n,P,DA,HA);
%将四个领域搜索个体拼接起来
Chrom1 = Chrom1(:,1:n*2);Chrom2 = Chrom2(:,1:n*2);
Chrom3 = Chrom3(:,1:n*2);Chrom4 = Chrom4(:,1:n*2);
Chrom = [Chrom1;Chrom2;Chrom3;Chrom4];
end

%% 领域搜索1
function Chrom = Serach1(Chrom,n)
[np,index] = size(Chrom);
Chrom_son_jc = Chrom(:,1:n);
%找到每个个体需要拆卸的步骤保存起来
Temp = zeros(np,n);
for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    len = length(indices);
    Temp(i,1:len) = Chromone(indices);   %需要拆卸的步骤
end

for i = 1:np
    a = unidrnd(np);
    while a == i  %取出的另一个体不能是当前个体
        a = unidrnd(np);
    end
    p1(1,:) = Chrom_son_jc(i,:);
    p1(2,:) = Chrom_son_jc(a,:);
    %随机生成两个交叉点，从而生成我们需要的交叉片段
    p = randperm(n);
    p = p(1:2);
    p = sort(p);
    b=p(1,1);c=p(1,2);
    n1(1,1) = c - b + 1;   %基因片段基因数
    p01 = p1(1,b:c);
    %交叉片段一在对面顺序提取
    for k = 1:n1(1,1)
        w = 1;
        g = 1;
        while g == 1
            if p1(2,w)==p01(1,k)
                p01(2,k) = w;  %在第二段中的位置
                g = 0;
            else
                w = w + 1;
            end
        end
    end
    p001 = [];
    [~,index] = sort(p01(2,:));
    for k = 1:length(index)
        p001(1,k) = p01(1,index(k)); %按照另一个体中的基因顺序排列
    end
    p1(1,b:c) = p001;
    Chrom_son_jc(i,1:n) = p1(1,:);
end
for i = 1:np
    Chrom_son_jc(i,n+1:n*2) = 0;
end
%根据前面存储的Temp，重新置0，1
for i = 1:np
    ChromOne = Chrom_son_jc(i,:);
    % 获取第一行不等于0的元素的逻辑索引
    nonZeroIndices = Temp(i, :) ~= 0;
    % 使用逻辑索引获取需要拆卸的步骤
    dis = Temp(i, nonZeroIndices);
    [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
    Chrom_son_jc(i,n+idx1) = 1;
end
Chrom = Chrom_son_jc;
end

%% 领域搜索2
function Chrom = Serach2(Chrom,n,P)
A1=[];
Chrom_son=Chrom;
[np,~] = size(Chrom);
% 找到or约束
P0 = P;
%生成or集合
a = 1;
%找到每个个体需要拆卸的步骤保存起来
Temp = zeros(np,n);
for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    len = length(indices);
    Temp(i,1:len) = Chromone(indices);   %需要拆卸的步骤
end
for k = 1:n  %遍历行
    if ismember(-1,P0(k,:)) == 1 %存在or约束,记录行数（前驱零件）
        A1(1,a) = k;
        a = a + 1;
    end
end
for k = 1 : n  %遍历列
    if ismember(-1,P0(:,k)) == 1 %存在or约束,记录列数（后驱零件）
        A1(1,a) = k;
        a = a + 1;
    end
end
for i = 1:np
    P0 = P;
    p01=[];p03=[];p2=[];
    p = Chrom(i,:);
    %随机生成2个变异点
    a = randi(n - 3);
    b = randi(n);
    while b - a <= 2   %保证变异段基因大于等于2
        b = randi(n);
    end
    %两个片段之间的基因为变异片段
    p01 = p(1,1:a);
    p03 = p(1,b:n);
    p1 = zeros(1,length(p01) + length(p03));
    p1(1,1:length(p01)) = p01;
    p1(1,length(p01) + 1:length(p01) + length(p03)) = p03;%需要消除影响的零件矩阵
    p2 = p(1,a + 1:b - 1);          %得到变异片段
    p3 = [];
    %消除变异部分外零件的影响
    p00 = [];
    for j = 1:length(p01) + length(p03)
        val = p1(1,j);
        if ismember(val,A1) == 1 && val~=12 && val ~= 22 %若选取的工件在or集合中
            c = find(P0(val,:) == -1); % 找到val or约束的位置
            [~,n2] = size(c);         % or约束对象个数
            %消除or关系影响
            for k = 1:n2
                p00 = P0(:,c(1,k));
                p00(p00 == -1) = 0;
                P0(:,c(1,k)) = p00;
            end
            %消除and关系影响
            p00 = P0(val,:);  %将val对其后继零件影响消除
            p00 = 0;
            P0(val,:) = p00;
        else
            p00 = P0(val,:);    %消除and影响
            p00 = 0;
            P0(val,:) = p00;
        end
    end
    %变异部分重新排序
    for j = 1:length(p2)
        A2 = [];
        %生成优先任务集合
        c  = 1;
        for k = 1:n
            if any(P0(:,k)) == 0 %k零件没有任何前驱约束
                A2(1,c) = k;
                c = c + 1;
            end
        end
        [~,n0] = size(A2);
        %A2只保留p2中的操作
        A20 = intersect(A2, p2);
        [~,n1] = size(A20);
        f = 1;%保证选出的新任务没有安排过
        while f == 1
            val = A20(1,randi(n1));
            if ismember(val,p3) == 0   %在p3中不存在，也就是没有被安排过
                f = 0;
            end
        end
        p00 = [];
        if ismember(val,A1) == 1%若选取的工件在or集合中
            %消除or关系影响
            c = find(P0(val,:) == -1);
            [~,n0] = size(c);%约束对象个数
            for k = 1:n0
                p00 = P0(:,c(1,k));
                p00(find(p1 == -1)) = 0;
                P0(:,c(1,k)) = p00;
            end
            p00 = P0(val,:);    %消除and影响
            p00 = 0;
            P0(val,:) = p00;
        else
            p00 = P0(val,:);    %消除and影响
            p00 = 0;
            P0(val,:) = p00;
        end
        p3(1,j) = val;
    end
    Chrom_son(i,a+1:b-1) = p3(1,:);
end
%根据前面存储的Temp，重新置0，1
for i = 1:np
    ChromOne = Chrom_son(i,:);
    % 获取第一行不等于0的元素的逻辑索引
    nonZeroIndices = Temp(i, :) ~= 0;
    % 使用逻辑索引获取需要拆卸的步骤
    dis = Temp(i, nonZeroIndices);
    [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
    Chrom_son(i,n+idx1) = 1;
end
Chrom = Chrom_son;
end

%% 领域搜索3
function Chrom = Serach3(Chrom,n,P)
[np,~] = size(Chrom);
Chrom_use = Chrom;
P0 = P;
%找到每个个体需要拆卸的步骤保存起来
Temp = zeros(np,n);
for i = 1:np
    Chromone = Chrom(i,:);
    % 将向量转换为矩阵
    matrix_Chrom = reshape(Chromone, 1, []);
    % 找到需要拆卸的步骤
    indices = find(matrix_Chrom(1, n+1:2*n) == 1);
    len = length(indices);
    Temp(i,1:len) = Chromone(indices);   %需要拆卸的步骤
end
for i  = 1:np
    Chrom_son = Chrom_use(i,1:n);
    len = length(Chrom_son);
    index = unidrnd(len);
    temp = Chrom_son(index);   % 取出变异点基因
    A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
    A2 = find(P0(temp, :) == -1 | P0(temp, :) == 1);  %后置与或集合
    % 找到A1中数字在Chrom_son中的索引位置
    indices = ismember(Chrom_son, A1);
    indexes_A1 = find(indices);
    % 找到A2中数字在Chrom_son中的索引位置
    indices = ismember(Chrom_son, A2);
    indexes_A2 = find(indices);
    if isempty(A1)
        if ~isempty(A2)                  %无前置有后置
            later_indices = indexes_A2(indexes_A2 > index);   %找到后置索引
            later_indices = min(later_indices);
            % 1到later_indices-1为插入区间
            random_number = randi([1, later_indices - 1]);  %插入位置
            if random_number < index
                Chrom_son(random_number+1:index) = Chrom_son(random_number:index-1);
                Chrom_son(random_number) = temp;
            else
                Chrom_son(index:random_number-1) = Chrom_son(index+1:random_number);
                Chrom_son(random_number) = temp;
            end
        end
    elseif isempty(A2)
        if ~isempty(A1)                  %无后置有前置
            former_indices = indexes_A1(indexes_A1 < index); %找到前置索引
            former_indices = max(former_indices);
            % former_indices+1到len是插入区间
            random_number = randi([former_indices + 1, len]); %插入位置
            if random_number < index
                Chrom_son(random_number+1:index) = Chrom_son(random_number:index-1);
                Chrom_son(random_number) = temp;
            else
                Chrom_son(index:random_number-1) = Chrom_son(index+1:random_number);
                Chrom_son(random_number) = temp;
            end
        end
    else
        former_indices = indexes_A1(indexes_A1 < index);
        later_indices = indexes_A2(indexes_A2 > index);
        former_indices = max(former_indices);
        later_indices = min(later_indices);
        % former_indices+1到later_indices-1是插入区间
        if any(former_indices ~= 0) && any(later_indices ~= 0)
            random_number = randi([former_indices + 1, later_indices - 1]); %插入位置
            if random_number < index
                Chrom_son(random_number+1:index) = Chrom_son(random_number:index-1);
                Chrom_son(random_number) = temp;
            else
                Chrom_son(index:random_number-1) = Chrom_son(index+1:random_number);
                Chrom_son(random_number) = temp;
            end
        end
    end
    Chrom_use(i,1:n) = Chrom_son;
end
%根据前面存储的Temp，重新置0，1
for i = 1:np
    ChromOne = Chrom_use(i,:);
    % 获取第一行不等于0的元素的逻辑索引
    nonZeroIndices = Temp(i, :) ~= 0;
    % 使用逻辑索引获取需要拆卸的步骤
    dis = Temp(i, nonZeroIndices);
    [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
    Chrom_use(i,n+idx1) = 1;
end
Chrom = Chrom_use;
end

%% 领域搜索4
function Chrom = Serach4(Chrom,n,P,DA,HA)
%重新选择是否拆卸
[np,~] = size(Chrom);
P0 = P;
indices_ca = find(DA == 1);  % 需求操作
indices_ha = find(HA == 1);  % 危险操作
indices = [indices_ca,indices_ha];    %所有有需求和危险的操作
len = length(indices);
A = indices;
dis = [];    %拆卸集合
A2 = [0]; %保证进入循环
while ~isempty(A2)
    len = length(indices);
    for i = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
        temp = indices(i);   %操作
        A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
        A1 = A1';
        A = [A,A1];
    end
    A = unique(A);   %必须拆卸的步骤
    dis = [dis,A];
    %获得前置步骤的前置步骤
    A2 = A(~ismember(A, indices));
    indices = A2;
    A = indices;
end
dis = unique(dis, 'rows', 'stable');   %必须拆卸的步骤
for j = 1:np
    notdis = setdiff(1:n, dis);               %其余步骤
    len = length(notdis);
    if len > 0
        random_number = randi([1, len]); % 生成 0 到列表长度的随机数
        indices = randperm(len, random_number); % 生成 random_number 个不重复的随机索引
        random_dis = notdis(indices); % 随机取出的拆卸操作
        B = random_dis;               %重新存储一次随机取出的拆卸操作
        random_dis_fm = random_dis;
        A2 = [0];    %保证进入循环
        while ~isempty(A2)
            len = length(random_dis);
            for i = 1:len   % 找出有需求和危险操作的所有前置拆卸操作
                temp = random_dis(i);   %操作
                A1 = find(P0(:, temp) == -1 | P0(:, temp) == 1);  %前置与或集合
                A1 = A1';
                A = [A,A1];
            end
            A = unique(A);   %必须拆卸的步骤
            dis = [dis,A];
            %获得前置步骤的前置步骤
            A2 = A(~ismember(A, random_dis));
            random_dis = A2;
            A = random_dis;
        end
        random_dis_fm = unique(random_dis_fm, 'rows', 'stable');   %必须拆卸的步骤
        ChromOne = Chrom(j,1:n);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
        [~, idx2] = ismember(random_dis_fm, ChromOne);  %随机的拆卸步骤
        Chrom(j,n+idx2) = 1;
    else
        ChromOne = Chrom(j,1:n);
        [~, idx1] = ismember(dis, ChromOne);     %必须拆卸步骤
        Chrom(j,n+idx1) = 1;
    end
end
end

%% 更新Q值表
function Q = flashQ(Q,i,randomIndex,Gamma,Alpha,pareto_frontier,Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh,IGD_t_sum,IGD_t_max,IGD_one_S2,IGD_one_sum,IGD_one_max)
%计算目标值
[np,~] = size(Chrom);
[Chrom,new_fit,PL_Fit] = Cal(Chrom,n,Td,C,DA,CA,HA,Pr,Cw,Cr,Cpt,Crt,Ch,Ept,Ert,Eh);  %计算目标值
%计算这一代的IGD和以及IGD最大值
new_IGD_t_sum = 0;
new_IGD_t_max = 0;
%把第一代个体的IGD值求和
for j = 1:np
    new_IGD_t_sum = new_IGD_t_sum + calculate_igd(pareto_frontier,new_fit(j,:));
    if calculate_igd(pareto_frontier,new_fit(j,:)) > new_IGD_t_max
        new_IGD_t_max = calculate_igd(pareto_frontier,new_fit(j,:));
    end
end
s1 = new_IGD_t_sum/IGD_one_sum;
% 将种群多样性进行归一化
new_IGD_t_S2 = 0;
for j = 1:np
    new_IGD_t_S2 = new_IGD_t_S2 + abs(calculate_igd(pareto_frontier,new_fit(j,:)) - new_IGD_t_sum/np);
end
s2 = new_IGD_t_S2/IGD_one_S2;     %种群多样性进行归一化结果
%将最佳适应度归一化
s3 = new_IGD_t_max / IGD_one_max;
news = 0.35*s1 + 0.35*s2 + 0.3*s3;
% if news > Max
%    disp(news);
% end
%计算奖励值
R = ((IGD_t_sum - new_IGD_t_sum)/IGD_t_sum)  + ((IGD_t_max - new_IGD_t_max)/IGD_t_max) ;
if news == 0
    Q(i,randomIndex) =  Q(i,randomIndex) + Alpha*(R + Gamma * max(Q(1,:) - Q(i,randomIndex)));
elseif news>1
    j = 20;
    Q(i,randomIndex) =  Q(i,randomIndex) + Alpha*(R + Gamma * max(Q(j,:)) - Q(i,randomIndex));
else
    for j = 1:20
        if news > 1.5/(20)*(j-1) && news <= 1.5/(20)*(j)
            Q(i,randomIndex) =  Q(i,randomIndex) + Alpha*(R + Gamma * max(Q(j,:)) - Q(i,randomIndex));
            break;
        end
    end
end
end



%% frontvalue = Newranking(ObjV)                         %非支配排序
function frontvalue = Newranking(ObjV)                         %非支配排序

fnum=0;                                                        %当前分配的前沿面编号
cz=false(1,size(ObjV,1));                              %记录个体是否已被分配编号
frontvalue=zeros(size(cz));                            %每个个体的前沿面编号
[functionvalue_sorted,newsite]=sortrows(ObjV);         %对种群按第一维目标值大小进行升序排序
while ~all(cz)                                         %开始迭代判断每个个体的前沿面,采用改进的deductive sort
    fnum=fnum+1;
    d=cz;
    for i=1:size(ObjV,1)
        if ~d(i)
            for j=i+1:size(ObjV,1) %判断i对应的所有集合里面的支配和非支配的解，被i支配则为1，不被i支配则为0
                if ~d(j)
                    k=1;
                    for m=2:size(ObjV,2) %判断是否支配，找到个体p不支配的个体，标记为k=0
                        if functionvalue_sorted(i,m)>functionvalue_sorted(j,m)%判断i,j是否支配，如果成立i,j互不支配
                            k=0;%i,j互不支配
                            break
                        end
                    end
                    if k
                        d(j)=true;%那么p所支配的个体k=1并记录在d里，则后面该个体已被支配就不能在这一层里进行判断
                    end
                end
            end
            frontvalue(newsite(i))=fnum;%实际位置的非支配层赋值
            cz(i)=true;
        end
    end
end
end

%% Selch=select1('rws',ChromNumber,Fitnv);  %根据Fitnv适应度值轮盘赌进行选择
function Chrom=select(NIND,frontvalue,FChrom,Objv)
%---计算拥挤距离/选出下一代个体==============================================zxs
fnum=0;                                                                %当前前沿面
while numel(frontvalue,frontvalue<=fnum+1)<NIND              %判断前多少个面的个体能完全放入下一代种群
    fnum=fnum+1;
end

newnum=numel(frontvalue,frontvalue<=fnum);                     %前fnum个面的个体数
Chrom(1:newnum,:)=FChrom(frontvalue<=fnum,:);                  %将前fnum个面的个体复制入下一代
popu=find(frontvalue==fnum+1);                                 %popu记录第fnum+1个面上的个体编号
distancevalue=zeros(size(popu));                               %popu各个体的拥挤距离
fmax=max(Objv(popu,:),[],1);                                   %popu每维上的最大值
fmin=min(Objv(popu,:),[],1);                                   %popu每维上的最小值
for i=1:size(Objv,2)                                           %分目标计算每个目标上popu各个体的拥挤距离
    [~,newsite]=sortrows(Objv(popu,i));                        %popu里对第一维排序之后的位置
    distancevalue(newsite(1))=inf;
    distancevalue(newsite(end))=inf;
    for j=2:length(popu)-1
        distancevalue(newsite(j))=distancevalue(newsite(j))+(Objv(popu(newsite(j+1)),i)-Objv(popu(newsite(j-1)),i))/(fmax(i)-fmin(i));
    end
end
popu=-sortrows(-[distancevalue;popu]')';                             %按拥挤距离降序排序第fnum+1个面上的个体
Chrom(newnum+1:NIND,:)=FChrom(popu(2,1:NIND-newnum),:);          %将第fnum+1个面上拥挤距离较大的前popnum-newnum个个体复制入下一代
end

%%非支配解分层排序
function [frontlabel]=quickrank(lastpopution,frontvalue)
lastpopution=[zeros((size(lastpopution,1)),1),lastpopution];
frontlabel=[];
ceng=0;
while ~isempty(lastpopution)
    %      size(lastpopution,1)
    %      pause(2)
    label=[];
    ceng=ceng+1;
    if size(lastpopution,1)>1
        for i1=1:size(lastpopution,1)
            n=0;
            pnum=[];
            qnum=[];
            for j1=1:size(lastpopution,1)
                if i1==j1
                    continue;%结束当前循环,不往下执行
                end
                chazhi=frontvalue(i1,:)-frontvalue(j1,:);
                if (length(find(chazhi(1,:)>=0))==3)&(length(find(chazhi(1,:)>0))>=1)    %%被别人支配
                    n=1;
                    break;
                end
                if  ((length(find(chazhi(1,:)<=0))==3)&(length(find(chazhi(1,:)<0))>=1))||(((length(find(chazhi(1,:)>0))>=1))&(length(find(chazhi(1,:)<0))>=1))
                    plabel=j1;%%%支配解
                    pnum=[pnum; plabel]; %%%支配解标签
                end
                if (length(find(chazhi(1,:)==0))==3)  %%相同点
                    qlabel=j1;  %%%相同点标签
                    qnum=[qnum;qlabel];  %%%相同点-------------------标记相同点的
                end
            end
            pnum1=length(pnum);  %66
            qnum1=length(qnum);
            zonggeshu=pnum1+qnum1;
            if  (n==0)&(qnum1==0)
                lastpopution(i1,1)=ceng;
                label1=i1;
                label=[label;label1];
            end
            if  (n==0)&(zonggeshu==(size(lastpopution,1)-1))&(qnum1>0)
                %qnum=[qnum;i1];
                lastpopution(i1,1)=ceng;
                label1=i1;
                label=[label;label1];
            end
        end
        zhongjian1=lastpopution(label,:);
        frontlabel=[frontlabel;zhongjian1];
        lastpopution(label,:)=[];
        frontvalue(label,:)=[];
    else
        lastpopution(1,1)=ceng;
        zhongjian1=lastpopution(1,:);
        frontlabel=[frontlabel;zhongjian1];
        lastpopution(1,:)=[];
        frontvalue(label,:)=[];
    end
end
end

%%只是为了得到非支配解
function on_inferior=feizhipei(lastpopution,fitvalue)
%一个解不受任何解支配就是非支配解
%目前是：最小时间，最大利润，最小能耗
%值相减，有正有负，则互不支配，存在都是正的就是存在被支配的解，均为0则去除重复点。
%目前对目标函数2进行了取倒数，处理，就都是求最小了。
simi_m=[];  %存储一系列相同值
on_inferior=[];  %存储非支配解
for i=1:size(lastpopution,1)
    if length(find(simi_m==i))==0
        label=1:size(lastpopution,1);
        fit_value=fitvalue(i,:);  %函数值
        label(i)=[];
        fit_value_1=repmat(fit_value,length(label),1);
        chazhi=fit_value_1-fitvalue(label,:); %差值
        on_inferior_label=1;
        simi_position=[]; %存储相同值的位置
        for j=1:length(label)
            if chazhi(1)<0.000000001&&chazhi(1)>-0.000000001
                chazhi(1)=0;
            end
            if chazhi(2)<0.000000001&&chazhi(2)>-0.000000001
                chazhi(2)=0;
            end
            if chazhi(3)<0.000000001&&chazhi(3)>-0.000000001
                chazhi(3)=0;
            end
            if length(find(chazhi(j,:)==0))==3
                simi_position=[simi_position,label(j)];
            end
            if (length(find(chazhi(j,:)>=0))==3)&&(length(find(chazhi(j,:)>0))>1)
                on_inferior_label=0;
                break
            end
        end
        if  on_inferior_label==1
            if length(simi_position)>0
                simi_m=[simi_m,simi_position];
            end
            on_inferior=[on_inferior,i];
        end
    end
end
on_inferior=lastpopution(on_inferior,:);
end

%%%拥挤度的计算
function population=crowded(on_inferior,popunumber)
functionvalue=on_inferior(:,1:3);
distancevalue=zeros(size(functionvalue,1),1);
fmax=max(on_inferior(:,1:3),[],1);                       %三个目标上的最大值
fmin=min(on_inferior(:,1:3),[],1);                          %三个目标上的最小值
%对数据进行归一化
for i=1:size(functionvalue,1)
    functionvalue(i,1)=(functionvalue(i,1)-fmin(1))/(fmax(1)-fmin(1));
    functionvalue(i,2)=(functionvalue(i,2)-fmin(2))/(fmax(2)-fmin(2));
    functionvalue(i,3)=(functionvalue(i,3)-fmin(3))/(fmax(3)-fmin(3));
end
for l=1:3     %分目标计算每个目标上popu各个体的拥挤距离
    [~,newsite]=sortrows(functionvalue(:,l));
    distancevalue(newsite(1))=inf;
    distancevalue(newsite(end))=inf;
    for m=2:(size(functionvalue,1)-1)
        distancevalue(newsite(m))=distancevalue(newsite(m))+(functionvalue(newsite(m+1),l)-functionvalue(newsite(m-1),l))/(fmax(l)-fmin(l));
    end
end
%         distancevalue
distance_value=zeros(size(functionvalue,1),1);
for l=1:3
    [~,newsite]=sortrows(functionvalue(:,l));
    for m=2:(size(functionvalue,1)-1)
        distance_value(newsite(m))= distance_value(newsite(m))+(abs((functionvalue(newsite(m+1),l)-functionvalue(newsite(m-1),l)))-distancevalue(newsite(m)))^2;
    end
end
distance_value=distance_value/3;
[~,newsite1]=sortrows(distance_value);  %升序
newsite1=fliplr(newsite1');  %反转
newsite1=newsite1(1:popunumber);
population=on_inferior(newsite1,4:end);
end