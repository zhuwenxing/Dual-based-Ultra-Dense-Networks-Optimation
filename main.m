%% Subgradient Method to slove resource allocation in small cell networks
% Li Zhou 2014/5/24

showflag = 0;
type = 1;
% problem constants
num_sm = 4; %small cell number
usernum_sm = 4; %number of users in each small cell
chnum = 32; % channel number
usernum_ma = usernum_sm; %number of users in macro cell, set it equal to small cell for simplicity
eta = 0; % energy efficiency factor
MAX_ITERS_ETA = 50;  % max iterations for searching eta
MAX_ITERS_SUB = 100; % max iterations for subgradient search
r_ma = 289; %radius of macro cell,m
r_sm = 40; %radius of small cell,m
totalsunum = num_sm*usernum_sm; %total su number
totalusernum = usernum_ma+totalsunum; %total user number
Pm = 40;%W, 46dBm = 40W
Ps = 10; %W, 30dBm = 10W
%refer to 'Power consumption modeling of different base station types in
%heterogeneous cellular networks, 2010'
% P = [40 10 10]; % maximum power
P = zeros(1,num_sm+1);
P(1) = 40;
P(2:num_sm+1) = 10;% initial value of maximum power
% Pc = [54.8 15 15]; % circuit power
Pc = zeros(1,num_sm+1);
Pc(1) = 54.8;
Pc(2:num_sm+1) = 15;% initial value of mu
% gamma = [2.63 5 5]; % power-amplifier inefficiency factor of BS i
gamma = zeros(1,num_sm+1);
gamma(1) = 2.63;
gamma(2:num_sm+1) = 5;% initial value of gamma
Gamma =  1.9972; % SINR gap
alpha = 1;% fairness parameter
% beta = [1 1 1]; % QoS weight
beta = zeros(1,num_sm+1);
beta(1) = 1;
beta(2:num_sm+1) = 1;% initial value of beta
w = zeros(1+num_sm,usernum_sm); % weight for each user
w = w + 1; % set initial value to be 1
Noise=10^((-174+10*log10(10e7/chnum))/10); % white noise
%% for monte carlo simulation
% MonteCarloTimes = 50;
HistoryBestTra = []; 
HistoryEPCap = [];
HistoryEPTra = [];
HistoryLag = [];
HistoryCap = [];
HistoryTra = [];
Historyeta = [];
Historypx = [];
HistoryBestSumP = [];
HistoryBestCap = [];
% for monte = 1:MonteCarloTimes 
%% create scenario
ap = CreateAPs(num_sm,showflag);
user = CreateUsers(num_sm,usernum_ma,usernum_sm,ap,showflag);
% compute pathloss and channel gain from the macro basestation
dis = zeros(1+num_sm,totalusernum); % distance to donor BS
PL = zeros(1+num_sm,totalusernum); % pathloss
Loss = zeros(1+num_sm,totalusernum); % pathloss factor
Gain = zeros(1+num_sm,totalusernum,chnum); %channel gain
g = zeros(1+num_sm,totalusernum,chnum); %channel gain with pathloss
for i = 1:1+num_sm
    for j = 1:totalusernum
        dis(i,j) = sqrt((user(1,j)-ap(1,i))^2+(user(2,j)-ap(2,i))^2);
        if i == 1
            PL(i,j) = MacroPathloss(dis(i,j));
        else
            PL(i,j) = pathloss(dis(i,j));
        end
        Loss(i,j) = 1/(10^(PL(i,j)/10));
        Gain(i,j,:) = raylrnd(1:chnum);
%         Gain(i,j,:) = 1;
        g(i,j,:) = Gain(i,j,:)*Loss(i,j);
    end
end
% save 'scenariodata.mat'
% load('scenariodata.mat');
tStart_eta = tic;
%% search eta
m = 1;
while m <= MAX_ITERS_ETA 
    %% initialization
    x = zeros(1+num_sm,usernum_sm,chnum);% channel allocation indicator
    x = x + 1;
    p = zeros(1+num_sm,usernum_sm,chnum); % power allocation policy
    p0 = zeros(1+num_sm,usernum_sm,chnum);
    pbest = zeros(1+num_sm,usernum_sm,chnum);
    p_vec = zeros(1,(1+num_sm)*usernum_sm*chnum);
    p0_vec = zeros(1,(1+num_sm)*usernum_sm*chnum);
    p0_vec = p0_vec + 2; % initial value of p if fsolve
    mu = zeros(1,num_sm+1);
    mu(1) = 40;
    mu(2:num_sm+1)=10;% initial value of mu
%     mu = [40 10 10]; 
    lambda = zeros(num_sm+1,usernum_sm,chnum);
    lambda_star = zeros(num_sm+1,chnum);
    lambda_star_best = zeros(num_sm+1,chnum);
    step = []; % save history steps
    Lagrangian = []; % save history Lagrangians
    TotalCapacity = []; % save history total capacity of the network
    TotalTradeoffValue = []; % save history total trade off value of the network
    Lbest = 3000;
    BestCapacity = 0;
    etabest = 0;
    rhobest = zeros(1,num_sm);
    record_p = [];
    record_px = [];
    record_fval = [];
    record_mu = [];   
    t = 1;
    %% subgradient search
    tStart_sub = tic;
    while t <= MAX_ITERS_SUB 
        fun = @(p_vec)powerfun1(p_vec,num_sm,usernum_sm,chnum,g,x,w,mu,eta,gamma,Gamma,Noise);
%         options = optimoptions(@fsolve,'Algorithm', 'Levenberg-Marquardt','MaxIter',5000,'MaxFunEvals',2000,'TolFun',1e-6,'TolX',1e-6);%,'Display','iter-detailed');
        options = optimset('Algorithm', 'Levenberg-Marquardt','MaxIter',5000,'MaxFunEvals',2000,'TolFun',1e-6,'TolX',1e-6);
        [p_vec,fval, exitflag, output] = fsolve(fun,p0_vec,options);
        fval;
        record_fval = [record_fval fval];
        p_vec = max(p_vec,0);
        tt = 0;
        for n = 1:chnum
            for k = 1:usernum_sm
                for i = 1:num_sm+1
                    tt = tt + 1;
                    p(i,k,n) = p_vec(tt);
                end
            end
        end
        % sum(sum(sum(p(:,:,:))));
        record_p = [record_p sum(sum(sum(p(:,:,:))))];
        %%  caculate lambda_star
        for i = 1:num_sm+1
            for n = 1:chnum
                for k = 1:usernum_sm
                    I = 0;
                    globe_k = (i-1)*usernum_sm+k;
                    for j = 1:num_sm+1
                        if j~=i
                            for q = 1:usernum_sm
                                I = I + p(j,q,n)*g(j,globe_k,n);
                            end
                        end
                    end
                    SINR = p(i,k,n)*g(i,globe_k,n)/(I+Noise);
                    lambda(i,k,n) = w(i,k)*log2(1+Gamma*SINR)-eta*gamma(i)*p(i,k,n)-mu(i)*p(i,k,n);
                end
                [lambda_star(i,n),k_star] = max(lambda(i,:,n));
                x(i,:,n) = zeros(1,usernum_sm);
                x(i,k_star,n) = 1;
            end
        end

       %% caculate the power allocation after channel allocation
        p = p.*x;
       %% make p fits the constraints
        for i = 1:num_sm+1
            if sum(sum(p(i,:,:))) > P(i)
                p(i,:,:) = p(i,:,:)*P(i)/sum(sum(p(i,:,:)));
            end
        end
        record_px = [record_px sum(sum(sum(p(:,:,:))))];
        TempCapacity = 0;
        for i = 1:num_sm+1
            for k = 1:usernum_sm
                for n = 1:chnum
                    %caculate f1
                    tempI = 0;
                    for j = 1
                        if j ~= i
                            for q = 1:usernum_sm
                                tempI = tempI + p(j,q,n)*g(j,(i-1)*usernum_sm+k,n);
                            end
                        end
                    end
                    tempI = tempI + Noise;% + p(i,k,n)*g(i,(i-1)*usernum_sm+k,n);
                    TempCapacity = TempCapacity + x(i,k,n)*log2(1+Gamma*p(i,k,n)*g(i,(i-1)*usernum_sm+k,n)/tempI)/chnum;
                end
            end
        end
        TotalCapacity = [TotalCapacity TempCapacity];
        TempTradeoffValue = TempCapacity;
        for i = 1:num_sm+1
            TempTradeoffValue = TempTradeoffValue - gamma(i)*sum(sum(x(i,:,:).*p(i,:,:)))-Pc(i);
        end
        TempTradeoffValue;
        TotalTradeoffValue = [TotalTradeoffValue TempTradeoffValue];
       %% caculate Lagrangian
        tempLagrangian = sum(mu.*P) - eta*sum(Pc) + sum(sum(lambda_star));
%         temp = 0;
%         for i = 1:num_sm+1
%             temp = temp + (eta*gamma(i)+mu(i))*sum(sum(p(i,:,:))) + eta*Pc(i);
%         end
%         tempLagrangian = TempCapacity + sum(mu.*P) - temp;
        Lagrangian = [Lagrangian tempLagrangian];
        %% subgradient update
        %     if t < 1000
        %         tempstep = 0.1/((P(1)-sum(sum(p(1,:,:))))^2+(P(2)-sum(sum(p(2,:,:))))^2+(P(3)-sum(sum(p(3,:,:))))^2);
        %     else
        %         tempstep = abs(Lbest - tempLagrangian)/((P(1)-sum(sum(p(1,:,:))))^2+(P(2)-sum(sum(p(2,:,:))))^2+(P(3)-sum(sum(p(3,:,:))))^2)
        %     end
%         tempstep = 0.1/((P(1)-sum(sum(p(1,:,:))))^2+(P(2)-sum(sum(p(2,:,:))))^2+(P(3)-sum(sum(p(3,:,:))))^2+1);
%         tempstep = 10/sqrt(t); %step1
        stepsize = zeros(1,num_sm+1);
        e = 0.1;
        for i = 1:num_sm+1
            stepsize(i) = e/((P(i)-sum(sum(p(i,:,:))))^2 + e);
            mu(i) = max(mu(i) - stepsize(i)*(P(i)-sum(sum(p(i,:,:)))),0.00001);
        end
%         step = [step tempstep];
%         mu(1) = max(mu(1) - tempstep*(P(1)-sum(sum(p(1,:,:)))),0.00001);
%         mu(2) = max(mu(2) - tempstep*(P(2)-sum(sum(p(2,:,:)))),0.00001);
%         mu(3) = max(mu(3) - tempstep*(P(3)-sum(sum(p(3,:,:)))),0.00001);
        record_mu = [record_mu; mu];
        %% Lagrangian update
%         if sum(sum(p(1,:,:)))<=P(1) && sum(sum(p(2,:,:)))<=P(2) && sum(sum(p(3,:,:)))<=P(3) % if the current point is feasible
            if tempLagrangian<Lbest  % update the beset Lagrangian
                Lbest = tempLagrangian;
                pbest = p;
                xbest = x;
                tbest = t;
                mubest = mu;
            end
            if TempCapacity > BestCapacity  % update the best capacity
                BestCapacity = TempCapacity;
                tCapbest = t;
            end
%         end
        % print time
        tElapsed_sub = toc(tStart_sub);
        tRemain_sub = tElapsed_sub*(MAX_ITERS_SUB/t-1);
        tTotal_sub = tElapsed_sub + tRemain_sub;
        fprintf('SUB : %d/%d  %d/%d  %f seconds have elapsed, %f seconds remain\n',t,MAX_ITERS_SUB,m,MAX_ITERS_ETA,tElapsed_sub,tRemain_sub+tTotal_sub*(MAX_ITERS_ETA-m));
        t = t + 1; 
    end
    HistoryLag = [HistoryLag;Lagrangian];
    HistoryCap = [HistoryCap;TotalCapacity];
    HistoryTra = [HistoryTra;TotalTradeoffValue];
    Historypx = [Historypx;record_px];
    %% caculate the best capacity
    BestCapacity = 0;
    for i = 1:num_sm+1
        for k = 1:usernum_sm
            for n = 1:chnum
                %caculate f1
                tempI = 0;
                for j = 1
                    if j ~= i
                        for q = 1:usernum_sm
                            tempI = tempI + pbest(j,q,n)*g(j,(i-1)*usernum_sm+k,n);
                        end
                    end
                end
                tempI = tempI + Noise;
                BestCapacity = BestCapacity + xbest(i,k,n)*log2(1+Gamma*pbest(i,k,n)*g(i,(i-1)*usernum_sm+k,n)/tempI)/chnum;
            end
        end
    end
    BestTradeoffValue = BestCapacity;% initial value
    BestSumP = 0;
    for i = 1:num_sm+1
        tempBestSumP = gamma(i)*sum(sum(xbest(i,:,:).*pbest(i,:,:)))+Pc(i);
        BestSumP = BestSumP + tempBestSumP;
        BestTradeoffValue = BestTradeoffValue - eta*tempBestSumP;
    end
    HistoryBestTra = [HistoryBestTra BestTradeoffValue];
    HistoryBestSumP = [HistoryBestSumP BestSumP];
    HistoryBestCap = [HistoryBestCap BestCapacity];
    %% caculate equal power capacity using xbest
    Ep = zeros(1+num_sm,usernum_sm,chnum);
    for i = 1:num_sm+1
        Ep(i,:,:) = Ep(i,:,:)+P(i)/chnum;
    end
    Ep = Ep.*x;
    EPCapacity = 0;
    for i = 1:num_sm+1
        for k = 1:usernum_sm
            for n = 1:chnum
                %caculate f1
                tempI = 0;
                for j = 1
                    if j ~= i
                        for q = 1:usernum_sm
                            tempI = tempI + Ep(j,q,n)*g(j,(i-1)*usernum_sm+k,n);
                        end
                    end
                end
                tempI = tempI + Noise;% + Gamma*Ep(i,k,n)*g(i,(i-1)*usernum_sm+k,n);
                EPCapacity = EPCapacity + xbest(i,k,n)*log2(1+Gamma*Ep(i,k,n)*g(i,(i-1)*usernum_sm+k,n)/tempI)/chnum;
            end
        end
    end
    HistoryEPCap = [HistoryEPCap EPCapacity];
    EPTradeoffValue = EPCapacity;
    for i = 1:num_sm+1
        EPTradeoffValue = EPTradeoffValue - eta*(gamma(i)*sum(sum(xbest(i,:,:).*Ep(i,:,:)))+Pc(i));
    end
    HistoryEPTra = [HistoryEPTra EPTradeoffValue];   
    %% update eta
    Historyeta = [Historyeta eta];
    eta = BestCapacity/BestSumP;
    %% print iteration and time
    tElapsed_eta = toc(tStart_eta);
    tRemain_eta = tElapsed_eta*(MAX_ITERS_ETA/m-1);
    fprintf('ETA : %d/%d    %f seconds have elapsed, %f seconds remain\n',m,MAX_ITERS_ETA,tElapsed_eta,tRemain_eta);
    m = m + 1;
end % for MAX_ITER_ETA
%% save data and plot figure
save('data6162.mat')
% figure
% h = mycdfplot(BestCapacityData,1);
% hold on
% h = mycdfplot(EPCapacityData,2);
% legend('Best','EP')
% 
% figure
% subplot(3,1,1)
% plot(Lagrangian,'r')
% hold on
% plot(tbest*ones(1,141),-70:1:70,'b')
% title('Lagrangian')
% subplot(3,1,2)
% plot(TotalCapacity,'b')
% hold on
% plot(EPCapacity*ones(1,length(TotalCapacity)),'k')
% plot(tbest*ones(1,10),1:1:10,'r')
% title('Capacity')
% subplot(3,1,3)
% plot(TotalTradeoffValue,'b')
% hold on
% plot(EPTradeoffValue*ones(1,length(TotalTradeoffValue)),'k')
% plot(tCapbest*ones(1,20),1:1:20,'r')
% title('TradeoffValue')
% figure 
% subplot(3,1,1)
% plot(step)
% title('step')
% subplot(3,1,2)
% plot(record_p)
% title('record-p')
% subplot(3,1,3)
% plot(record_px)
% title('record-px')
% figure 
% subplot(3,1,1)
% plot(record_mu(:,1))
% title('mu(1)')
% subplot(3,1,2)
% plot(record_mu(:,2))
% title('mu(2)')
% subplot(3,1,3)
% plot(record_mu(:,3))
% title('mu(3)')