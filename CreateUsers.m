function user = CreateUsers( num_sm,usernum_ma,usernum_sm,ap,showflag )
%CREATEUSERS Summary of this function goes here
%   Detailed explanation goes here
% generate macro cell users
r_ma = 289; %radius of macro cell,m
r_sm = 40; %radius of small cell,m
x_mu = zeros(1,usernum_ma);
y_mu = zeros(1,usernum_ma);
x_su = zeros(num_sm,usernum_sm);
y_su = zeros(num_sm,usernum_sm);
mu = []; %mu coordinate
su = []; %su coordinate
user = []; %user coordinate
x_sm = ap(1,2:num_sm+1);
y_sm = ap(2,2:num_sm+1);
if usernum_ma ~=0
    for i = 1:usernum_ma
        [x_mu(i),y_mu(i)] = random_generate(r_ma,0,0);
        flag = 0;
        while flag == 0
            flag = 1;
            for j = 1:num_sm
                if distance(x_mu(i),y_mu(i),x_sm(j),y_sm(j))<r_sm*1.05
                    [x_mu(i),x_mu(i)] = random_generate(r_ma,0,0);
                    flag = 0;
                end
            end
        end
        mu = [mu [x_mu(i);y_mu(i)]];
        if showflag == 1
            plot(x_mu(i),y_mu(i),'*');   
            str =  sprintf('%d ',i);
% 	        text(x_mu(i),y_mu(i),str)
        end
    end
    user = [user mu];
end
% generate small cell users
for i = 1:num_sm
    str =  sprintf('%d ',i);
% 	text(ap(1,i+1),ap(2,i+1),str)
    for j = 1:usernum_sm
        [x_su(i,j),y_su(i,j)] = random_generate(r_sm,x_sm(i),y_sm(i));
        if showflag == 1
            plot(x_su(i,j),y_su(i,j),'*');  
            str =  sprintf('%d ',usernum_ma+(i-1)*usernum_sm+j);
% 	        text(x_su(i,j),y_su(i,j),str)
            hold on
        end
        su = [su [x_su(i,j);y_su(i,j)]];
    end
end
if showflag == 1
    xlabel('x/m');ylabel('y/m')
end
user = [user su];

end

