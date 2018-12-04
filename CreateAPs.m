function ap = CreateAPs( num_sm,showflag )
%CREATEAP Summary of this function goes here
%   Detailed explanation goes here
% generate cells
% showflag = 1;
r_ma = 289; %radius of macro cell,m
r_sm = 70; %radius of small cell,m
x_sm = zeros(1,num_sm);
y_sm = zeros(1,num_sm);

ap = []; % BS coordinate
flag = 0;
if showflag == 1
    circle(r_ma,0,0);
    hold on
end
[x_sm(1),y_sm(1)] = random_generate(r_ma,0,0);
if showflag == 1
    circle(r_sm,x_sm(1),y_sm(1));
end
for i = 2:num_sm
    [x_sm(i),y_sm(i)] = random_generate(r_ma,0,0);
    flag = 0;
    while flag == 0
        flag = 1;
        for j = 1:i-1
            if distance(x_sm(i),y_sm(i),x_sm(j),y_sm(j))<r_sm*0.5 || distance(x_sm(i),y_sm(i),0,0)<r_sm
                [x_sm(i),y_sm(i)] = random_generate(r_ma,0,0);
                flag = 0;
            end
        end
    end
    if showflag == 1
        circle(r_sm,x_sm(i),y_sm(i));
    end
end
ap = [0 x_sm;0 y_sm];
end

