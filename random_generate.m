function [x, y]=random_generate(r,x0,y0)
rdis = rand(1)*r;
angle = rand(1)*2*3.1415926;
x = rdis*cos(angle)+x0;
y = rdis*sin(angle)+y0;