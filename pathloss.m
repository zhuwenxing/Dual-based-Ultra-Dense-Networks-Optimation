function pl = pathloss(dis)
pl = 140.7+37.6*log10(dis/1000);
% pl = 128.1*(dis/1000)^(-3.76);