function [handleCDF,stats] = mycdfplot(x,linestyle)
%CDFPLOT Display an empirical cumulative distribution function.
%   CDFPLOT(X) plots an empirical cumulative distribution function (CDF) 
%   of the observations in the data sample vector X. X may be a row or 
%   column vector, and represents a random sample of observations from 
%   some underlying distribution.
%
%   H = CDFPLOT(X) plots F(x), the empirical (or sample) CDF versus the
%   observations in X. The empirical CDF, F(x), is defined as follows:
%
%   F(x) = (Number of observations <= x)/(Total number of observations)
%
%   for all values in the sample vector X. If X contains missing data
%   indicated by NaN's (IEEE arithmetic representation for
%   Not-a-Number), the missing observations will be ignored.
%
%   H is the handle of the empirical CDF curve (a Handle Graphics 'line'
%   object). 
%
%   [H,STATS] = CDFPLOT(X) also returns a statistical summary structure
%   with the following fields:
%
%      STATS.min    = minimum value of the vector X.
%      STATS.max    = maximum value of the vector X.
%      STATS.mean   = sample mean of the vector X.
%      STATS.median = sample median (50th percentile) of the vector X.
%      STATS.std    = sample standard deviation of the vector X.
%
%   In addition to qualitative visual benefits, the empirical CDF is 
%   useful for general-purpose goodness-of-fit hypothesis testing, such 
%   as the Kolmogorov-Smirnov tests in which the test statistic is the 
%   largest deviation of the empirical CDF from a hypothesized theoretical 
%   CDF.
%
%   See also QQPLOT, KSTEST, KSTEST2, LILLIETEST.

% Copyright 1993-2011 The MathWorks, Inc.


% Get sample cdf, display error message if any
[yy,xx,~,~,eid] = cdfcalc(x);
if isequal(eid,'VectorRequired')
    error(message('stats:cdfplot:VectorRequired'));
elseif isequal(eid,'NotEnoughData')
    error(message('stats:cdfplot:NotEnoughData'));
end

% Create vectors for plotting
k = length(xx);
n = reshape(repmat(1:k, 2, 1), 2*k, 1);
% xCDF    = [-Inf; xx(n); Inf];
% yCDF    = [0; 0; yy(1+n)];
xCDF    = xx(n);
yCDF    = yy(1+n);
%
% Now plot the sample (empirical) CDF staircase.
%

% pnum = 4;
% len = length(xCDF);
% slen = floor(len/pnum);  
% xxCDF = zeros(1,slen+1);
% yyCDF = zeros(1,slen+1);
% xxCDF(1) = xCDF(1);
% yyCDF(1) = yCDF(1);
% for i = 2:slen
% 	xxCDF(i) = xCDF(pnum*i);
% 	yyCDF(i) = yCDF(pnum*i);
% end
% xxCDF(slen+1) = xCDF(len);
% yyCDF(slen+1) = yCDF(len);

% c = polyfit(xCDF, yCDF, 3);
% yCDF = polyval(c, xCDF, 1);
% num = 500;
% for i = num+1:length(yCDF)-num
%     sum = 0;
%     for j = -num:num
%         sum = sum+yCDF(i+j);
%     end
%     yCDF(i)= sum/(2*num+1);
% end
xxCDF = zeros(1,length(xCDF)/2);
yyCDF = zeros(1,length(yCDF)/2);
for i = 2:2:length(xCDF)
    xxCDF(i/2) = xCDF(i);
    yyCDF(i/2) = yCDF(i);
end
% c = polyfit(xxCDF, yyCDF, 5);
% yyCDF = polyval(c, xxCDF, 1);
% for i = 1:length(xxCDF)
%     if xxCDF(i)<88
%         yyCDF(i)=0;
%     end
%     if yyCDF(i)<0
%         yyCDF(i)=0;
%     end
%     if yyCDF(i)>1
%         yyCDF(i)=1;
%     end
% end


switch linestyle
    case 1
%         hCDF = plot(xCDF , yCDF,'r');
        hold on
%         hCDF = plot(xxCDF , yyCDF, 'ro-');
        hCDF = line_fewer_markers(xxCDF , yyCDF,12, 'k-d','MarkerSize', 6);
    case 2
%         hCDF = plot(xCDF , yCDF,'b');
        hold on
%         hCDF = plot(xxCDF , yyCDF, 'b*-');
        hCDF = line_fewer_markers(xxCDF , yyCDF,12, 'r-s','MarkerSize', 6);
    case 3
%         hCDF = plot(xCDF , yCDF,'c');
        hold on
%         hCDF = plot(xxCDF , yyCDF, 'c^-');
        hCDF = line_fewer_markers(xxCDF , yyCDF,12, 'b-v','MarkerSize', 6);
    case 4
%         hCDF = plot(xCDF , yCDF,'g');
        hold on
%         hCDF = plot(xxCDF , yyCDF, 'g+-');
        hCDF = line_fewer_markers(xxCDF , yyCDF,12, '-o','MarkerSize', 6,'Color',[0.1 0.5 0.3]);
%         for i = 1:length(xxCDF)
%             if mod(i,125) == 0
%                 text(xxCDF(i)-0.8 , yyCDF(i), '*','Fontsize',18);
%             end
%         end
    case 5
%         hCDF = plot(xCDF , yCDF,'k');
        hold on
%         hCDF = plot(xxCDF , yyCDF, 'ks-');
        hCDF = line_fewer_markers(xxCDF , yyCDF,12, 'm-*','MarkerSize', 6);
    case 6
        hold on
        hCDF = line_fewer_markers(xxCDF , yyCDF,12, '->','MarkerSize', 6,'Color',[0.8 0.5 0.3]);
end


if (nargout>0), handleCDF=hCDF; end
grid  ('on')
% xlabel(getString(message('stats:cdfplot:LabelX')))
% ylabel(getString(message('stats:cdfplot:LabelFx')))
% title (getString(message('stats:cdfplot:Title')))
xlabel('Network Spectral Efficiency(b/s/Hz)' );
ylabel('CDF');
grid off

%
% Compute summary statistics if requested.
%

if nargout > 1
   stats.min    =  min(x);
   stats.max    =  max(x);
   stats.mean   =  mean(x);
   stats.median =  median(x);
   stats.std    =  std(x);
end
