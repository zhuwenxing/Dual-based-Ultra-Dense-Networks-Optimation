function F = powerfun1( p_vec,S,M,N,g,x,w,mu,eta,gamma,Gamma,Noise )
% S = 2; %small cell number
% M = 2; %number of users in a cell
% N = 2; % channel number
p = zeros(1+S,M,N);
F = [];
%% p mapping
tt = 1;
for n = 1:N
    for k = 1:M
        for i = 1:S+1       
            p(i,k,n) = p_vec(tt);
            tt = tt+1;
        end
    end
end

%% caculate F
for i = 1:S+1
    for k = 1:M
        for n = 1:N
            %caculate f1
            tempI = 0;
            for j = 1
                if j ~= i
                    for q = 1:M
                        tempI = tempI + p(j,q,n)*g(j,(i-1)*M+k,n);
                    end
                end
            end
            tempI = tempI + Noise + Gamma*p(i,k,n)*g(i,(i-1)*M+k,n);
            f1 = Gamma*w(i,k)*x(i,k,n)*g(i,(i-1)*M+k,n)/tempI;
            %caculate f2
            f2 = 0;
            for j = 1:S+1
                if j ~= i
                    tempI = 0;
                    for c = 1:S+1
                        if c ~= j
                            for q = 1:M
                                tempI = tempI + p(c,q,n)*g(c,(i-1)*M+k,n);
                            end
                        end
                    end
                    tempI = tempI + Noise;
                    Den = tempI^2 + Gamma*p(j,k,n)*g(j,(i-1)*M+k,n)*tempI;
                    f2 = f2 + Gamma*w(j,k)*x(j,k,n)*p(j,k,n)*g(j,(i-1)*M+k,n)*g(i,(i-1)*M+k,n)/Den;
                end
            end 
            F = [F ;f1 - f2 - reallog(2)*mu(i) - reallog(2)*eta*gamma(i)];  
        end
    end
end

end

