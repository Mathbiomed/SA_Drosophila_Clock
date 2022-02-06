function dx = ldmodel(t,x,AT,K,a3,b1,b2,b3,KPA,KA,Kd,KP,period,inten,ddinten,smoothL,smoothD,delay)
    % The model describing TTFL under LD
    dx=[K*transcription(x(3),AT,KPA,KA,Kd,KP)-b1*x(1);x(1)-ldlight(t,period,inten,ddinten,smoothL,smoothD,delay)*x(2);a3*x(2)-b3*x(3)];
end