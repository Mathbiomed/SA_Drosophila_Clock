function dx = model(t,x,AT,K,a3,b1,b2,b3,KPA,KA,Kd,KP)
    % The model describing TTFL under DD
    dx=[K*transcription(x(3),AT,KPA,KA,Kd,KP)-b1*x(1);x(1)-b2*x(2);a3*x(2)-b3*x(3)];
end