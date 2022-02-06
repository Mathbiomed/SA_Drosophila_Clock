function tf=isrhythmM(AT,K,a3,b1,b2,b3,tranRatio,KPA,KA,Kd,KP,period,inten,ddinten,smoothL,smoothD,delay)
    % The function to discriminate whether the simulated per mRNA time series is rhythmic, or not
    tspan=0:0.01:period;
    msl = ldlight(tspan,period,inten,ddinten,smoothL,smoothD,0); 
    mt=find(msl==min(msl), 1, 'first');
    startL = period - 0.01*mt;
    if period ~= 0 
        sample = 15;tspan=0:period/240:(31+sample)*period;periLength=240;
        [~,proteins]= ode23(@(t,x) ldmodel(t,x,AT,K,a3,b1,b2,b3,KPA,KA,Kd,KP,period,inten,ddinten,smoothL,smoothD,delay+startL),tspan,[0 0 0]);
        proM = proteins((29+sample)*periLength:(30+sample)*periLength,1);
        proM = proM/max(proM);

        if max(proM)-min(proM) > 0.01
            tf=true;
            % If the simulated per mRNA time series fluctuates more than 1%, then it is rhythmic
        else
            tf=false;
        end
    else
        tf=false;
    end
end