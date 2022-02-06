function [startL,trsctimeld,amplitudeR,periodR,relampR,maxR,thphaseR,peakphaseR,costs]=nldmeasure(AT,K,a3,b1,b2,b3,tranRatio,KPA,KA,Kd,KP,period,inten,ddinten,smoothL,smoothD,delay)
    % The function calculating amplitude, period, phase, ... with parameters under LD
    tspan=0:0.01:period;
    msl = ldlight(tspan,period,inten,ddinten,smoothL,smoothD,0); 
    mt=find(msl==min(msl), 1, 'first');
    startL = period - 0.01*mt;
    if period ~= 0 
        sample = 15; periodR = []; ampR = []; lmm = []; tspan=0:0.01:(31+sample)*period;periLength=ceil(period/0.01);zt2=round(periLength/12);
        [times,proteins]= ode23(@(t,x) ldmodel(t,x,AT,K,a3,b1,b2,b3,KPA,KA,Kd,KP,period,inten,ddinten,smoothL,smoothD,delay+startL),tspan,[0 0 0]);
        proR = proteins(30*periLength:(30+sample)*periLength,2)+tranRatio*proteins(30*periLength:(30+sample)*periLength,3);
        expR=[];
        for i=1:sample
            expR = [expR; [proteins((29+i)*periLength+zt2,2)+tranRatio*proteins((29+i)*periLength+zt2,3),...
                proteins((29+i)*periLength+3*zt2,2)+tranRatio*proteins((29+i)*periLength+3*zt2,3),...
                proteins((29+i)*periLength+5*zt2,2)+tranRatio*proteins((29+i)*periLength+5*zt2,3),...
                proteins((29+i)*periLength+7*zt2,2)+tranRatio*proteins((29+i)*periLength+7*zt2,3),...
                proteins((29+i)*periLength+9*zt2,2)+tranRatio*proteins((29+i)*periLength+9*zt2,3),...
                proteins((29+i)*periLength+11*zt2,2)+tranRatio*proteins((29+i)*periLength+11*zt2,3)]];
        end
        proEA = transcription(proteins((30+sample)*periLength:(30+sample)*periLength,3),AT,KPA,KA,Kd,KP);
        trsctimeld = 24*length(proEA(proEA>0.1))*0.01/period;
        [maxes, maxT]= findpeaks(proR); 
        [mins, minT]= findpeaks(-proR);mins=-mins;
        % Extract peak and trough points from the simulated time series

        maxamp = max(maxes) - min(mins);
        endR = min(length(maxes), length(mins));
        maxes = maxes(1:endR); mins = mins(1:endR);maxT = maxT(1:endR); minT = minT(1:endR);
        lmm=sortrows([[maxT' minT'];[maxes' mins']]');       
        fakepeak = [];
        difference=diff(lmm(:,2));
        for x = 1:2:length(difference) 
            if abs(difference(x)) < 0.1*maxamp 
                fakepeak = [fakepeak x x + 1];
            end
        end
        lmm = removerows(lmm,'ind',fakepeak);
        % Distinguish fake peak points
        
        if ~isempty(lmm)
            if lmm(1,2) > lmm(2,2) 
                peaks = lmm(1:2:end,2); 
                ths = lmm(2:2:end,2); 
                iMaxes = 0.01*lmm(1:2:end,1); 
                iMins = 0.01*lmm(2:2:end,1);
            else
                peaks = lmm(2:2:end,2);
                ths = lmm(1:2:end,2);  
                iMaxes = 0.01*lmm(2:2:end,1);
                iMins = 0.01*lmm(1:2:end,1);
            end
        else
            iMaxes = []; 
            iMins = [];
        end
        % Delete fake peak points and construct an array of peak points
        
        if length(iMaxes) >= 10 
            ampR = abs(diff(lmm(:,2))); 
            perR = diff(iMaxes);
            % With the extracted peak and trough points, calculate amplitude and period of each cycles of oscillation
            peakValue=zeros(1,sample);thValue=zeros(1,sample);peakIndex=zeros(1,sample);thIndex=zeros(1,sample);
            for i=1:sample
                [peakV,peakI]=max(expR(i,:));
                [thV,thI]=min(expR(i,:));
                peakValue(i)=peakV;
                peakIndex(i)=peakI;
                thValue(i)=thV;
                thIndex(i)=thI;
            end

            amplitudeR = mean(peakValue-thValue); periodR = mean(perR); 
            relampR = mean((peakValue-thValue)./peakValue); 
            thphaseR = 4*mean(thIndex)-2; 
            peakphaseR = 4*mean(peakIndex)-2;
            maxR = mean(peakValue); 
            costs = (std(ampR)/mean(ampR)+std(perR)/mean(perR));
            % Calculate amplitude, period, ..., and the entrain cost 
        else
            amplitudeR=0;periodR=0;relampR=0;thphaseR=0;peakphaseR=0;maxR=proR(end);costs=2;
        end
    else
        amplitudeR=0;periodR=0;relampR=0;thphaseR=0;peakphaseR=0;maxR=0;costs=2;
        % If entrain period is 0, then the entrain cost has its maximum value 
    end
end