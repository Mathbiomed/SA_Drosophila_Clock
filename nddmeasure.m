function [trsctimedd,amplitudeR,periodR,relampR,levelR,maxR,costf]=nddmeasure(AT,K,a3,b1,b2,b3,tranRatio,KPA,KA,Kd,KP)
    % The function calculating amplitude, period, level, ... with parameters under DD
    syms r
    steadyR=double(AT*vpasolve(transcription(r*AT,AT,KPA,KA,Kd,KP) == r*AT*b1*b2*b3/(K*a3),r));
    steadyRT=steadyR*(tranRatio+b3/a3);
    initM=0;initPc=0;initR=0;
    osc=1;
    tspan = 0:0.1:60;
    [times,proteins]= ode23(@(t,x) model(t,x,AT,K,a3,b1,b2,b3,KPA,KA,Kd,KP),tspan,[initM initPc initR]);
    initM=proteins(end,1);initPc=proteins(end,2);initR=proteins(end,3); 
    proR = proteins(:,2)+tranRatio*proteins(:,3); 
    proEA = transcription(proteins(:,3),AT,KPA,KA,Kd,KP);tt=times;
    [peakR, peakT] = findpeaks(proR,'MinPeakHeight',1.1*steadyRT); 
    [thR, thT] = findpeaks(-proR,'MinPeakHeight',-0.9*steadyRT);thR=-thR;    
    endR = min(length(peakR), length(thR)); 
    peakR = peakR(1:endR);thR = thR(1:endR);
    peakT = peakT(1:endR);thT = thT(1:endR);
    % Extract peak and trough points from the simulated time series
    
    iter=1;
    while endR <= 5 && iter <= 20
        iter=iter+1;
        if max(abs(proR(end-100:end) - steadyRT)) < 0.01
            osc=0;trsctimedd=0;amplitudeR=0;periodR=0;relampR=0;levelR=steadyRT;maxR=levelR;costf=50;
            break
        % If the maximum point stay near the steady state, stop to solve ODE and then define amplitude, peiod, and ... as 0
        else
            [times,proteins]= ode23(@(t,x) model(t,x,AT,K,a3,b1,b2,b3,KPA,KA,Kd,KP),tspan,[initM initPc initR]);
            initM=proteins(end,1);initPc=proteins(end,2);initR=proteins(end,3);
            proR = [proR; proteins(:,2)+tranRatio*proteins(:,3)];
            proEA = [proEA; transcription(proteins(:,3),AT,KPA,KA,Kd,KP)];
            tt=[tt; times+tt(end)]; 
            [peakR, peakT] = findpeaks(proR,'MinPeakHeight',1.1*steadyRT); 
            [thR, thT] = findpeaks(-proR,'MinPeakHeight',-0.9*steadyRT);thR=-thR; 
            endR = min(length(peakR), length(thR)); 
            peakR=peakR(1:endR);thR=thR(1:endR);
            peakT=peakT(1:endR);thT=thT(1:endR);
        end
    end
    if endR <= 5
        osc=0;
        % If the number of the detected peak is less than 5, there is no rhythm
    end
    if osc == 1
        riffle=zeros(1,length(peakR)+length(thR));
        riffle(1:2:end)=peakR;
        riffle(2:2:end)=thR;
        ampR = abs(diff(riffle));
        % Calculate differences of peaks and trough to construct the array of amplitudes for each cycles
        perR = peakT;
        decayR = ampR(2:end)./ampR(1:end-1);
        % Calculate decaying rate of rhythms by comparing the amplitude of each cycles
        
        if decayR(end) < 0.95 || decayR(end-1) < 0.95 || isempty(perR)
            amplitudeR=0;periodR=0;relampR=0;trsctimedd=0;maxR = peakR(end);levelR = steadyRT; 
            costf = max(min((2-decayR(end-1)-decayR(end)),2),0);
        else
            amplitudeR = ampR(end); 
            periodR = 0.1*(perR(end) - perR(end-1)); 
            maxR = peakR(end);levelR = maxR - amplitudeR/2; 
            relampR = amplitudeR/maxR; 
            oneEA = proEA(end-floor(periodR/0.1):end); 
            trsctimedd = 24*length(oneEA(oneEA>0.1))*0.1/periodR;
            costf = max(min((2-decayR(end-1)-decayR(end)),2),0);
            % Calculate amplitude, period, ..., and the fluctuation cost
        end
    else
        amplitudeR=0;periodR=0;relampR=0;trsctimedd=0;maxR = steadyRT;levelR = steadyRT; 
        costf = 2;
        % If there is no rhythm, then amplitude, period, ... are 0, and the fluctuation cost has its maximum value
    end
%     disp('DD profile')
%     disp(proR(end-15:end)')
end