function [ldcost,stabcost,ldfitcost,lddiffcost,lampcost,lpercost,lrelampcost,lphasecost,ltranTcost, propertylamp, propertylrelamp, propertylperiod, propertylpeakphase,propertylthphase]=newldsa(lnAT,lnK,lnb1,lnb2,dnAT,diffmutAT,dnK,dna3,dnb1,dnb2,dnb3,tranRatio,KPA,wtKA,wtKd,wtKP,mutKA,mutKd,mutKP,period,lnintenD,lnintenL,dnintenD,dnintenL,lnsmoothL,lnsmoothD,dnsmoothL,dnsmoothD,delay)
    % The function calculating costs with parameters under LD
    reallnwtld = [1139022.849, 242089.3524, 123922.7424, 66102.53051, 706435.8671, 983584.1296];
    reallnmutld = [660432.8507, 377446.7034,97211.61756, 70867.41019, 334259.0076,458278.6519];
    realdnwtld = [767800.0708, 887875.5591, 700344.972,1154179.688, 1908693.994, 1891455.992];
    realdnmutld = [1666175.844,847737.7061, 649666.6085, 476710.7156, 1221862.2, 1391018.067];

    realpeaklnwtld=2;realpeaklnmutld=2;realpeakdnwtld=18;realpeakdnmutld=2;
    realthlnwtld=14;realthlnmutld=14;realthdnwtld=10;realthdnmutld=14;
    normln = max([max(reallnwtld),max(reallnmutld)]);normdn = max([max(realdnwtld),max(realdnmutld)]);
    realamplnwtld=max(reallnwtld/normln)-min(reallnwtld/normln);realamplnmutld=max(reallnmutld/normln)-min(reallnmutld/normln);
    realampdnwtld=max(realdnwtld/normdn)-min(realdnwtld/normdn);realampdnmutld=max(realdnmutld/normdn)-min(realdnmutld/normdn);
    realrelamplnwtld=1-min(reallnwtld)/max(reallnwtld);realrelamplnmutld=1-min(reallnmutld)/max(reallnmutld);
    realrelampdnwtld=1-min(realdnwtld)/max(realdnwtld);realrelampdnmutld=1-min(realdnmutld)/max(realdnmutld);
    
    mutlnAT=lnAT*diffmutAT;mutdnAT=dnAT*diffmutAT;
    if period ~= 0
        % With 'nldmeasure' function, calcuate amplitude, period, phase, ... of the four models if the entrain period is not 0
        relampdnMUTld=0;leveldnMUTld=0;amplitudednMUTld=0;perioddnMUTld=0;peakphasednMUTld=0;thphasednMUTld=0;maxdnMUTld=0;
        [startL,~,amplitudeR,periodR,relampR,maxR,thphaseR,peakphaseR,costs]=nldmeasure(mutdnAT,dnK,dna3,dnb1,dnb2,dnb3,tranRatio,KPA,mutKA,mutKd,mutKP,period,dnintenL,dnintenD,dnsmoothL,dnsmoothD,delay);
        stabcost = costs;
        startLdnMUTld = startL;
        relampdnMUTld=relampR;amplitudednMUTld=amplitudeR;perioddnMUTld=periodR;peakphasednMUTld=peakphaseR;thphasednMUTld=thphaseR;maxdnMUTld=maxR;

        relamplnWTld=0;levellnWTld=0;amplitudelnWTld=0;periodlnWTld=0;peakphaselnWTld=0;thphaselnWTld=0;maxlnWTld=0;
        [startL,trsctimeld,amplitudeR,periodR,relampR,maxR,thphaseR,peakphaseR,costs]=nldmeasure(lnAT,lnK,1,lnb1,lnb2,1,tranRatio,KPA,wtKA,wtKd,wtKP,period,lnintenL,lnintenD,lnsmoothL,lnsmoothD,delay);
        stabcost = stabcost + costs;
        startLlnWTld = startL;
        trscTlnWTld=trsctimeld;relamplnWTld=relampR;amplitudelnWTld=amplitudeR;periodlnWTld=periodR;peakphaselnWTld=peakphaseR;thphaselnWTld=thphaseR;maxlnWTld=maxR;

        relampdnWTld=0;leveldnWTld=0;amplitudednWTld=0;perioddnWTld=0;peakphasednWTld=0;thphasednWTld=0;maxdnWTld=0;
        [startL,trsctimeld,amplitudeR,periodR,relampR,maxR,thphaseR,peakphaseR,costs]=nldmeasure(dnAT,dnK,dna3,dnb1,dnb2,dnb3,tranRatio,KPA,wtKA,wtKd,wtKP,period,dnintenL,dnintenD,dnsmoothL,dnsmoothD,delay);
        stabcost = stabcost + costs;
        startLdnWTld = startL;
        trscTdnWTld=trsctimeld;relampdnWTld=relampR;amplitudednWTld=amplitudeR;perioddnWTld=periodR;peakphasednWTld=peakphaseR;thphasednWTld=thphaseR;maxdnWTld=maxR;

        relamplnMUTld=0;levellnMUTld=0;amplitudelnMUTld=0;periodlnMUTld=0;peakphaselnMUTld=0;thphaselnMUTld=0;maxlnMUTld=0;
        [startL,~,amplitudeR,periodR,relampR,maxR,thphaseR,peakphaseR,costs]=nldmeasure(mutlnAT,lnK,1,lnb1,lnb2,1,tranRatio,KPA,mutKA,mutKd,mutKP,period,lnintenL,lnintenD,lnsmoothL,lnsmoothD,delay);
        stabcost = (stabcost + costs);
        startLlnMUTld = startL;
        relamplnMUTld=relampR;amplitudelnMUTld=amplitudeR;periodlnMUTld=periodR;peakphaselnMUTld=peakphaseR;thphaselnMUTld=thphaseR;maxlnMUTld=maxR;
        
        stabcot=10*stabcost;
        lampcost = 5*(min(abs(1-(amplitudednWTld/max(amplitudednMUTld, 0.0001))/(realampdnwtld/realampdnmutld)),1) +min(abs(1 - (amplitudelnWTld/max(amplitudelnMUTld,0.0001))/(realamplnwtld/realamplnmutld)),1));
        lrelampcost = 2.5*(abs(1- relampdnMUTld/realrelampdnmutld)+abs(1 - relamplnWTld/realrelamplnwtld)+abs(1 - relamplnMUTld/realrelamplnmutld)+abs(1 - relampdnWTld/realrelampdnwtld));
        lpercost = 100*(abs(1 - periodlnWTld/period) + abs(1 - periodlnMUTld/period) + abs(1 - perioddnMUTld/period) + abs(1 - perioddnWTld/period));
        lphasecost = 0.8*((1-abs(abs(realthlnwtld-thphaselnWTld)/12-1))+(1-abs(abs(realthlnmutld-thphaselnMUTld)/12-1))+(1-abs(abs(realthdnwtld-thphasednWTld)/12-1))+(1-abs(abs(realthdnmutld-thphasednMUTld)/12-1))+...
            (1-abs(abs(realpeaklnwtld-peakphaselnWTld)/12-1))+(1-abs(abs(realpeaklnmutld-peakphaselnMUTld)/12-1))+(1-abs(abs(realpeakdnwtld-peakphasednWTld)/12-1))+(1-abs(abs(realpeakdnmutld-peakphasednMUTld)/12-1)));
        lddiffcost = tanh(abs(log10(lnintenD/dnintenD))/2) +tanh(abs(log10((lnintenD+lnintenL)/(dnintenD +dnintenL)))/2)+tanh(abs(log10(lnsmoothD/dnsmoothD))/2)+tanh(abs(log10(lnsmoothL/dnsmoothL))/2);
        ldfitcost = max(min(lampcost,10),1)*lampcost + max(min(lphasecost,10),1)*lphasecost + max(min(lrelampcost,10),1)*lrelampcost + max(min(lpercost,10),1)*lpercost;
        ldcost = stabcost + ldfitcost + lddiffcost;
        % Calculate the fitting cost with the measured amplitudes, periods,and relative amplitude
    else
        relamplnWTld=0;levellnWTld=0;amplitudelnWTld=0;periodlnWTld=0;peakphaselnWTld=0;thphaselnWTld=0;maxlnWTld=0;
        relamplnMUTld=0;levellnMUTld=0;amplitudelnMUTld=0;periodlnMUTld=0;peakphaselnMUTld=0;thphaselnMUTld=0;maxlnMUTld=0;
        relampdnWTld=0;leveldnWTld=0;amplitudednWTld=0;perioddnWTld=0;peakphasednWTld=0;thphasednWTld=0;maxdnWTld=0;
        relampdnMUTld=0;leveldnMUTld=0;amplitudednMUTld=0;perioddnMUTld=0;peakphasednMUTld=0;thphasednMUTld=0;maxdnMUTld=0;
        stabcost = 40;lampcost = 10;lrelampcost = 10;lpercost = 400;lphasecost = 4.8;
        lddiffcost = 4;
        ldfitcost = max(min(lampcost,10),1)*lampcost + max(min(lphasecost,10),1)*lphasecost + max(min(lrelampcost,10),1)*lrelampcost + max(min(lpercost,10),1)*lpercost;
        ldcost = stabcost + ldfitcost + lddiffcost;
        % If entrain period is 0, then the fitting cost has its maximum value 
    end
    ltranTcost=0;
    propertylamp=[amplitudelnWTld amplitudelnMUTld amplitudednWTld amplitudednMUTld];
    propertylrelamp=[relamplnWTld relamplnMUTld relampdnWTld relampdnMUTld];
    propertylperiod=[periodlnWTld periodlnMUTld perioddnWTld perioddnMUTld];
    propertylpeakphase=[peakphaselnWTld peakphaselnMUTld peakphasednWTld peakphasednMUTld];
    propertylthphase=[thphaselnWTld thphaselnMUTld thphasednWTld thphasednMUTld];    
end