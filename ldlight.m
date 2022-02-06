function pcdeg=ldlight(t,period,inten,ddinten,smoothL,smoothD,delay)
    % The function describing the day-night change of the degradation rate of the cytoplasmic PER
    pcdeg=0.5*inten*((tanh((t - delay - period*floor((t-delay)/period))/smoothL)-tanh((t-delay-period*floor((t-delay)/period)-period/2)/smoothD))+(1+tanh((t-delay-period*floor((t-delay)/period)-period)/smoothL)))+ddinten;
end