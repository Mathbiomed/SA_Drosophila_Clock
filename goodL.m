function tf=goodL(smoothL,smoothD)
    % The function to discriminate whether the light function can be accepted
    t=0:0.01:1;
    tf=true;
    light=0.5*((tanh((t - floor(t))/smoothL)-tanh((t - floor(t)-1/2)/smoothD))+(1+tanh((t-floor(t)-1)/smoothL)));
    if (numel(findpeaks(light))~=1) || (max(light) < 0.99)
        tf=false;
        % If the light function is not unimodal, or the maximum value is
        % much less than 1, it cannot be accepted
    end
end