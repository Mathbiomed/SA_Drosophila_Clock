function EA=transcription(R,AT,KPA,KA,Kd,KP)
    % The equation describing the transcription rate of per mRNA
    EA=(KP*AT+KPA*((AT-R-Kd + sqrt((AT - R - Kd).^2 + 4*AT*Kd))/2))./((AT -KA + 2*KPA)*((AT-R-Kd+sqrt((AT-R-Kd).^2 + 4*AT*Kd))/2) + (AT + KPA).*R + (KA + KP - KPA - AT)*AT + KA*(KP + KPA));
end