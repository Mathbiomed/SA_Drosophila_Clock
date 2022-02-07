# SA_Drosophila_Clock
Matlab codes for the mathematial model of the Drosophila circadian clock and its parameter esitmation using the simulated annealing (SA) method. See Jeong et al, Systematic modeling-driven experiments identify distinct molecular clockworks underlying hierarchically organized pacemaker neurons, PNAS (2022) for details.

## Code Description
1. simSA.m

> The main function for this package. If one specify the input of random number generator, then the code provides the estimated parameters of the Drosophila circadian clock models using the SA method.

Since all the below functions are automatically run by simSA.m function, users do not need to run nor edit the below functions separately.

2. nddmeasure.m

> This function measures amplitude, period, level, and relative amplitude from the simulated PER time series under DD. 

3. nldmeasure.m

> This function measures amplitude, period, relative amplitude, and peak and trough phases from the simulated PER time series under LD.

4. newddsa.m

> This function calculates all costs including the fluctuation, fitting and regularization costs with the measured amplitude, period, level, and relative amplitude using the function 'nddmeasure.m'.

5. newldsa.m

> This function calculates all costs including the entrain, fitting and regularization costs with the measured amplitude, period, relative amplitude, and peak and trough phases using the function 'nldmeasure.m'.

6. transcription.m

> This equation describes the transcription rate of the per mRNA with respect to the molar ratio between nuclear PER and CLK. 

7. ldlight.m

> This function describes the day-night change of the degradation rate of cytoplasmic PER under LD. 
As the light destabilizes the cytoplasmic PER, the function increases under L phase and decreases under D phase. To simulate such day-night change under LD, the function is construced by the sum of tanh terms.

8. goodL.m

> This function check whether tanh terms of the function 'ldlight.m' are smoothly combined.

9. model.m

> Mathematical model describing the core TTFL of the Drosophila circadian clock under DD.

10. ldmodel.m

> Mathematical model describing the core TTFL of the Drosophila circadian clock under LD.

11. isrhythmM.m

> This function check whether the model can simulate sustained per mRNA rhythms under LD.
