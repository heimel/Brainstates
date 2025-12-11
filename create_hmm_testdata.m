% idle, investigatory 1, investigatory 2
trans = [...
    0.90 0.10 0.00;
    0.00 0.90 0.10;
    0.10 0.00 0.90 ];

emissions = [...
    0.80 0.10 0.10 0.00 ; ...
    0.10 0.80 0.10 0.00 ; ...
    0.00 0.10 0.80 0.10 ];

seq = hmmgenerate(10000,trans,emissions);


trans_guess = [0.5 0.25 0.25; 0.25 0.5 0.25;0.25 0.25 0.5];
emis_guess = rand(3,4);
emis_guess = emis_guess ./ sum(emis_guess,2); % make stochastic matrix

[trans_estimate,emis_estimate] = hmmtrain(seq,trans_guess,emis_guess)


hmmdecode(seq,trans_estimate,emis_estimate)



