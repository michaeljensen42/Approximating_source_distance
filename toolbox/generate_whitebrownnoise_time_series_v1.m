function [whitenoise, noise] = generate_whitebrownnoise_time_series_v1(multiplier, trials)

    srate = 10000;
    t_tot= (trials*3) + 3;

    whitenoise = multiplier*randn(srate.*t_tot, 1);

    brownnoise = cumsum(whitenoise);

    noise = ieeg_highpass(brownnoise, srate);

end