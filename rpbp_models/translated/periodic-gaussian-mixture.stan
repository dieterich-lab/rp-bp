functions {
}

data {
    int<lower=0> T;         // number of observations in each frame
    vector[T] x_1;            // the observations in frame 1
    vector[T] x_2;            // the observations in frame 2
    vector[T] x_3;            // the observations in frame 3
}

transformed data {
    // set the hyperparameters based on the data
    real background_location_prior_location;
    real background_location_prior_scale;

    real background_scale_prior_location;
    real background_scale_prior_scale;

    real signal_location_prior_location;
    real signal_location_prior_scale;

    real signal_scale_prior_location;
    real signal_scale_prior_scale;

    // vectors to hold observations of each type
    vector[T] signal;
    vector[2*T] background;

    signal = x_1;
    background = append_row(x_2, x_3);

    // we just use the emprical values for hyperparameters
    background_location_prior_location = mean(background);
    background_scale_prior_location = variance(background);

    // the scale cannot be 0, so adjust if necessary
    background_scale_prior_location = fmax(background_scale_prior_location, 0.1);

    background_location_prior_scale = sqrt(background_location_prior_location);
    background_scale_prior_scale = sqrt(background_scale_prior_location);

    // and fix these scales, in case they are 0 (or very, very small)
    background_location_prior_scale = fmax(background_location_prior_scale, 0.1);
    background_scale_prior_scale = fmax(background_scale_prior_scale, 0.1);


    // signal hyperparameters
    signal_location_prior_location = mean(signal);
    signal_scale_prior_location = variance(signal);

    // the scale cannot be 0, so adjust if necessary
    signal_scale_prior_location = fmax(signal_scale_prior_location, 0.1);

    signal_location_prior_scale = sqrt(signal_location_prior_location);
    signal_scale_prior_scale = sqrt(signal_scale_prior_location);

    // and fix these scales, in case they are 0 (or very, very small)
    signal_location_prior_scale = fmax(signal_location_prior_scale, 0.1);
    signal_scale_prior_scale = fmax(signal_scale_prior_scale, 0.1);
}

parameters {
    real background_location;
    real<lower=0> background_scale;

    real signal_location;
    real<lower=0> signal_scale;
}

transformed parameters {
}

model {
    background_location ~ cauchy(background_location_prior_location, background_location_prior_scale);
    background_scale ~ cauchy(background_scale_prior_location, background_scale_prior_scale);
    background ~ normal(background_location, background_scale);

    signal_location ~ cauchy(signal_location_prior_location, signal_location_prior_scale);
    signal_scale ~ cauchy(signal_scale_prior_location, signal_scale_prior_scale);
    signal ~ normal(signal_location, signal_scale);

}

generated quantities {
}
