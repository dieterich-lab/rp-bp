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

    // and delta is based on the observed difference
    real delta_location_l;
    real delta_location_h;

    // vectors to hold observations of each type
    vector[T] signal;
    vector[2*T] background;

    signal <- x_1;
    background <- append_row(x_2, x_3);
    
    // we just use the emprical values to model "background"
    background_location_prior_location <- mean(background);
    background_scale_prior_location <- variance(background);

    // the scale cannot be 0, so adjust if necessary
    background_scale_prior_location <- fmax(background_scale_prior_location, 0.1);

    background_location_prior_scale <- sqrt(background_location_prior_location);
    background_scale_prior_scale <- sqrt(background_scale_prior_location);

    // and fix these scales, in case they are 0
    background_location_prior_scale <- fmax(background_location_prior_scale, 0.1);
    background_scale_prior_scale <- fmax(background_scale_prior_scale, 0.1);

    // the deltas are picked such that the location and scale reflect
    // the difference present in the empirical distributions
    delta_location_l <- mean(signal) - background_location_prior_location;
    delta_location_l <- fmax(delta_location_l, 1);

    delta_location_h <- fmax(2*delta_location_l, 10);


    //print("signal:", signal);
    //print("background:", background);

    //print("background_location_prior_location:", background_location_prior_location);
    //print("background_scale_prior_location:", background_scale_prior_location);

    //print("delta_location_h:", delta_location_h);
}

parameters {
    real background_location;
    positive_ordered[1] background_scale;   // use built-in tyes to keep this positive
}

model {
    background_location ~ cauchy(background_location_prior_location, background_location_prior_scale);
    background_scale ~ cauchy(background_scale_prior_location, background_scale_prior_scale);
    background ~ normal(background_location, background_scale[1]);

    // signal is based on a sort of mixture model
    for (n in 1:T) {
        vector[2] p;

        p[1] <- cauchy_log(signal[n], background_location + delta_location_l, background_scale[1]);
        p[2] <- cauchy_log(signal[n], background_location + delta_location_h, background_scale[1]);

        increment_log_prob(max(p));
    }
}

generated quantities {
    real delta_l;
    real delta_h;

    delta_l <- delta_location_l;
    delta_h <- delta_location_h;
}
