functions {
}

data {
    int<lower=0> T;         // number of observations in each frame
    int<lower=0> nonzero_x_1;   // the number of 0 observations in x_1
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
    
    // we will only worry about the nonzero part of x_1 as the signal
    vector[nonzero_x_1] signal;
    vector[2*T] background;

    // we will use the overall proportion of 0s to model the missing signal
    real theta_zero_data;        // the proportion of 0s in the overall signal
    real log_p_zero_signal_data; // the log probability of observing the zeros in the signal

    // first, copy over the nonzeros to signal
    {
        int n;
        n <- 1;
        for (t in 1:T) {
            if (x_1[t] != 0) {
                signal[n] <- x_1[t];
                n <- n+1;
            }
        }
    }

    background <- append_row(x_2, x_3);

    // now, count the zeros in the entire profile
    {
        real zeros;
        int x_1_zeros;

        x_1_zeros <- T - nonzero_x_1;
        zeros <- x_1_zeros;

        for (t in 1:2*T) {
            if (background[t] == 0) {
                zeros <- zeros + 1;
            }
        }

        // we will assume the zeros in the signal follow a similar distribution
        // to the rest of the profile. also add some pseudocounts.
        theta_zero_data <- (zeros+1) / (3*T + 1);

        log_p_zero_signal_data <- binomial_log(x_1_zeros, T, theta_zero_data);
    }
    
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


    #print("signal:", signal);
    #print("background:", background);

    #print("background_location_prior_location:", background_location_prior_location);
    #print("background_scale_prior_location:", background_scale_prior_location);
    #print("background_location_prior_scale:", background_location_prior_scale);
    #print("background_scale_prior_scale:", background_scale_prior_scale);

    #print("delta_location_h:", delta_location_h);
}

parameters {
    real background_location;
    positive_ordered[1] background_scale;   // use built-in tyes to keep this positive
}

model {
    background_location ~ cauchy(background_location_prior_location, background_location_prior_scale);
    background_scale ~ cauchy(background_scale_prior_location, background_scale_prior_scale);

    #print("lp after hyperparameters", get_lp());

    // for some reason, this does not seem to work...?
    //background ~ normal(background_location, background_scale[1]);
    //print("lp after distribution: ", get_lp());
    
    {
        // but this does...?
        real p;
        p <- normal_log(background, background_location, background_scale[1]);
        increment_log_prob(p);
    }
    

    //print("lp after background: ", get_lp());

    // and take into account the 0s
    increment_log_prob(log_p_zero_signal_data);
    
    #print("lp after zero signal", get_lp());

    // signal is based on a sort of mixture model
    for (n in 1:nonzero_x_1) {
        vector[2] p;

        p[1] <- cauchy_log(signal[n], background_location + delta_location_l, background_scale[1]); // low
        p[2] <- cauchy_log(signal[n], background_location + delta_location_h, background_scale[1]); // high

        increment_log_prob(max(p));
    }

    #print("lp after data: ", get_lp());
}

generated quantities {
    real delta_l;
    real delta_h;
    real theta_zero;
    real log_p_zero_signal;

    delta_l <- delta_location_l;
    delta_h <- delta_location_h;
    theta_zero <- theta_zero_data;
    log_p_zero_signal <- log_p_zero_signal_data;
}
