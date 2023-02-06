functions {
}

data {
    int<lower=0> T;         // number of observations in each frame
    vector[T] x_1;            // the observations in frame 1
    vector[T] x_2;            // the observations in frame 2
    vector[T] x_3;            // the observations in frame 3

    real very_high_prior_location;
}

transformed data {
    // set the hyperparameters based on the data
    real low_location_prior_location;
    real low_location_prior_scale;

    real low_scale_prior_location;
    real low_scale_prior_scale;

    // and delta is based on the observed difference
    real delta_location;

    // vectors to hold observations of each type
    vector[1] very_high;
    vector[T-1] high;
    vector[2*T] low;

    very_high = head(x_1, 1);
    high = tail(x_1, T-1);
    low = append_row(x_2, x_3);

    // we just use the emprical values to model "low"

    // we pick the lower of the two
    low_location_prior_location = mean(x_2);
    if (mean(x_3) < low_location_prior_location) {
        low_location_prior_location = mean(x_3);
    }

    low_scale_prior_location = variance(x_2);
    if (variance(x_3) < low_scale_prior_location) {
        low_scale_prior_location = variance(x_3);
    }

    low_location_prior_scale = sqrt(low_location_prior_location);
    low_scale_prior_scale = sqrt(low_scale_prior_location);

    // the prior scales cannot be 0
    low_location_prior_scale = fmax(low_location_prior_scale, 0.1);
    low_scale_prior_scale = fmax(low_scale_prior_scale, 0.1);

    // the deltas are picked such that the location and scale reflect
    // the difference present in the empirical distributions
    delta_location = mean(high) - low_location_prior_location;
    delta_location = fmax(delta_location, 3*low_location_prior_location);

    /*
    print("very_high:", very_high);
    print("high:", high);
    print("low:", low);

    print("low_location_prior_location:", low_location_prior_location);
    print("low_scale_prior_location:", low_scale_prior_location);
    print("delta_location:", delta_location);
    */
}

parameters {
    real low_location;
    positive_ordered[1] low_scale;
}

model {
    // low does not require anything special
    low_location ~ cauchy(low_location_prior_location, low_location_prior_scale);
    low_scale ~ cauchy(low_scale_prior_location, low_scale_prior_scale);
    low ~ cauchy(low_location, low_scale[1]);

    // high and very_high are simply increases based on the low parameters
    high ~ cauchy(low_location+delta_location, low_scale[1]);
    very_high ~ cauchy(very_high_prior_location, low_scale[1]);
}
