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
    real location_prior_location;
    real location_prior_scale;

    real scale_prior_location;
    real scale_prior_scale;

    vector[3*T] background;
    vector[2*T] temp;

    temp = append_row(x_1, x_2);
    background = append_row(temp, x_3);

    location_prior_location = mean(temp);
    scale_prior_location = variance(temp);

    location_prior_scale = sqrt(location_prior_location);
    scale_prior_scale = sqrt(scale_prior_location);

    // the prior scales cannot be 0
    location_prior_scale = fmax(location_prior_scale, 0.1);
    scale_prior_scale = fmax(scale_prior_scale, 0.1);
}

parameters {
    real location;
    positive_ordered[1] scale;
}

transformed parameters {
}

model {

    // the hyperparameters
    location ~ cauchy(location_prior_location, location_prior_scale);
    scale ~ cauchy(scale_prior_location, scale_prior_scale);

    background ~ cauchy(location, scale[1]);
}
