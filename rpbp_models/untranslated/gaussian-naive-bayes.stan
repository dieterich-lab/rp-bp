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

    // vectors to hold observations of each type
    vector[3*T] background;
    vector[2*T] temp;

    temp = append_row(x_1, x_2);
    background = append_row(temp, x_3);

    // we just use the emprical values to model "background"
    background_location_prior_location = mean(background);
    background_scale_prior_location = variance(background);

    // the scale cannot be 0, so adjust if necessary
    background_scale_prior_location = fmax(background_scale_prior_location, 0.1);

    background_location_prior_scale = sqrt(background_location_prior_location);
    background_scale_prior_scale = sqrt(background_scale_prior_location);

    // and fix these scales, in case they are 0
    background_location_prior_scale = fmax(background_location_prior_scale, 0.1);
    background_scale_prior_scale = fmax(background_scale_prior_scale, 0.1);

    //print("background:", background);

    //print("background_location_prior_location:", background_location_prior_location);
    //print("background_scale_prior_location:", background_scale_prior_location);
}

parameters {
    real background_location;
    real<lower=0> background_scale;
}

model {
    background_location ~ cauchy(background_location_prior_location, background_location_prior_scale);
    background_scale ~ cauchy(background_scale_prior_location, background_scale_prior_scale);
    background ~ normal(background_location, background_scale);
}
