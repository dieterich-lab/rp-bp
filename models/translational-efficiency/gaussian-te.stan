functions {
}

data {
    int<lower=1>                            num_ribo_replicates;
    real<lower=0>                           ribo_variance_loc;
    real<lower=0>                           ribo_variance_scale;
    vector<lower=0>[num_ribo_replicates]    ribo_abundances;

    int<lower=1>                            num_rna_replicates;
    real<lower=0>                           rna_variance_loc;
    real<lower=0>                           rna_variance_scale;
    vector<lower=0>[num_rna_replicates]     rna_abundances;
}

transformed data {
    real ribo_abundance_loc;
    real rna_abundance_loc;

    ribo_abundance_loc <- mean(ribo_abundances);
    rna_abundance_loc <- mean(rna_abundances);
}

parameters {
    real<lower=0> ribo_abundance_mean;
    real<lower=0> ribo_abundance_var;
    
    real<lower=0> rna_abundance_mean;
    real<lower=0> rna_abundance_var;
}

transformed parameters {
}

model {
    ribo_abundance_var ~ cauchy(ribo_variance_loc, ribo_variance_scale);
    ribo_abundance_mean ~ normal(ribo_abundance_loc, ribo_abundance_var);
    ribo_abundances ~ normal(ribo_abundance_mean, ribo_abundance_var);
    
    rna_abundance_var ~ cauchy(rna_variance_loc, rna_variance_scale);
    rna_abundance_mean ~ normal(rna_abundance_loc, rna_abundance_var);
    rna_abundances ~ normal(rna_abundance_mean, rna_abundance_var);
}

generated quantities {
    real log_translational_efficiency;

    log_translational_efficiency <- log(ribo_abundance_mean) - log(rna_abundance_mean);
}
