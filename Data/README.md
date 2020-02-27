# Data Files Included in Repo

`priors_s-0.0.1_t-1-1_r-10-10_1e6.pickle`

nchrom: 12
demes: 30
theta: 1
Nm: 2
num samples: 1000000
replicates per sample: 25
eta (mutation scaling): Normal(0, 0.1)
rho (sex ratio): Beta(10, 10)
tau (vertical transmission rate): Beta(1, 1)
prior seed: 3
random seed: 42

##columns

`fst_mean`: Mean F_ST over all sites, calculated as G_ST
`fst_sd`: F_ST standard deviation over all sites
`pi_h`: Nucleotide diversity, calculated as sum of heterozygosities across all sites.
`eta`, `tau`, `rho`: Above parameter values for each simulation.
