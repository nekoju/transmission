# Changelog

## 0.0.2

- Removed README.rst from repo, added generation of same to setup.py.
- Added changelog.md.
- Added requirements.txt.
- Moved simulate_prior_stats to transmission/cli_tools.py.

## 0.0.3

- Added Dockerfile

## 0.0.4

- Tested prior generation with Docker and Singularity.
- Resolved a problem wherein host_Ne was actually host N.

## 0.0.4a

- Previous release had a typo in _sim().

## 0.0.4b

- Fixed math behind parameter estimation, N/Ne problem in _sim().

## 0.0.4c

- Fixed above math again.

## 0.0.5a

- Changed logit_bounds in Abc() for rho to reflect (0, 1) range.

## 0.0.5b

- Got rid of beta_nonstandard().

## 0.0.5c

- Fixed beta to use 'size' rather than 'n'.

## 0.0.6

- Added demos.
- Added num_sites, theta_w to _sim().

## 0.0.6a

- Fixed a small syntax error.

## 0.0.7

- Changed host_theta to reflect proposal.
- Added host_source to _sim() and generate_priors()
- Changed test_sim_fst() to have 0.10 rtol.
- Added tests for "nei" and "tajima" Sample.pi() methods, as well as
  Sample.num_mutants() and Sample.polymorphic()
- Updated README.md to include information on running a Jupyter notebook
  server from the Docker image.

## 0.0.8

- Tested Sample.segsites()

## 0.0.11

- Worked on fixing math and generating new priors for demo.
- Enable multiprocessing within module.

## 0.0.12

- Allow use of only a subset of simulated populations for MetaSample.
- Change package name to txmn.
- Fix to properly calculate migration rates from host_Nm with different
  migration models.
