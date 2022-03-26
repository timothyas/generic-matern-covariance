# A Practical Gaussian Covariance Formulation for Applications with Anisotropic, Nonstationary, and Multivariate Fields

## TODO

- Performance note, see `llcreader_with_smooth.ipynb`:
    - I made a slight rewrite in `smooth_filtervar3d` to print out smoothfld on
      a per time stamp basis
    - This really cleans up the `open_smoothdataset`, now found in
      `new_smooth_store`
    - But, the classic `mds_store` really performs poorly for llc data
    - So to get initial plots, I used `read_and_store` to do this the eager way
    - But, it looks like `llcreader` can be used much more efficiently than
      `mds_store`. By default, this chunks `face` by 3 rather than 1, has
      ~4.2 tasks per chunk, `mds_store` does 5 ... but `llcreader` can instantly
      (?) chunk all z levels at once, leading to far fewer chunks. Cant do this
      with `mds_store`...
    - Before running more experiments and going through:
        1. Run experiment
        2. Read in the hacked way
        3. Create a new zarr store
      It's a good idea to see if this llcreader implementation can "just work".


- [x] Run some ECCO cases to compute empirical correlations, mapped to isotropic
  space
    - [x] Make `SampleDriver`
    - [x] Get ECCO-type setup to compile
    - [x] Check binaries for why this is producing NaNs
    - [x] First draft done. See sample plot of 100 sample correlation
- [ ] Compare to WC01
    - [ ] Take a dive into the stability criterion
    - [ ] How do correlation lengths change with number of iterations? Even if
      things are
- [ ] Further Matern exploration and tightening up:
    - [ ] How inexact can we be?
    - [ ] Higher powers
    - [ ] Intuition plots:
        - [ ] Z vs Longitude plot of samples for a single slice
        - [ ] X-Y plots near the equator
    - [ ] Correlation plot:
        - [ ] More samples
        - [ ] interquartile range (or something) rather than +/-
          std (crossing 0 issue)

- [ ] Think about an OI application


## Introduction

- Open this up to generically address the geosciences, then at the last minute
  focus on oceanography (and specifically z-coordinate model)?
- focusing on static (not dynamic) covariance formulation
- Can probably de-emphasize multivariate part, since that comes somewhat
  naturally and is nice but not the main focus of this

1. Gaussian distributions are used all the time, and a major challenge is
   specifying the covariance
2. For a number of applications in the geosciences, need anisotropic,
   nonstationary, multivariate ... and continental boundaries
3. Typically done with a WC01 approach, we make a modification and this is
   useful

## Review

(keep as is for now)

## Discretization

- Go straight from 3.1 to 3.3? If so, hyperparameter section could be ingested
  here
- Incorporate boundary conditions from the get-go?

1. Define correlation matrix C = Ainv Dz
2. Covariance comes naturally a la Weaver and Courtier
3. Specification of Phi, delta/rho, and length scales

## Description of domains used

1. Describe the PIG domain
2. Describe the ECCO/LLC domain/grid
3. Describe data used for OI problem, or other application...

## Comparison to theory

1. Show samples from both applications, give some intuition for length scales
    - 3D case: show some slices in all dimensions, maybe showing fewer length
      scales
2. Show correlation coefficient from theory and experiments
    - Might be useful to put all dimensions together, so we just have 2 panels:
      2D anisotropic case and 3D nonstationary, anisotropic case

## Application to OI problem?


## Discussion

### Comparison to WC01

- Having access to A or Ainv makes it so we can take samples or apply inverse in
  regularization ... allows for sampling based (ref) or optimization based (ref)
  inverse problems...

