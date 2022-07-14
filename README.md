# A Practical Gaussian Covariance Formulation for Applications with Anisotropic, Nonstationary, and Multivariate Fields

## Writing TODOs

- geophysical inverse problems to variational data assimilation ??
- possible abbreviations:
    - Lindgren 2011
    - Variational DA -> VDA
    - elliptic PDE, PDE, SPDE ,... etc


### Intro
- somewhere mention the Gaussian, probabilistic stuff
- remove the double statements of "we can use a high tolerance yo"
- disambiguate cov and corr ...

- get more citations of WC01 or diffusion like examples in use e.g. carrier and
  Ngodock ...? Met Office?


### Review
- Statement about "matern field" is where the covariance analytical function
  holds ... but do I actually use that term?

### Methods
- Should discuss how we get random samples, and boundary conditions

### Results

### Discussion
- In conclusions, note that time could additionally be handled as a separate
  dimension, and now "grid cells" has the interpretation of "time steps", once
  we normalize by deltat ... there are definitely other citations to mention
  here (some by the Lindgren et al crew)
- Since we have access to B and Binverse, could use this for obs operator ...



## TODO

- [x] Run some ECCO cases to compute empirical correlations, mapped to isotropic
  space
    - [x] Make `SampleDriver`
    - [x] Get ECCO-type setup to compile
    - [x] Check binaries for why this is producing NaNs
    - [x] First draft done. See sample plot of 100 sample correlation
- [ ] Compare to WC01
    - [x] Take a dive into the stability criterion
    - [x] How do correlation lengths change with number of iterations?
    - [ ] Do we still get the factor of ~3 when the diffusivity is constant
      (stationary)?
- [ ] Further Matern exploration and tightening up:
    - [ ] How inexact can we be? 10-2 looks good, checking 10-1
    - [ ] Report average numbers
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

