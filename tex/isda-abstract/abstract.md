# A Practical Formulation for an Anisotropic and Nonstationary Matern Class Correlation Operator
## Timothy A. Smith

### Cooperative Institute for Research in Environmental Sciences (CIRES) at the University of Colorado Boulder, Boulder, CO, USA
### Physical Sciences Laboratory (PSL), National Oceanic and Atmospheric Administration (NOAA), Boulder, CO, USA
### Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin, Austin, TX, USA

A key component of variational data assimilation methods is the specification of a univariate spatial correlation operator, which appears in the background-error covariance.
For oceanographic applications, complex boundaries are handled naturally with a filtering approach based on the application of an elliptic, Laplacian-like operator.
Here we show how an elliptic operator can be formulated to capture a very general Matern-type correlation structure.
We show how nonstationarity (variations in space) and anisotropy (variation in each direction) can be encoded into the operator via a simple change of variables based on user-defined normalization length scales.
The change of variables defines a mapping between the computational domain and a space where the stationary and isotropic analytical Matern correlation function is valid.

In addition to the mapping, two other hyperparameters separately control the correlation length scale (i.e.\ range) and shape.
As a practical use-case, we apply the operator to a global ocean model.
We show that when the normalizing length scales correspond to the local grid scale, the range parameter has an intuitive interpretation as the number of neighboring grid cells at which correlation drops to 0.14.

Finally, the correlation model is shown to be computationally efficient in two regards.
First, the necessary linear solve can be performed with a high tolerance (~ 0.001) while still achieving the correct statistics, requiring few iterations to converge.
Secondly, the operator's exponent, which controls the correlation shape and how many times the inverse Laplacian is applied, is linearly related to the diagonal elements of its matrix representation.
As a result, convergence properties actually improve with higher exponents.
Thus, the framework provides flexibility in controlling correlation shape.
