
# Segmentation models 

Segmentation models applied to overdispersed, correlated and multivariate count data. Models are Sums and Shares and Poisson log-Normal mixture models.


This project presents mixture-model-based clustering approachs to analyze multidimensional, overdispersed and correlated count data. We gather these models under the name of segmentation models. These models are applied here to clustering simulated count data.
Models are Sums and Shares (described in [Jones and Marchand (2019)](https://doi.org/10.1016/j.jmva.2018.11.011)) and Poisson log-Normal mixture models (described in [Chiquet, Mariadassou and Robin (2021)](https://doi.org/10.3389/fevo.2021.588292)). Note that we propose also a "smoothed" version of these models, their acronyms are preceded by an "s".
The following models are present in this project:

 * **Poisson sums and multinomial shares mixture model (`PoiMult` and `sPoiMult`)**. The Poisson distribution modelize the sum of the counts at each time slot, multinomial distribution is used for the joint distribution of counts.
 * **Negative binomial sums and Pólya shares mixture model (`NegPol` and `sNegPol`)**. The negative binomial distribution modelize the sum of the counts at each time slot, dirichlet-multinomial distribution is used for the joint distribution of counts.
 * **Poisson-Lognormal mixture model**. This strategy consists in using a two-layer hierarchical model, with one observation layer mod-
eling the count data (Poisson distribution) and one hidden layer that estimates the dependencies between the different count series throught a covariance matrix (Gaussian distribussion). We consider two cases: either the covariance matrix is diagonal (`PLNdiag` and `sPLNdiag`) or it is full (`PLNfull` and `sPLNfull`).


## Simulated data
Count data are collected at multiple points, at each time slot of each days. The goal is then to associate each day with the dynamics of one segment among S possible segments with a certain probability. 

The presented models can be applied to any set of multivariate count data  from which we have informations about days and time slots of collection. In this project we work with simulated data. The goal is to evaluate the capacities of the models to
correctly classify days coming from time series subject to controlled global or local regime changes. The data generation protocol is as follows: we create 5 series of counts, subject to S = 3 regime changes during a period of `nombreJ` days. Each day is subdivided into `nombreT` time slots. These series are generated from a `model_sim` model (PoiMult, NegPol or PLNdiag) with parameters specified a priori. For the PLNdiag simulation models the variance values are set to `sigma_pln`. For the NegPol simulation model the `size_NegPol` argument can be specified for the size of the negative binomial distribution.

It is possible for the user to specify a priori the segments and their characteristics. In particular the user can specify the sequence of segments over the whole period with the argument `suite_jours`. For example if this one is (1, 2, 3, 1) it means that there are three segments with different dynamics with an alternation from segment 1 to segment 2 and to segment 3 then a return of segment 1. The user can specify the size of these segments with the argument `prop_jours`. Continuing the previous example, if `prop_jours` is (.1, .6, .2, .1), then the segment sizes will be 10%, 60%, 20% and 10% respectively for the whole period.
Regime changes are generated by increasing the generated counts in the first segment and decreasing them in the third segment. For each serie this increase/decrease can be specified by the user with the argument `diff_dat`. For example if this one is set to (10, 5, 2, 10, 10), the increase/decrease will be around 10% for the first serie, 5% for the second, 2% for the third and 10% for the fourth and fifth series.


## Models usage
Once the simulated count data are generated it is possible to run one of the 8 previously specified mixture models to detect segments. The model used is specified with `model_train` argument. It is possible to run this model `Ntries` times in order to test its capacities to detect segments multiple times on a given situation. At each try, a new count dataset is generated and the model has to detect S =1 to S =  `nb_S` segments. Running multiple times the model allow to make some statistics as the mean misclassification rate of the days and the percentage of trials where the simulation found 3 segments according to the BIC criterion. Arguments are the following: 
|  Arguments   | Details  |
|  ----  | ----  |
| nb_S  | maximum number of segments  |
| nb_try | number of successfull tries for a model at a given run |
| max_try | maximum number of tries before stopping|
| paral | TRUE or FALSE for the NegPol and sNegPol models|
| nbcores | if paral is TRUE for NegPol and sNegPol models|


## Requirements
Inside the project run the following command in order to retrieve all needed packages 

```r
packrat::restore()
```

Then restart R.

## Run tests

You may specify the desired parameters for simulated data and model. All these parameters can be modified in the parameters file. Then the file train_model.r should be executed in a R environment. Results are collected in the "results" folder. "results_simu.txt" displays some elements about the simulated data, misclassification rate and percentage of well-found trials. Some visuals of segmentations are also displayed.

