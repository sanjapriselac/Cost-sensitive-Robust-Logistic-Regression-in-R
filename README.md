# Cost-sensitive-Robust-Logistic-Regression-in-R

Author: Sanja Priselac 

The code for implementing robust logistic regression for imbalanced data sets 
based on the Bianco-Yohai estimator. The imbalance learning problem is 
addressed by including the cost-sensitive features in the method.

The implementation is based on the iterative algorithm introduced by 
Croux and Haesbroeck (doi:10.1016/S0167-9473(03)00042-2). As a simple initial solution, 
the authors propose a weighted logistic regression based on the leverage detection 
of the MCD estimator. Since the MCD estimator cannot be computed often in practice, 
the implementation also includes the PCDist algorithm (doi:10.1016/j.ins.2012.10.017) 
for leverage detection. The weights are also used for the weighted 
Bianco-Yohai estimator. 

The code is implemented in the R programming language. 
It is a modification of the existing functions *BYlogreg*, *glmrobBY* and *glmrob* 
from the R package *robustbase*.

The performance of the implemented cost-sensitive Bianco-Yohai estimator 
in original and weighted form is analyzed in comparison to non-cost-sensitive and 
non-robust estimators and can be found in the *Evaluation* folder.

### Folder structure ###
* glmrobBy_UPDATE.R : the adoptation of the functions *BYlogreg*, *glmrobBY* and *glmrob*
* Evaluation : 
	* SimulationData.R : involves artificial data formation with 
		four configurations for outlier types (bad leverage points, 
		mislabeled observations, and the combination), two settings 
		for the number of explanatory variables (p=2 and p=10), and 
		different imbalance proportions (20%, 10%, 5%, 1% of positives). 
		The estimated parameters of the six classifiers are compared with 
		the true parameter (bias, MSE).
	* giveMeCredit.R: the evaluation of the six classifiers based on the Gini index. 
		Data can be found here: https://www.kaggle.com/c/GiveMeSomeCredit/data
* additionalFunctions.R : Additional functions for the model evaluation and simulation settings 
* Results_and_Explanatory_Analysis : 
	* MT_SimulationResults.pdf: Tables of simulation results - the bias and MSE and the mean and 
				   standard deviation of the Gini index evaluated on based on 500 simulation runs.
	* MT_GiveMeCreditDescription.pdf: Description of the explanatory variables and their univariate distribution plots (original and log-transformed).

