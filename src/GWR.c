#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Rdynload.h"
/*
Copyright (C) 2018 Akio Onogi

Released under the MIT license
http://opensource.org/licenses/mit-license.php
*/

/*
digamma function from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.C. But now cannot be found
*/
double Digamma(const double x)
{
	double c = 8.5;
	double d1 = -0.5772156649;
	double r;
	double s = 0.00001;
	double s3 = 0.08333333333;
	double s4 = 0.0083333333333;
	double s5 = 0.003968253968;
	double value;
	double y;

	y = x;
	value = 0.0;
	if (y <= s)
	{
		value = d1 - 1.0 / y;
		return (value);
	}
	while (y < c)
	{
		value = value - 1.0 / y;
		y = y + 1.0;
	}
	r = 1.0 / y;
	value = value + log(y) - 0.5 * r;
	r = r * r;
	value = value - r * (s3 - r * (s4 - r * s5));

	return (value);
}

/* C = t(A)%*%B
A: m x n
B: m x k
C: n x k
The order of elements in each matrix is [row1column1, row2column1, ..., row1column2, row2column2, ...]
*/
void Innerproduct_tAB(double *A, double *B, int n, int m, int k, double *C)
{
	int nn, mm, kk, pos1;

	for (nn = 0; nn < n; nn++)
	{
		for (kk = 0; kk < k; kk++)
		{
			pos1 = kk * n + nn;
			C[pos1] = 0.0;
			for (mm = 0; mm < m; mm++)
			{
				C[pos1] += A[nn*m + mm] * B[kk*m + mm];
			}
		}
	}
}

/* C = t(A)%*%D%*%A where D is diagonal
A: m x n
D: m x m
C: n x n
The order of elements in each matrix is [row1column1, row2column1, ..., row1column2, row2column2, ...]
*/
void Innerproduct_tADA(double *A, double *D, int n, int m, double *C)
{
	int nn1, nn2, mm, pos1, pos2;

	for (nn1 = 0; nn1 < n; nn1++)
	{
		for (nn2 = nn1; nn2 < n; nn2++)
		{
			pos1 = nn2 * n + nn1;
			pos2 = nn1 * n + nn2;
			C[pos1] = 0.0;
			for (mm = 0; mm < m; mm++)
			{
				C[pos1] += A[nn1 * m + mm] * A[nn2 * m + mm] * D[mm];
			}
			C[pos2] = C[pos1];
		}
	}
}

/* Genome-wide regression based on variational Bayesian algorithms

This script was originally developed in VIGoR (Onogi and Iwata 2016) and modified in a stand alone C application Heading (Onogi et al. 2016)

The modifications from Heading are

(1) function type is changed from int to void.
(2) All arguments are now pointers.
(3) Objects stored in Xstruct, Ystruct, and Hstruct are now stored in each object independently.
(4) Genomic BLUP (GBLUP) was added as Priortype 3 and Methodcode 8. The inverse of the covariance matrix is defined by iA.
(5) X_use (X to Y) and Y_use (Y to X) are used to link Y and X (genotype), that is, design matries. But note that Both cannot be used to link X with multiple Y.
(6) EM algorithms were removed
(7) Quit initializing objects (arguments)
(8) Polygenic effects can be added (except for GBLUP).

Missing values in Y should be eliminated before this function.
Also Y should be standardized before this function

The inverse of Tau0 is the residual variance.
This is used as the variance of the prior distribution when the model parameters are inferred.

*/
static void GenomeWideRegression(int *Priortype, int *Methodcode, int *CondResidual, int *P, int *F, int *Ny, int *Nx, int *Niterations,
	double *Y_stobs, int *YtoX, double *Y_variance, double *Y_expErrors,
	double *X_covariate, double *X_x2, double *X_expEffect, double *X_varEffect, double *X_exp2Effect, double *X_expGamma, double *X_exp2Gamma,
	double *X_expTau2, double *X_expInTau2, double *X_expEta2, double *X_expSigma2, double *X_s2,
	double *Q_covariate, double *Q_expEffect, double *Q_varEffect, double *Q_exp2Effect, double *Q_x2,
	double *H_deltaShape, double *H_deltaRate, double *H_etaShape, double *H_etaRate, double *H_v, double *H_s2, double *H_pi, double *H_c,
	double *expDelta2, double *Tau0, int *Order,
	double *Eval, double *Evec, double *tEvec, double *U_expBv, double *U_varBv, int *XtoY, double *U_expSigma2, double *U_varSigma2, double *Threshold, int *Polygenic)
{

	/* Count iterations */
	int		ite;

	/* For repeat statement */
	int		record, locus, target;

	/* Temporary objects */
	double	temp = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0;

	/* New values */
	double	prop, prop2;

	/* for update of Re2 */
	double	sumVarB;

	/* Product of hyperparameters (when Priortype is 2) */
	double	vS2 = 0.0;

	/* For variable selection (wBSR, SSVS, MIX, and BayesC) */
	double	logPi = 0.0, log1minusPi = 0.0;

	/* Probability of gamma (for wBSR, SSVS, MIX, and BayesC) */
	double	ProbInclude = 0.0, ProbExclude = 0.0, ProbIncludeConstant = 0.0, ProbExcludeConstant = 0.0;

	/* Shape (a) and rate (b) parameters of posterior gamma distributions (for BL and EBL) */
	double	a1 = 0.0, b1 = 0.0, a2 = 0.0, b2 = 0.0;

	/* for update of Delta2 */
	double	sumEta2InTau2 = 0.0;

	/* Used when CondResidual == 1 */
	double	sumTau2B2 = 0.0;

	/* Sum of marker effects. Used when BayesC, SSVS, and MIX */
	double	sumGammaB2[2] = { 0.0, 0.0 };

	/* Used when MIX and BayesC (sumGamma, sum of gamma), and when Mix (vsS2, product of hyperparameters). */
	double	sumGamma[2] = { 0.0, 0.0 }, vcS2 = 0.0;

	/* Used when SSVS */
	double	invC = 0.0, OneMinusInvC = 0.0, logC = 0.0;

	/* check convergence */
	double	Check1, Check2;

	/* Limits of estimates */
	double	Lowesteffect = 1.0e-150;

	/* Object for GBLUP */
	int		level, level2;
	double	S2 = 0.0, vN = 0.0, *prop3 = NULL, *Yr = NULL, *D = NULL, vS2_poly = 0.0;

	switch (Priortype[0]) {
	case 2: /* wBSR, BayesB, BayesC, SSVS, and MIX */
		vS2 = H_v[0] * H_s2[0];
		switch (Methodcode[0]) {
		case 3: /* wBSR */
			if (H_pi[0]<1.0) { logPi = log(H_pi[0]); log1minusPi = log(1.0 - H_pi[0]); } /* when BayesB  */
			break;
		case 4: /* BayesC */
			if (H_pi[0]<1.0) /* BayesC */
			{
				logPi = log(H_pi[0]); log1minusPi = log(1.0 - H_pi[0]);
				sumGamma[0] = 0.0; sumGamma[1] = 0.0;
				sumGammaB2[0] = 0.0; sumGammaB2[1] = 0.0;
				for (locus = 0; locus < P[0]; locus++)
				{
					sumGamma[0] += X_expGamma[locus];
					sumGammaB2[0] += X_expGamma[locus] * X_exp2Effect[locus];
				}
			}
			else
			{	/* BRR */
				sumGammaB2[0] = 0.0; sumGammaB2[1] = 0.0;
			}
			break;
		case 5: /* SSVS */
			logPi = log(H_pi[0]); log1minusPi = log(1.0 - H_pi[0]);
			invC = 1.0 / H_c[0]; OneMinusInvC = 1.0 - invC; logC = log(H_c[0]);
			break;
		case 6: /* MIX */
			logPi = log(H_pi[0]); log1minusPi = log(1.0 - H_pi[0]);
			sumGamma[0] = 0.0; sumGamma[1] = 0.0;
			vcS2 = H_c[0] * H_v[0] * H_s2[0];
			for (locus = 0; locus < P[0]; locus++)
			{
				sumGamma[0] += X_expGamma[locus];
			}
			break;
		case 7: /* BayesB */
			if (H_pi[0]<1.0) { logPi = log(H_pi[0]); log1minusPi = log(1.0 - H_pi[0]); } /* when BayesB  */
			break;
		}
		break;
	case 3:/* GBLUP */
		vS2_poly = H_v[0] * H_s2[0];
		vN = H_v[0] + (double)Nx[0];
		S2 = 1.0;
		prop3 = (double*)calloc(Nx[0], sizeof(double));
		Yr = (double*)calloc(Nx[0], sizeof(double));
		D = (double*)calloc(Nx[0], sizeof(double));
		break;
	}

	if (Polygenic[0])
	{
		vS2_poly = 0.0;
		vN = -2.0 + (double)Nx[0];
		S2 = 1.0;
		prop3 = (double*)calloc(Nx[0], sizeof(double));
		Yr = (double*)calloc(Nx[0], sizeof(double));
		D = (double*)calloc(Nx[0], sizeof(double));
	}

	/* Calcualte residual errors */
	switch (Priortype[0])
	{
	default:
		for (record = 0; record<Ny[0]; record++)
		{
			Y_expErrors[record] = Y_stobs[record];
			for (locus = 0; locus<F[0]; locus++)
				Y_expErrors[record] -= Q_covariate[locus*Ny[0] + record] * Q_expEffect[locus];

			if (Methodcode[0] == 3 || Methodcode[0] == 4 || Methodcode[0] == 7)
			{
				for (locus = 0; locus<P[0]; locus++)
					Y_expErrors[record] -= X_covariate[locus*Nx[0] + YtoX[record]] * X_expEffect[locus] * X_expGamma[locus];
			}
			else
			{
				for (locus = 0; locus<P[0]; locus++)
					Y_expErrors[record] -= X_covariate[locus*Nx[0] + YtoX[record]] * X_expEffect[locus];
			}
		}
		if (Polygenic[0])
		{
			for (record = 0; record<Ny[0]; record++)
			{
				Y_expErrors[record] -= U_expBv[YtoX[record]];
			}
		}
		break;
	case 3:/*GBLUP*/
		for (record = 0; record<Ny[0]; record++)
		{
			Y_expErrors[record] = Y_stobs[record];
			for (locus = 0; locus<F[0]; locus++)
				Y_expErrors[record] -= Q_covariate[locus*Ny[0] + record] * Q_expEffect[locus];

			Y_expErrors[record] -= U_expBv[YtoX[record]];
		}
		break;
	}

	/* start optimization */
	for (ite = 1; ite <= Niterations[0]; ite++)
	{
		/* To check convergence */
		Check1 = 0.0; Check2 = 0.0;

		/* For update of residual variance. Only used when variational Bayesian analysis */
		sumVarB = 0.0;

		/* update of fixed effects */
		for (target = 0; target<F[0]; target++)
		{
			for (record = 0, temp = 0.0; record<Ny[0]; record++)
				temp += Q_covariate[target*Ny[0] + record] * (Y_expErrors[record] + Q_covariate[target*Ny[0] + record] * Q_expEffect[target]);

			temp *= Tau0[0];
			temp2 = 1.0 / (Q_x2[target] * Tau0[0]);

			prop = temp * temp2;
			prop2 = prop * prop + temp2;
			for (record = 0; record<Ny[0]; record++)
				Y_expErrors[record] += (Q_covariate[target*Ny[0] + record] * (Q_expEffect[target] - prop));

			sumVarB += (Q_x2[target] * temp2);

			Check1 += pow((prop - Q_expEffect[target]), 2.0);
			Check2 += pow(prop, 2.0);
			Q_expEffect[target] = prop;
			Q_exp2Effect[target] = prop2;
			Q_varEffect[target] = temp2;
		}

		/* update of B */
		switch (Priortype[0]) {
		case 1: /* BL and EBL */
			for (locus = 0; locus<P[0]; locus++)
			{
				target = Order[locus];

				for (record = 0, temp = 0.0; record<Ny[0]; record++)
					temp += X_covariate[target*Nx[0] + YtoX[record]] * (Y_expErrors[record] + X_covariate[target*Nx[0] + YtoX[record]] * X_expEffect[target]);

				temp *= Tau0[0];
				if (CondResidual[0]) { temp3 = Tau0[0]; }
				else { temp3 = 1.0; }
				temp2 = 1.0 / (X_x2[target] * Tau0[0] + X_expTau2[target] * temp3);

				/* prop: E[B], prop2: E[B^2], temp2: V[B] */
				prop = temp * temp2;	if (fabs(prop)<Lowesteffect) { prop = Lowesteffect; }
				prop2 = prop * prop + temp2;

				for (record = 0; record<Ny[0]; record++)
					Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * (X_expEffect[target] - prop));
				sumVarB += (X_x2[target] * temp2);

				Check1 += pow((prop - X_expEffect[target]), 2.0);
				Check2 += pow(prop, 2.0);

				X_expEffect[target] = prop;
				X_exp2Effect[target] = prop2;
				X_varEffect[target] = temp2;
			}
			break;

		case 2: /* wBSR, BRR, BayesC, SSVS, and MIX */
			switch (Methodcode[0]) {
			case 4: /* BayesC */
				if (H_pi[0]<1.0)
				{	/* BayesC */
					sumGammaB2[1] = 0.0; sumGamma[1] = 0.0;
					//ProbIncludeConstant = Digamma(0.5*(H_v[0] + sumGamma[0])) - 0.5*log(0.5*(sumGammaB2[0] + vS2)) + logPi - log1minusPi;
					ProbIncludeConstant = 0.5 * Digamma(0.5*(H_v[0] + sumGamma[0])) - 0.5*log(0.5*(sumGammaB2[0] + vS2)) + logPi;
				}
				else
				{	/* BRR */
					sumGammaB2[0] = 0.0;
				}
				break;
			case 5: /* SSVS */
				sumGammaB2[0] = 0.0; sumGammaB2[1] = 0.0;
				break;
			case 6: /* MIX */
				sumGammaB2[0] = 0.0; sumGammaB2[1] = 0.0;
				break;
			}
			for (locus = 0; locus<P[0]; locus++)
			{
				target = Order[locus];
				if (Methodcode[0] == 3 || Methodcode[0] == 4 || Methodcode[0] == 7) /* wBSR, BayesB, and BayesC */
				{
					for (record = 0, temp = 0.0; record<Ny[0]; record++)
						temp += X_covariate[target*Nx[0] + YtoX[record]] * (Y_expErrors[record] + X_covariate[target*Nx[0] + YtoX[record]] * X_expEffect[target] * X_expGamma[target]);
					if (Methodcode[0] == 3) temp *= X_expGamma[target];
				}
				else
				{	/* SSVS, and MIX */
					for (record = 0, temp = 0.0; record<Ny[0]; record++)
						temp += X_covariate[target*Nx[0] + YtoX[record]] * (Y_expErrors[record] + X_covariate[target*Nx[0] + YtoX[record]] * X_expEffect[target]);
				}
				temp *= Tau0[0];

				switch (Methodcode[0]) {
				case 3:	/* wBSR */
					temp2 = 1.0 / (X_x2[target] * Tau0[0] * X_exp2Gamma[target] + 1.0 / X_s2[target]);
					break;
				case 4: /* BayesC */
					temp2 = 1.0 / (X_x2[target] * Tau0[0] + 1.0 / X_s2[0]);
					break;
				case 5: /* SSVS */
					temp2 = 1.0 / (X_x2[target] * Tau0[0] + (X_expGamma[target] * OneMinusInvC + invC) / X_s2[0]);
					break;
				case 6: /* MIX */
					temp2 = 1.0 / (X_x2[target] * Tau0[0] + X_expGamma[target] / X_s2[0] + (1.0 - X_expGamma[target]) / X_s2[1]);
					break;
				case 7: /* BayesB */
					temp2 = 1.0 / (X_x2[target] * Tau0[0] + 1.0 / X_s2[target]);
					break;
				}

				/* prop: E[B], prop2: E[B^2], temp2: V[B] */
				prop = temp * temp2;
				prop2 = pow(prop, 2.0) + temp2;

				switch (Methodcode[0]) {
				case 3:	/* wBSR */
					for (record = 0; record<Ny[0]; record++)
						Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * X_expGamma[target] * (X_expEffect[target] - prop));
					if ((int)H_pi[0] == 1) /* when BayesA */
						sumVarB += (X_x2[target] * temp2);
					break;
				case 4:	/* BayesC */
					if (H_pi[0]<1.0)
					{	/* BayesC */
						/* update Gamma */
						ProbInclude = 0.5*temp2*temp*temp + 0.5*log(temp2);
						ProbInclude += ProbIncludeConstant;
						if (ProbInclude>20.0) ProbInclude = 20.0; /* to avoid overflow */
						ProbInclude = exp(ProbInclude);
						//ProbExclude = 1.0;
						//temp3 = ProbInclude / (ProbInclude + ProbExclude);
						temp3 = ProbInclude / (ProbInclude + 1.0 - H_pi[0]);

						/* update residuals */
						for (record = 0; record<Ny[0]; record++)
						{
							Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * X_expGamma[target] * X_expEffect[target]);
							Y_expErrors[record] -= (X_covariate[target*Nx[0] + YtoX[record]] * temp3 * prop);
						}
						X_expGamma[target] = temp3;
						X_exp2Gamma[target] = pow(X_expGamma[target], 2.0) + X_expGamma[target] * (1.0 - X_expGamma[target]);
						sumVarB += (X_x2[target] * X_expGamma[target] * (prop2 - X_expGamma[target] * prop*prop));
						sumGammaB2[1] += prop2 * X_expGamma[target];
						sumGamma[1] += X_expGamma[target];
						/*---- Note -----------------------------------------------------------------------------------------------------
						sumGamma[0] and sumGammaB2[0] are not updated here, because these are the posterior parameters of X_expSigma[0].
						These values are replaced by sumGamma[1] and sumGammaB2[1] at the update of X_expSigma[0].
						---------------------------------------------------------------------------------------------------------------*/
					}
					else
					{	/* BRR */
						for (record = 0; record<Ny[0]; record++)
							Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * (X_expEffect[target] - prop));
						sumVarB += (X_x2[target] * temp2);
						sumGammaB2[0] += prop2;
					}
					break;
				case 5: /* SSVS */
					for (record = 0; record<Ny[0]; record++)
						Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * (X_expEffect[target] - prop));
					sumVarB += (X_x2[target] * temp2);
					sumGammaB2[0] += prop2 * X_expGamma[target];
					sumGammaB2[1] += prop2 * (1.0 - X_expGamma[target]);
					break;
				case 6:	/* MIX (same as SSVS) */
					for (record = 0; record<Ny[0]; record++)
						Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * (X_expEffect[target] - prop));
					sumVarB += (X_x2[target] * temp2);
					sumGammaB2[0] += prop2 * X_expGamma[target];
					sumGammaB2[1] += prop2 * (1.0 - X_expGamma[target]);
					break;
				case 7:	/* BayesB */
					if (H_pi[0]<1.0)
					{	/* BayesB */
						/* update Gamma */
						ProbIncludeConstant = 0.5 * Digamma(0.5*(H_v[0] + X_expGamma[target])) - 0.5*log(0.5*(X_expGamma[target] * X_exp2Effect[target] + vS2)) + logPi;
						ProbInclude = 0.5*temp2*temp*temp + 0.5*log(temp2);
						ProbInclude += ProbIncludeConstant;
						if (ProbInclude>20.0) ProbInclude = 20.0; /* to avoid overflow */
						ProbInclude = exp(ProbInclude);
						temp3 = ProbInclude / (ProbInclude + 1.0 - H_pi[0]);

						/* update residuals */
						for (record = 0; record<Ny[0]; record++)
						{
							Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * X_expGamma[target] * X_expEffect[target]);
							Y_expErrors[record] -= (X_covariate[target*Nx[0] + YtoX[record]] * temp3 * prop);
						}
						X_expGamma[target] = temp3;
						X_exp2Gamma[target] = pow(X_expGamma[target], 2.0) + X_expGamma[target] * (1.0 - X_expGamma[target]);
						sumVarB += (X_x2[target] * X_expGamma[target] * (prop2 - X_expGamma[target] * prop*prop));
					}
					else
					{	/* BayesA */
						for (record = 0; record<Ny[0]; record++)
							Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * (X_expEffect[target] - prop));
						sumVarB += (X_x2[target] * temp2);
					}
					break;
				}
				Check1 += pow((prop - X_expEffect[target]), 2.0);
				Check2 += pow(prop, 2.0);
				X_expEffect[target] = prop;
				X_exp2Effect[target] = prop2;
				X_varEffect[target] = temp2;
			}
			break;
		}

		/* update of Tau2 or Sigma2*/
		switch (Priortype[0]) {
		case 1: /* BL and EBL */
			for (locus = 0, sumEta2InTau2 = 0.0, sumTau2B2 = 0.0; locus<P[0]; locus++)
			{
				target = Order[locus];

				if (CondResidual[0]) { temp = Tau0[0]; }
				else { temp = 1.0; }
				prop = sqrt(expDelta2[0] * X_expEta2[target] / (X_exp2Effect[target] * temp));
				Check1 += pow((prop - X_expTau2[target]), 2.0);
				Check2 += pow(prop, 2.0);

				if (CondResidual[0])
					sumTau2B2 += X_exp2Effect[target] * prop;

				X_expTau2[target] = prop;

				prop = 1.0 / X_expTau2[target] + 1.0 / (expDelta2[0] * X_expEta2[target]);
				sumEta2InTau2 += prop * X_expEta2[target];

				X_expInTau2[target] = prop;
			}
			break;

		case 2:	/* wBSR, BRR, SSVS, MIX, and BayesC */
			switch (Methodcode[0]) {
			case 3:	/* wBSR */
				for (locus = 0; locus<P[0]; locus++)
				{
					target = Order[locus];
					temp = X_exp2Effect[target] + vS2;
					X_expSigma2[target] = temp / (H_v[0] - 1.0);
					prop = temp / (H_v[0] + 1.0);

					Check1 += pow((prop - X_s2[target]), 2.0);
					Check2 += pow(prop, 2.0);

					X_s2[target] = prop;
				}
				break;
			case 4: /* BayesC */
				if (H_pi[0]<1.0)
				{	/* BayesC */
					sumGammaB2[0] = sumGammaB2[1];
					sumGamma[0] = sumGamma[1];
					temp = sumGammaB2[0] + vS2;
					X_expSigma2[0] = temp / (H_v[0] + sumGamma[0] - 2.0);
					prop = temp / (H_v[0] + sumGamma[0]);
				}
				else
				{	/* BRR */
					temp = sumGammaB2[0] + vS2;
					X_expSigma2[0] = temp / (H_v[0] + (double)P[0] - 2.0);
					prop = temp / (H_v[0] + (double)P[0]);
				}
				Check1 += pow((prop - X_s2[0]), 2.0);
				Check2 += pow(prop, 2.0);
				X_s2[0] = prop;
				break;
			case 5: /* SSVS */
				temp = sumGammaB2[0] + sumGammaB2[1] * invC + vS2;
				X_expSigma2[0] = temp / (H_v[0] + (double)P[0] - 2.0);
				prop = temp / (H_v[0] + (double)P[0]);

				Check1 += pow((prop - X_s2[0]), 2.0);
				Check2 += pow(prop, 2.0);
				X_s2[0] = prop;
				break;
			case 6: /* MIX */;
				temp = sumGammaB2[0] + vS2;
				X_expSigma2[0] = temp / (H_v[0] + sumGamma[0] - 2.0);
				prop = temp / (H_v[0] + sumGamma[0]);

				Check1 += pow((prop - X_s2[0]), 2.0);
				Check2 += pow(prop, 2.0);
				X_s2[0] = prop;

				temp = sumGammaB2[1] + vcS2;
				X_expSigma2[1] = temp / (H_v[0] + (double)P[0] - sumGamma[0] - 2.0);
				prop = temp / ((double)P[0] - sumGamma[0]);

				Check1 += pow((prop - X_s2[1]), 2.0);
				Check2 += pow(prop, 2.0);
				X_s2[1] = prop;
				break;
			case 7: /* BayesB */
				if (H_pi[0]<1.0)
				{	/* BayesB */
					for (locus = 0; locus<P[0]; locus++)
					{
						target = Order[locus];
						temp = X_expGamma[target] * X_exp2Effect[target] + vS2;
						X_expSigma2[target] = temp / (H_v[0] + X_expGamma[target] - 2.0);
						prop = temp / (H_v[0] + X_expGamma[target]);
						Check1 += pow((prop - X_s2[target]), 2.0);
						Check2 += pow(prop, 2.0);
						X_s2[target] = prop;
					}
				}
				else
				{	/* BayesA */
					for (locus = 0; locus<P[0]; locus++)
					{
						target = Order[locus];
						temp = X_exp2Effect[target] + vS2;
						X_expSigma2[target] = temp / (H_v[0] - 1.0);
						prop = temp / (H_v[0] + 1.0);
						Check1 += pow((prop - X_s2[target]), 2.0);
						Check2 += pow(prop, 2.0);
						X_s2[target] = prop;
					}
				}
				break;
			}
			break;
		}

		/* update of Delta2 (for BL and EBL)*/
		if (Priortype[0] == 1)
		{
			a2 = (double)P[0] + H_deltaShape[0];
			b2 = 0.5 * sumEta2InTau2 + H_deltaRate[0];
			prop = a2 / b2;

			Check1 += pow((prop - expDelta2[0]), 2.0);
			Check2 += pow(prop, 2.0);
			expDelta2[0] = prop;
		}

		/* update of Eta2 (for EBL)*/
		if (Methodcode[0] == 2)
		{
			for (locus = 0; locus<P[0]; locus++)
			{
				target = Order[locus];
				prop = (1.0 + H_etaShape[0]) / (0.5 * expDelta2[0] * X_expInTau2[target] + H_etaRate[0]);

				Check1 += pow((prop - X_expEta2[target]), 2.0);
				Check2 += pow(prop, 2.0);
				X_expEta2[target] = prop;
			}
		}

		/* Update of Gamma (for wBSR, SSVS, and MIX) */
		switch (Methodcode[0]) {
		case 3: /* wBSR */
			if (H_pi[0]<1.0)
			{ /*when BayesB */
				for (locus = 0; locus<P[0]; locus++)
				{
					target = Order[locus];

					for (record = 0, ProbInclude = 0.0, ProbExclude = 0.0; record<Ny[0]; record++)
					{
						ProbInclude += pow((Y_expErrors[record] + (X_expGamma[target] - 1.0) * X_expEffect[target] * X_covariate[target*Nx[0] + YtoX[record]]), 2.0);
						ProbExclude += pow((Y_expErrors[record] + X_expGamma[target] * X_expEffect[target] * X_covariate[target*Nx[0] + YtoX[record]]), 2.0);
					}
					/* In this DVR-Regression model, Y_variances[record] is added both to ProbInclude and ProbExclude. But bcecause these are cancelled out, these are not added explicitly.*/

					ProbInclude += X_x2[target] * X_varEffect[target];
					ProbInclude *= -0.5 * Tau0[0];	ProbExclude *= -0.5 * Tau0[0];
					ProbInclude += logPi;				ProbExclude += log1minusPi;

					temp = ProbInclude;
					if (temp < ProbExclude) temp = ProbExclude;
					ProbInclude -= temp;	ProbExclude -= temp;

					ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
					prop = ProbInclude / (ProbInclude + ProbExclude);
					prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

					/* update sumVarB and expErrors */
					sumVarB += (X_x2[target] * (prop2*X_exp2Effect[target] - pow(prop*X_expEffect[target], 2.0)));
					for (record = 0; record<Ny[0]; record++)
						Y_expErrors[record] += (X_covariate[target*Nx[0] + YtoX[record]] * X_expEffect[target] * (X_expGamma[target] - prop));

					/* for convergence check */
					Check1 += pow((prop - X_expGamma[target]), 2.0);
					Check2 += pow(prop, 2.0);

					X_expGamma[target] = prop;
					X_exp2Gamma[target] = prop2; /* used when variational Bayes */
				}
			}
			break;

		case 5: /* SSVS */
			for (locus = 0; locus<P[0]; locus++)
			{
				target = Order[locus];
				temp4 = X_exp2Effect[target];
				ProbInclude = -0.5 * temp4 / X_s2[0] + logPi;
				ProbExclude = -0.5 * temp4 / X_s2[0] * invC + log1minusPi - 0.5 * logC;

				temp = ProbInclude;
				if (temp < ProbExclude) temp = ProbExclude;
				ProbInclude -= temp;	ProbExclude -= temp;

				ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
				prop = ProbInclude / (ProbInclude + ProbExclude);
				prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

				/* for convergence check */
				Check1 += pow((prop - X_expGamma[target]), 2.0);
				Check2 += pow(prop, 2.0);

				X_expGamma[target] = prop;
				X_exp2Gamma[target] = prop2;
			}
			break;

		case 6: /* MIX */
			ProbIncludeConstant = Digamma(0.5*(H_v[0] + sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[0] + vS2)) + logPi;
			ProbExcludeConstant = Digamma(0.5*(H_v[0] + (double)P[0] - sumGamma[0])) - 0.5 * log(0.5*(sumGammaB2[1] + vS2)) + log1minusPi;
			/*---- Note -----------------------------------------------------------------------------------------------------
			sumGamma[0] and sumGammaB2 are not updated here, because these are the posterior parameters of X_expSigma[0].
			sumGamma[0] is replaced by sumGamma[1] after the update of all gamma.
			sumGammaB2 are updated when marker effects are updated.
			---------------------------------------------------------------------------------------------------------------*/
			for (locus = 0, sumGamma[1] = 0.0; locus<P[0]; locus++)
			{
				target = Order[locus];
				temp4 = X_exp2Effect[target];
				ProbInclude = ProbIncludeConstant - 0.5 * temp4 / X_s2[0];
				ProbExclude = ProbExcludeConstant - 0.5 * temp4 / X_s2[1];

				temp = ProbInclude;
				if (temp < ProbExclude) temp = ProbExclude;
				ProbInclude -= temp;	ProbExclude -= temp;

				ProbInclude = exp(ProbInclude);	ProbExclude = exp(ProbExclude);
				prop = ProbInclude / (ProbInclude + ProbExclude);
				prop2 = pow(prop, 2.0) + prop * (1.0 - prop);

				/* for convergence check */
				Check1 += pow((prop - X_expGamma[target]), 2.0);
				Check2 += pow(prop, 2.0);

				X_expGamma[target] = prop;
				X_exp2Gamma[target] = prop2;
				sumGamma[1] += X_expGamma[target];
			}
			sumGamma[0] = sumGamma[1];
			break;
		}

		if (Priortype[0] == 3 || Polygenic[0] == 1)
		{
			/* update of BV (GBLUP) */
			temp = S2 * Tau0[0];
			for (level = 0; level < Nx[0]; level++)
			{
				D[level] = S2 / (Eval[level] + temp);
				Yr[level] = Y_expErrors[XtoY[level]] + U_expBv[level];
			}

			Innerproduct_tADA(tEvec, D, Nx[0], Nx[0], U_varBv);
			Innerproduct_tAB(U_varBv, Yr, Nx[0], Nx[0], 1, prop3);
			for (level = 0; level < Nx[0]; level++)
			{
				prop3[level] *= Tau0[0];
				Y_expErrors[XtoY[level]] += (U_expBv[level] - prop3[level]);
				Check1 += pow((prop3[level] - U_expBv[level]), 2.0);
				Check2 += pow(prop3[level], 2.0);
				U_expBv[level] = prop3[level];
				sumVarB += U_varBv[level * Nx[0] + level];
			}

			/* update of Sigma2 (GBLUP) */
			temp2 = vS2_poly;
			for (level = 0; level < Nx[0]; level++)
			{
				temp = 0.0;
				for (level2 = 0; level2 < Nx[0]; level2++)
				{
					temp += U_expBv[level2] * Evec[level * Nx[0] + level2];
				}
				temp2 += Eval[level] * temp * temp;
				temp2 += (Eval[level] * S2) / (S2 * Tau0[0] + Eval[level]);
			}
			S2 = temp2 / vN;
			prop = vN * S2 / (vN - 2.0);
			Check1 += pow((prop - U_expSigma2[0]), 2.0);
			Check2 += pow(prop, 2.0);
			U_expSigma2[0] = prop;
			U_varSigma2[0] = 2.0 * pow(vN, 2.0) * pow(S2, 2.0) / (pow(vN - 2.0, 2.0) * (vN - 4.0));
		}

		/* update of Re */
		for (record = 0, temp = 0.0; record<Ny[0]; record++)
		{
			temp += (pow(Y_expErrors[record], 2.0) + Y_variance[record]);
		}
		/* "+Y_variance[record]" was added */

		if (Priortype[0] == 1 && CondResidual[0] == 1)
		{
			a1 = (double)(P[0] + Ny[0])*0.5;
			b1 = 0.5*(temp + sumVarB + sumTau2B2);
		}
		else
		{
			a1 = (double)Ny[0] * 0.5;
			b1 = 0.5*(temp + sumVarB);
		}
		prop = a1 / b1;
		Tau0[1] = a1 / pow(b1, 2.0);

		Check1 += pow((prop - Tau0[0]), 2.0);
		Check2 += pow(prop, 2.0);
		Tau0[0] = prop;

		Niterations[1] = ite;
		temp = Check1 / Check2;
		if (ite == Niterations[0] || temp<Threshold[0]) break;
	}	/* ite */

	if (Priortype[0] == 3 || Polygenic[0] == 1)
	{
		free(prop3); free(Yr); free(D);
	}
}

static const
R_CMethodDef cMethods[] = {
	{ "GenomeWideRegression", (DL_FUNC)&GenomeWideRegression, 50 }, {NULL, NULL, 0}
};

void R_init_GenomeBasedModel(DllInfo *info)
{
	/* Register the .C and .Call routines.
	No .Fortran() or .External() routines,
	so pass those arrays as NULL.
	*/
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
