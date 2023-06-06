#include <Rcpp.h>
using namespace Rcpp;

/*
 Copyright (C) 2018 Akio Onogi

Released under the MIT license
http://opensource.org/licenses/mit-license.php
*/

//Caluculate the exponential of likelihood
double Loglikelihood(NumericVector output, NumericVector y, NumericVector ve, int ne, double missing, int nve, NumericVector residualgroup, NumericVector loglikelihoodvec,
                     NumericVector y_uncertainty, NumericVector ve_proposal, NumericVector loglikelihoodvec_proposedve)
{
  int	env, pos;
  double loglikelihood = 0.0, sl;

  for (env = 0; env < nve; ++env)
  {
    loglikelihoodvec(env) = 0.0;
    loglikelihoodvec_proposedve(env) = 0.0;
  }
  for (env = 0; env < ne; ++env)
  {
    if (y(env) != missing)
    {
      pos = (int)residualgroup(env);
      sl = pow((y(env) - output(env)), 2.0);
      loglikelihoodvec(pos) -= sl * (0.5 / (ve(pos) + y_uncertainty(env)));
      loglikelihoodvec_proposedve(pos) -= sl * (0.5 / (ve_proposal(pos) + y_uncertainty(env)));
    }
  }
  for (env = 0; env < nve; ++env)
  {
    loglikelihood += loglikelihoodvec(env);
  }
  return(loglikelihood);
}

//log density of normal distribution
double	LogDensityKernel_Normal(double v, double mu, double var)
{
  return(-0.5 * pow((v - mu), 2.0) / var);
}

//[[Rcpp::export]]
int ModelParameterInference(NumericMatrix Input, NumericVector Freevec, NumericMatrix Y, double Missing, NumericMatrix Y_uncertainty,
                            int Nl, int Ne, int Np, int Nr,
                            NumericMatrix SampledPara, NumericMatrix SampledVe,
                            NumericMatrix MeanOfPrior, NumericVector VarOfPrior, NumericVector Ve, int Nve, NumericVector Residualgroup,
                            NumericVector SdforParameters, NumericVector SdforVe, NumericVector Loglike,
                            int Burnin, int Thi, int Totalite, IntegerVector Order,
                            NumericVector Lowerlimit, NumericVector Upperlimit, NumericVector AccParameters, NumericVector AccVe, Function Model,
                            NumericVector NeNl, NumericVector Transformation,
                            NumericMatrix nRandom_para, NumericMatrix nRandom_ve, NumericMatrix uRandom_para, NumericMatrix uRandom_ve)
{

  // For repeat statement
  int		mcmc, line, env, ii, para;

  // Model parameters
  NumericMatrix	Parameters(Np, Nl);

  // Output of the model
  NumericMatrix	Output_current(Ne,Nl), Output_proposed(Ne,Nl);

  // Sampling
  int		CumNs, CumNsab, DoSampling, AfterBurnin, Nsab;

  // MH update
  double  MHprob, Sd, Mu, Sigma, Proposedloglike, Currentloglike, temp;
  NumericVector  proposal(Nve), proposedpara(Np), currentpara(Np), ProposedloglikeVec(Nve), CurrentloglikeVec(Nve), SumCurrentloglikeVec(Nve), SumProposedloglikeVec(Nve), Temp(Ne);
  NumericVector  ProposedloglikeVec_ProposedVe(Nve), CurrentloglikeVec_ProposedVe(Nve);
  int   target, tf;

  // Initialization
  for (line = 0; line < Nl; ++line)
  {
    for (para = 0; para < Np; ++para)
    {
      Parameters(para, line) = nRandom_para(0, para * Nl + line) * sqrt(VarOfPrior(para)) + MeanOfPrior(para, line);
      if(Parameters(para, line)<=Lowerlimit(para)) Parameters(para, line) = Lowerlimit(para) + 1e-5;
      if(Parameters(para, line)>=Upperlimit(para)) Parameters(para, line) = Upperlimit(para) - 1e-5;
    }
    Temp = Model(Input(_,line), Freevec, Parameters(_, line));
    for(env = 0; env < Ne; ++env)
    {
      Output_current(env,line) = Temp (env);
    }
  }

  CumNs = -1;
  CumNsab = -1;
  Nsab = (Totalite - Burnin) / Thi;

  for (mcmc = 1; mcmc <= Totalite; ++mcmc)
  {
    //Sample or Not
    if (mcmc % Thi == 0)
    {
      DoSampling = 1; CumNs++;
    }
    else
    {
      DoSampling = 0;
    }
    if (DoSampling&&mcmc>Burnin)
    {
      AfterBurnin = 1; CumNsab++;
    }
    else
    {
      AfterBurnin = 0;
    }

    //propose residual variances
    for (env = 0; env < Nve; ++env)
    {
      proposal(env) = nRandom_ve(mcmc-1, env) * SdforVe(env) + Ve(env);
    }

    //Update of model parameters
    for (para = 0; para < Np; ++para)
    {
      for (env = 0; env < Nve; ++env)
      {
        SumCurrentloglikeVec(env) = 0.0;
        SumProposedloglikeVec(env) = 0.0;
      }
      for (line = 0; line < Nl; ++line)
      {
        target = Order(line);

        Sd = SdforParameters(para);
        Mu = MeanOfPrior(para, target);
        Sigma = VarOfPrior(para);
        for (ii = 0; ii < Np; ++ii)
        {
          currentpara(ii) = Parameters(ii, target);
          proposedpara(ii) = Parameters(ii, target);
        }

        MHprob = 0.0;
        proposedpara(para) = nRandom_para(mcmc, para * Nl + target) * Sd + currentpara(para);

        if(proposedpara(para)>Lowerlimit(para)&&proposedpara(para)<Upperlimit(para))
        {
          Temp = Model(Input(_,target), Freevec, proposedpara);
          for(env = 0; env < Ne; ++env)
          {
            Output_proposed(env,target) = Temp (env);
          }

          Currentloglike = Loglikelihood(Output_current(_,target), Y(_,target), Ve, Ne, Missing, Nve, Residualgroup, CurrentloglikeVec, Y_uncertainty, proposal, CurrentloglikeVec_ProposedVe);
          Proposedloglike = Loglikelihood(Output_proposed(_,target), Y(_,target), Ve, Ne, Missing, Nve, Residualgroup, ProposedloglikeVec, Y_uncertainty, proposal, ProposedloglikeVec_ProposedVe);

          MHprob += Proposedloglike;
          MHprob -= Currentloglike;

          MHprob += LogDensityKernel_Normal(proposedpara(para), Mu, Sigma);
          MHprob -= LogDensityKernel_Normal(currentpara(para), Mu, Sigma);

          tf = (int)Transformation(para);
          switch (tf){
          case 1:
            MHprob -= currentpara(para);
            MHprob += proposedpara(para);
            break;
          case 2:
            temp = exp(currentpara(para))/(1.0 + exp(currentpara(para)));
            MHprob -= log((temp*(1-temp)));
            temp = exp(proposedpara(para))/(1.0 + exp(proposedpara(para)));
            MHprob += log((temp*(1-temp)));
            break;
          }

          if (MHprob>log(uRandom_para(mcmc-1, para * Nl + target)))
          {
            Parameters(para, target) = proposedpara(para);
            AccParameters(para) += 1.0;
            for (env = 0; env < Ne; ++env)
            {
              Output_current(env,target) = Output_proposed(env, target);
            }
            for (env = 0; env < Nve; ++env)
            {
              SumCurrentloglikeVec(env) += ProposedloglikeVec(env);
              SumProposedloglikeVec(env) += ProposedloglikeVec_ProposedVe(env);
            }
          }
          else
          {
            if(para == Np - 1)
            {
              for (env = 0; env < Nve; ++env)
              {
                SumCurrentloglikeVec(env) += CurrentloglikeVec(env);
                SumProposedloglikeVec(env) += CurrentloglikeVec_ProposedVe(env);
              }
            }
          }
        }
        else
        {
          if(para == Np - 1)
          {
            Currentloglike = Loglikelihood(Output_current(_,target), Y(_,target), Ve, Ne, Missing, Nve, Residualgroup, CurrentloglikeVec, Y_uncertainty, proposal, CurrentloglikeVec_ProposedVe);
            for (env = 0; env < Nve; ++env)
            {
              SumCurrentloglikeVec(env) += CurrentloglikeVec(env);
              SumProposedloglikeVec(env) += CurrentloglikeVec_ProposedVe(env);
            }
          }
        }

        if (AfterBurnin)
        {
          SampledPara(para * Nsab + CumNsab, target) = Parameters(para,target);
        }
      }
    }

    //Update of Ve
    for(env = 0; env<Nve; ++env)
    {
      Currentloglike = SumCurrentloglikeVec(env) - 0.5 * NeNl(env) * log(Ve(env));
      if (proposal(env)>0.0)
      {
        Proposedloglike = SumProposedloglikeVec(env) - 0.5 * NeNl(env) * log(proposal(env));

        MHprob = Proposedloglike - Currentloglike;
        MHprob -= log(proposal(env));
        MHprob += log(Ve(env));

        if (MHprob>log(uRandom_ve(mcmc-1, env)))
        {
          Ve(env) = proposal(env);
          AccVe(env) += 1.0;
          if (DoSampling) { Loglike(CumNs) = Proposedloglike; }
        }
        else
        {
          if (DoSampling) { Loglike(CumNs) = Currentloglike; }
        }
      }
      else
      {
        if (DoSampling) { Loglike(CumNs) = Currentloglike; }
      }
      if (AfterBurnin)
      {
        SampledVe(CumNsab,env) = Ve(env);
      }
    }
  }	//mcmc

  return 0;
}

//[[Rcpp::export]]
int ModelParameterInference_PassMatrix(NumericMatrix Input, NumericVector Freevec, NumericMatrix Y, double Missing, NumericMatrix Y_uncertainty,
                                       int Nl, int Ne, int Np, int Nr,
                                       NumericMatrix SampledPara, NumericMatrix SampledVe,
                                       NumericMatrix MeanOfPrior, NumericVector VarOfPrior, NumericVector Ve, int Nve, NumericVector Residualgroup,
                                       NumericVector SdforParameters, NumericVector SdforVe, NumericVector Loglike,
                                       int Burnin, int Thi, int Totalite, IntegerVector Order,
                                       NumericVector Lowerlimit, NumericVector Upperlimit, NumericVector AccParameters, NumericVector AccVe, Function Model,
                                       NumericVector NeNl, NumericVector Transformation,
                                       NumericMatrix nRandom_para, NumericMatrix nRandom_ve, NumericMatrix uRandom_para, NumericMatrix uRandom_ve)
{

  // For repeat statement
  int		mcmc, line, env, ii, para;

  // Model parameters
  NumericMatrix	Parameters(Np, Nl);

  // Output of the model
  NumericMatrix	Output_current(Ne,Nl), Output_proposed(Ne,Nl);

  // Sampling
  int		CumNs, CumNsab, DoSampling, AfterBurnin, Nsab;

  // MH update
  double	MHprob, Sd, Mu, Sigma, Proposedloglike, Currentloglike, temp;
  NumericVector  proposal(Nve), ProposedloglikeVec(Nve), CurrentloglikeVec(Nve), SumCurrentloglikeVec(Nve), SumProposedloglikeVec(Nve);
  NumericVector  ProposedloglikeVec_ProposedVe(Nve), CurrentloglikeVec_ProposedVe(Nve);
  NumericMatrix  proposedpara(Np, Nl), currentpara(Np, Nl);
  int   target, tf;

  // Initialization
  for (line = 0; line < Nl; ++line)
  {
    for (para = 0; para < Np; ++para)
    {
      Parameters(para, line) = nRandom_para(0, para * Nl + line) * sqrt(VarOfPrior(para)) + MeanOfPrior(para, line);
      if(Parameters(para, line)<=Lowerlimit(para)) Parameters(para, line) = Lowerlimit(para) + 1e-5;
      if(Parameters(para, line)>=Upperlimit(para)) Parameters(para, line) = Upperlimit(para) - 1e-5;
    }
  }
  Output_current = Model(Input, Freevec, Parameters);

  CumNs = -1;
  CumNsab = -1;
  Nsab = (Totalite - Burnin) / Thi;
  for (mcmc = 1; mcmc <= Totalite; ++mcmc)
  {
    //Sample or Not
    if (mcmc % Thi == 0)
    {
      DoSampling = 1; CumNs++;
    }
    else
    {
      DoSampling = 0;
    }
    if (DoSampling&&mcmc>Burnin)
    {
      AfterBurnin = 1; CumNsab++;
    }
    else
    {
      AfterBurnin = 0;
    }

    //propose residual variances
    for (env = 0; env < Nve; ++env)
    {
      proposal(env) = nRandom_ve(mcmc-1, env) * SdforVe(env) + Ve(env);
    }

    //Update of model parameters
    for (para = 0; para < Np; ++para)
    {
      for (env = 0; env < Nve; ++env)
      {
        SumCurrentloglikeVec(env) = 0.0;
        SumProposedloglikeVec(env) = 0.0;
      }
      for (line = 0; line < Nl; ++line)
      {
        target = Order(line);

        Sd = SdforParameters(para);
        for (ii = 0; ii < Np; ++ii)
        {
          currentpara(ii, target) = Parameters(ii, target);
          proposedpara(ii, target) = Parameters(ii, target);
        }
        proposedpara(para, target) = nRandom_para(mcmc, para * Nl + target) * Sd + currentpara(para, target);
        if(proposedpara(para, target)<Lowerlimit(para)||proposedpara(para, target)>Upperlimit(para))
          proposedpara(para, target) = currentpara(para, target);
      }
      Output_proposed = Model(Input, Freevec, proposedpara);

      for (line = 0; line < Nl; ++line)
      {
        target = Order(line);
        Mu = MeanOfPrior(para, target);
        Sigma = VarOfPrior(para);
        MHprob = 0.0;

        if(proposedpara(para, target)!=currentpara(para, target))
        {
          Currentloglike = Loglikelihood(Output_current(_,target), Y(_,target), Ve, Ne, Missing, Nve, Residualgroup, CurrentloglikeVec, Y_uncertainty, proposal, CurrentloglikeVec_ProposedVe);
          Proposedloglike = Loglikelihood(Output_proposed(_,target), Y(_,target), Ve, Ne, Missing, Nve, Residualgroup, ProposedloglikeVec, Y_uncertainty, proposal, ProposedloglikeVec_ProposedVe);

          MHprob += Proposedloglike;
          MHprob -= Currentloglike;

          MHprob += LogDensityKernel_Normal(proposedpara(para, target), Mu, Sigma);
          MHprob -= LogDensityKernel_Normal(currentpara(para, target), Mu, Sigma);

          tf = (int)Transformation(para);
          switch (tf){
          case 1:
            MHprob -= currentpara(para, target);
            MHprob += proposedpara(para, target);
            break;
          case 2:
            temp = exp(currentpara(para, target))/(1.0 + exp(currentpara(para, target)));
            MHprob -= log((temp*(1-temp)));
            temp = exp(proposedpara(para, target))/(1.0 + exp(proposedpara(para, target)));
            MHprob += log((temp*(1-temp)));
            break;
          }

          if (MHprob>log(uRandom_para(mcmc-1, para * Nl + target)))
          {
            Parameters(para, target) = proposedpara(para, target);
            AccParameters(para) += 1.0;
            for (env = 0; env < Ne; ++env)
            {
              Output_current(env, target) = Output_proposed(env, target);
            }
            for (env = 0; env < Nve; ++env)
            {
              SumCurrentloglikeVec(env) += ProposedloglikeVec(env);
              SumProposedloglikeVec(env) += ProposedloglikeVec_ProposedVe(env);
            }
          }
          else
          {
            if(para == Np - 1)
            {
              for (env = 0; env < Nve; ++env)
              {
                SumCurrentloglikeVec(env) += CurrentloglikeVec(env);
                SumProposedloglikeVec(env) += CurrentloglikeVec_ProposedVe(env);
              }
            }
          }
        }
        else
        {
          if(para == Np - 1)
          {
            Currentloglike = Loglikelihood(Output_current(_,target), Y(_,target), Ve, Ne, Missing, Nve, Residualgroup, CurrentloglikeVec, Y_uncertainty, proposal, CurrentloglikeVec_ProposedVe);
            for (env = 0; env < Nve; ++env)
            {
              SumCurrentloglikeVec(env) += CurrentloglikeVec(env);
              SumProposedloglikeVec(env) += CurrentloglikeVec_ProposedVe(env);
            }
          }
        }

        if (AfterBurnin)
        {
          SampledPara(para * Nsab + CumNsab, target) = Parameters(para, target);
        }
      }//line
    }//para

    //Update of Ve
    for(env = 0; env<Nve; ++env)
    {
      Currentloglike = SumCurrentloglikeVec(env) - 0.5 * NeNl(env) * log(Ve(env));
      if (proposal(env)>0.0)
      {
        Proposedloglike = SumProposedloglikeVec(env) - 0.5 * NeNl(env) * log(proposal(env));

        MHprob = Proposedloglike - Currentloglike;
        MHprob -= log(proposal(env));
        MHprob += log(Ve(env));

        if (MHprob>log(uRandom_ve(mcmc-1, env)))
        {
          Ve(env) = proposal(env);
          AccVe(env) += 1.0;
          if (DoSampling) { Loglike(CumNs) = Proposedloglike; }
        }
        else
        {
          if (DoSampling) { Loglike(CumNs) = Currentloglike; }
        }
      }
      else
      {
        if (DoSampling) { Loglike(CumNs) = Currentloglike; }
      }
      if (AfterBurnin)
      {
        SampledVe(CumNsab,env) = Ve(env);
      }
    }
  }	//mcmc

  return 0;
}

static const
  R_CallMethodDef callMethods[] = {
    { "ModelParameterInference", (DL_FUNC)&ModelParameterInference, 34}, {NULL, NULL, 0},
    { "ModelParameterInference_PassMatrix", (DL_FUNC)&ModelParameterInference_PassMatrix, 34}, {NULL, NULL, 0}
  };

void R_init_GenomeBasedModel(DllInfo *info)
{
  /* Register the .C and .Call routines.
   No .Fortran() or .External() routines,
   so pass those arrays as NULL.
   */
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

