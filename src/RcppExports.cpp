// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ModelParameterInference
int ModelParameterInference(NumericMatrix Input, NumericVector Freevec, NumericMatrix Y, double Missing, NumericMatrix Y_uncertainty, int Nl, int Ne, int Np, int Nr, NumericMatrix SampledPara, NumericMatrix SampledVe, NumericMatrix MeanOfPrior, NumericVector VarOfPrior, NumericVector Ve, int Nve, NumericVector Residualgroup, NumericVector SdforParameters, NumericVector SdforVe, NumericVector Loglike, int Burnin, int Thi, int Totalite, IntegerVector Order, NumericVector Lowerlimit, NumericVector Upperlimit, NumericVector AccParameters, NumericVector AccVe, Function Model, NumericVector NeNl, NumericVector Transformation, NumericMatrix nRandom_para, NumericMatrix nRandom_ve, NumericMatrix uRandom_para, NumericMatrix uRandom_ve);
RcppExport SEXP _GenomeBasedModel_ModelParameterInference(SEXP InputSEXP, SEXP FreevecSEXP, SEXP YSEXP, SEXP MissingSEXP, SEXP Y_uncertaintySEXP, SEXP NlSEXP, SEXP NeSEXP, SEXP NpSEXP, SEXP NrSEXP, SEXP SampledParaSEXP, SEXP SampledVeSEXP, SEXP MeanOfPriorSEXP, SEXP VarOfPriorSEXP, SEXP VeSEXP, SEXP NveSEXP, SEXP ResidualgroupSEXP, SEXP SdforParametersSEXP, SEXP SdforVeSEXP, SEXP LoglikeSEXP, SEXP BurninSEXP, SEXP ThiSEXP, SEXP TotaliteSEXP, SEXP OrderSEXP, SEXP LowerlimitSEXP, SEXP UpperlimitSEXP, SEXP AccParametersSEXP, SEXP AccVeSEXP, SEXP ModelSEXP, SEXP NeNlSEXP, SEXP TransformationSEXP, SEXP nRandom_paraSEXP, SEXP nRandom_veSEXP, SEXP uRandom_paraSEXP, SEXP uRandom_veSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Input(InputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Freevec(FreevecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type Missing(MissingSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_uncertainty(Y_uncertaintySEXP);
    Rcpp::traits::input_parameter< int >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< int >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< int >::type Np(NpSEXP);
    Rcpp::traits::input_parameter< int >::type Nr(NrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SampledPara(SampledParaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SampledVe(SampledVeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type MeanOfPrior(MeanOfPriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type VarOfPrior(VarOfPriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< int >::type Nve(NveSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Residualgroup(ResidualgroupSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SdforParameters(SdforParametersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SdforVe(SdforVeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Loglike(LoglikeSEXP);
    Rcpp::traits::input_parameter< int >::type Burnin(BurninSEXP);
    Rcpp::traits::input_parameter< int >::type Thi(ThiSEXP);
    Rcpp::traits::input_parameter< int >::type Totalite(TotaliteSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Order(OrderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lowerlimit(LowerlimitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Upperlimit(UpperlimitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type AccParameters(AccParametersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type AccVe(AccVeSEXP);
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NeNl(NeNlSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Transformation(TransformationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nRandom_para(nRandom_paraSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nRandom_ve(nRandom_veSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type uRandom_para(uRandom_paraSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type uRandom_ve(uRandom_veSEXP);
    rcpp_result_gen = Rcpp::wrap(ModelParameterInference(Input, Freevec, Y, Missing, Y_uncertainty, Nl, Ne, Np, Nr, SampledPara, SampledVe, MeanOfPrior, VarOfPrior, Ve, Nve, Residualgroup, SdforParameters, SdforVe, Loglike, Burnin, Thi, Totalite, Order, Lowerlimit, Upperlimit, AccParameters, AccVe, Model, NeNl, Transformation, nRandom_para, nRandom_ve, uRandom_para, uRandom_ve));
    return rcpp_result_gen;
END_RCPP
}
// ModelParameterInference_PassMatrix
int ModelParameterInference_PassMatrix(NumericMatrix Input, NumericVector Freevec, NumericMatrix Y, double Missing, NumericMatrix Y_uncertainty, int Nl, int Ne, int Np, int Nr, NumericMatrix SampledPara, NumericMatrix SampledVe, NumericMatrix MeanOfPrior, NumericVector VarOfPrior, NumericVector Ve, int Nve, NumericVector Residualgroup, NumericVector SdforParameters, NumericVector SdforVe, NumericVector Loglike, int Burnin, int Thi, int Totalite, IntegerVector Order, NumericVector Lowerlimit, NumericVector Upperlimit, NumericVector AccParameters, NumericVector AccVe, Function Model, NumericVector NeNl, NumericVector Transformation, NumericMatrix nRandom_para, NumericMatrix nRandom_ve, NumericMatrix uRandom_para, NumericMatrix uRandom_ve);
RcppExport SEXP _GenomeBasedModel_ModelParameterInference_PassMatrix(SEXP InputSEXP, SEXP FreevecSEXP, SEXP YSEXP, SEXP MissingSEXP, SEXP Y_uncertaintySEXP, SEXP NlSEXP, SEXP NeSEXP, SEXP NpSEXP, SEXP NrSEXP, SEXP SampledParaSEXP, SEXP SampledVeSEXP, SEXP MeanOfPriorSEXP, SEXP VarOfPriorSEXP, SEXP VeSEXP, SEXP NveSEXP, SEXP ResidualgroupSEXP, SEXP SdforParametersSEXP, SEXP SdforVeSEXP, SEXP LoglikeSEXP, SEXP BurninSEXP, SEXP ThiSEXP, SEXP TotaliteSEXP, SEXP OrderSEXP, SEXP LowerlimitSEXP, SEXP UpperlimitSEXP, SEXP AccParametersSEXP, SEXP AccVeSEXP, SEXP ModelSEXP, SEXP NeNlSEXP, SEXP TransformationSEXP, SEXP nRandom_paraSEXP, SEXP nRandom_veSEXP, SEXP uRandom_paraSEXP, SEXP uRandom_veSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Input(InputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Freevec(FreevecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type Missing(MissingSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_uncertainty(Y_uncertaintySEXP);
    Rcpp::traits::input_parameter< int >::type Nl(NlSEXP);
    Rcpp::traits::input_parameter< int >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< int >::type Np(NpSEXP);
    Rcpp::traits::input_parameter< int >::type Nr(NrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SampledPara(SampledParaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SampledVe(SampledVeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type MeanOfPrior(MeanOfPriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type VarOfPrior(VarOfPriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< int >::type Nve(NveSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Residualgroup(ResidualgroupSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SdforParameters(SdforParametersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SdforVe(SdforVeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Loglike(LoglikeSEXP);
    Rcpp::traits::input_parameter< int >::type Burnin(BurninSEXP);
    Rcpp::traits::input_parameter< int >::type Thi(ThiSEXP);
    Rcpp::traits::input_parameter< int >::type Totalite(TotaliteSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Order(OrderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lowerlimit(LowerlimitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Upperlimit(UpperlimitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type AccParameters(AccParametersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type AccVe(AccVeSEXP);
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NeNl(NeNlSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Transformation(TransformationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nRandom_para(nRandom_paraSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nRandom_ve(nRandom_veSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type uRandom_para(uRandom_paraSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type uRandom_ve(uRandom_veSEXP);
    rcpp_result_gen = Rcpp::wrap(ModelParameterInference_PassMatrix(Input, Freevec, Y, Missing, Y_uncertainty, Nl, Ne, Np, Nr, SampledPara, SampledVe, MeanOfPrior, VarOfPrior, Ve, Nve, Residualgroup, SdforParameters, SdforVe, Loglike, Burnin, Thi, Totalite, Order, Lowerlimit, Upperlimit, AccParameters, AccVe, Model, NeNl, Transformation, nRandom_para, nRandom_ve, uRandom_para, uRandom_ve));
    return rcpp_result_gen;
END_RCPP
}