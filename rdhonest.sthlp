{smcl}
{cmd:help rdhonest}{right: ({browse "https://github.com/SaiChrisZHANG/RDHonest-vStata":Stata Journal No.})}
{hline}
{* *! version 0.2 SZhang 13Aug2022}{...}

{viewerjumpto "Syntax" "examplehelpfile##syntax"}{...}
{viewerjumpto "Description" "examplehelpfile##description"}{...}
{viewerjumpto "Options" "examplehelpfile##options"}{...}
{viewerjumpto "Remarks" "examplehelpfile##remarks"}{...}
{viewerjumpto "Examples" "examplehelpfile##examples"}{...}
{title:Title}

{phang}
{bf:rdhonest} {hline 2} Honestly and nearly-optimally calculate estimators, standard errors and bias-aware CIs.


{marker syntax}{...}
{title:Syntax}

{phang}
For sharp RD estimation: {p_end}

{p 8 14 2}
{cmdab:rdhonest} {depvar}
{it:{help varname:runvar}}
[{it:{help varlist:covar}}]
{ifin}
[{cmd:,} {it:options}]

{phang}
For fuzzy RD estimation: {p_end}

{p 8 14 2}
{cmdab:rdhonest} {depvar}
{cmd:(}{it:{help varname:treat}}{cmd: =}
    {it:{help varname:runvar}}{cmd:)} 
[{it:{help varlist:covar}}]
{ifin}
[{cmd:,} {it:options}]

{phang}
{it:runvar} is the running variable.{p_end}

{phang}
{it:treat} is the treatment variable (used in fuzzy RD estimation).{p_end}

{phang}
{it:covar} is the list of covariates.

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt c(#)}}set regression discontinuity cutoff; default is {cmd:c(0)}{p_end}
{p2coldent:* {opt m(#)}}set bound on the second derivative of the conditional mean function; if not specified, use Armstrong and Kolesar (2020) rule of thumb to calculate smoothness bound{p_end}
{synopt:{opt kernel(string)}}set kernel function used in the local regression; default is {cmd:{ul:tri}angular}; other options are {cmd:{ul:epa}nechnikov} and {cmd:{ul:uni}form}{p_end}
{synopt:{opt opt_criterion(string)}}set the optimality criterion that bandwidth is designed to optimize; default is {cmd:MSE}; other options are {cmd:OCI} and {cmd:FLCI}{p_end}
{synopt:{opt h(#)}}set the bandwidth; if not specified, optimal bandwidth computed according to {cmd:opt_criterion()} is used{p_end}
{synopt:{opt se_method(string)}}set method used to compute standard errors; default is {cmd:NN}, other options are {cmd:EHW}, ignored if {cmd:{ul:pvar}iance()} is specified{p_end}
{synopt:{opt alpha(#)}}set the worst-case coverage (1-alpha) of the confidence interval; default alpha is {cmd:c(0.05)}{p_end}

{syntab:Supporting}
{synopt:{opt beta(#)}}determine quantile of excess length to optimize if {cmd:opt_criterion()} is specifed as {cmd:OCI}; default is {cmd:c(0.8)}{p_end}
{synopt:{opt j(#)}}set number of nearest neighbors used if {cmd:se_method()} is specified as {cmd:NN}; default is {cmd:c(3)}{p_end}
{synopt:{opt pvar:iance(varlist)}}supply preliminary conditional variance{p_end}
{p2coldent:* {opt t0(#)}}set initial estimate of the treatment effect for calculating the optimal bandwidth; default is {cmd:c(0)}{p_end}

{syntab:Other}
{synopt:{opt cl:uster(varlist)}}supply list of clustering variables{p_end}
{synopt:{opt w:eight(varname)}}supply the weighting variable{p_end}

{syntab:Reporting}
{synopt:{opt savew:gtest(newvar)}}save the estimated weights as a new variable named {help newvar}; when combined with by, only save the estimated weights for the last by-group{p_end}
{synopt:{opt noparam:eter}}do not display parameters used for estimation{p_end}
{synopt:{opt iterl:og}}display the iteration log{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}{cmd:m(#)} requires two scalar inputs for Fuzzy RD estimation, {cmd:m(#1 #2)}, where {cmd:#1} is for the first-stage regression and {cmd:#2} is for the reduced-form regression.{p_end}
{p 4 6 2}{cmd:t0(#)} is only relevant for Fuzzy RD estimations.{p_end}
{p 4 6 2}{cmd:by} is allowed; see {helpb by} for details.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rdhonest} calculates honest and nearly-optimal one- and two-sided confidence intervals in fuzzy and sharp regression discontinuity designs based on local linear regression. 


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt c(#)} sets the regression discontinuity cutoff. The default value for the cutoff is {cmd:0}.

{phang}
{opt m(#)} sets bound on the second derivative of the conditional mean function. If not specified, a global quartic regression is used for a rule-of-thumb estimation, following Armstrong and Kolesar (2020).

{pmore}
For fuzzy RD, {cmd:m(#1 #2)} sets the bound for both the first-stage ({cmd:#1}) and the reduced-form regression ({cmd:#2}).

{phang}
{opt kernel(string)} sets the kernel for the local regression (in the context of standard sharp RD). The following options are allowed:

{p 12 16 4}
{opt tri:angular} kernel (default): k(u) = (1-|u|)_{+};

{p 12 16 4}
{opt epa:nechnikov} kernel: k(u) = 0.75(1-u^2)_{+}; and

{p 12 16 4}
{opt uni:form} kernel: k(u) = (|u|<1)/2.

{phang}
{opt opt_criterion(string)} sets the optimality criterion used to compute the optimal bandwidth. The following options are allowed:

{p 12 16 4}
{opt MSE} (default): the finite-sample maximum mean-squared error;

{p 12 16 4}
{opt OCI}: the given quantile of excess length of one-sided confidence intervals; and

{p 12 16 4}
{opt FLCI}: the length of two-sided fixed-length confidence intervals.

{phang}
{opt h(#)} sets the bandwidth. If not specified, optimal bandwidth is computed according to criterion given by {opt opt_criterion()}.

{phang}
{opt se_method(string)} sets the method used to compute standard errors when running local regressions. The following options are allowed:

{p 12 16 4}
{opt NN} (default): Nearest neighbor method;

{p 12 16 4}
{opt EHW}: Eicker-Huber-White, with residuals from local regressions (local polynomial estimators only).

{pmore}
If {opt pvar:iance(varlist)} is specified, any input for {opt se_method()} will be ignored. Instead, {cmd:"supplied_var"} will be stored in {helpb e()}.

{phang}
{opt alpha(#)} determines confidence level, 1-{cmd:alpha}, for constructing/optimizing confidence intervals, default value is {cmd:0.05}.

{dlgtab:Supporting}

{phang}
{opt beta(#)} sets the quantile of excess length, default value {cmd:0.8}. Only relevant when the one-sided confidence interval criterion ({cmd:OCI}) is specified for {opt opt_criterion()}, ignored otherwise.

{phang}
{opt j(#)} sets the number of nearest neighbors used to compute standard errors when {cmd: NN} is specified for {opt se_method()}, default value is {cmd:3}.

{phang}
{opt pvar:iance(varlist)} supplies conditional variance for estimation. When specified, any input for {opt se_method()} is ignored, {cmd:"supplied_var"} will be stored in {helpb e()}.

{phang}
{opt t0(#)} sets the initial estimate of the treatment effect for calculating the optimal bandwidth, default value is {cmd:0}. Only relevant for Fuzzy RD.

{dlgtab:Other}

{phang}
{opt cl:uster(varlist)} supplies the list of clustering variables.{p_end}

{phang}
{opt w:eight(varname)} supplies the weighting variable.{p_end}

{dlgtab:Reporting}

{phang}
{opt savew:gtest(newvar)} saves the estimation weights as a new variable named {help newvar}; when combined with {helpb by}, only the estimation weights for the last by-group are saved, 0 will be saved for other by-groups.{p_end}

{phang}
{opt noparam:eter} specifies that the parameters used for estimation are not displayed.{p_end}

{phang}
{opt iterl:og} specifies that the iteration log is displayed.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}{cmd:Sharp} RD examples{p_end}

{phang2}Setup{p_end}
{phang3}{cmd:. webuse lee08}{p_end}

{phang2}Bandwidth optimized, with uniform kernel{p_end}
{phang3}{cmd:. rdhonest voteshare margin, m(0.1) kernel("uni") }{p_end}

{phang2}Bandwidth optimized, smoothness bound estimated, with triangular kernel; estimated weights saved as {it:wgt}{p_end}
{phang3}{cmd:. rdhonest voteshare margin, kernel("tri") savew(wgt)}{p_end}

{phang2}Bandwidth optimized, with triangular kernel; display iteration log, not display parameters used{p_end}
{phang3}{cmd:. rdhonest voteshare margin, m(0.1) kernel("tri") noparam iterl}{p_end}

{pstd}{cmd:Fuzzy} RD examples{p_end}

{phang2}Setup{p_end}
{phang3}{cmd:. webuse rcp}{p_end}

{phang2}Bandwidth optimized, with uniform kernel{p_end}
{phang3}{cmd:. rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") t0(0)}{p_end}

{phang2}Bandwidth optimized, smoothness bound estimated, with triangular kernel; estimation weights saved as {it:wgt}{p_end}
{phang3}{cmd:. rdhonest cn (retired=elig_year), kernel("tri") t0(0) savew(wgt)}{p_end}

{phang2}Bandwidth optimized, with triangular kernel; display iteration log, not display parameters used{p_end}
{phang3}{cmd:. rdhonest cn (retired=elig_year), m(4 0.4) kernel("tri") t0(0) noparam iterl}{p_end}


{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:rdhonest} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(M)}}smoothness bound for sharp RD estimations, or reduced-form smoothness bound for fuzzy RD estimations{p_end}
{synopt:{cmd:e(M_fs)}}frst-stage smoothness bound for Fuzzy RD estimations{p_end}
{synopt:{cmd:e(bandwidth)}}bandwidth, supplied or optimized{p_end}
{synopt:{cmd:e(alpha)}}confidence level{p_end}
{synopt:{cmd:e(beta)}}quantile of excess length for evaluating minimax efficiency of one-sided CIs{p_end}
{synopt:{cmd:e(est)}}sharp RD point estimation, or fuzzy RD reduced-form point estimation{p_end}
{synopt:{cmd:e(est_fs)}}fuzzy RD first-stage point estimation{p_end}
{synopt:{cmd:e(se)}}standard errors of {cmd:e(est)}{p_end}
{synopt:{cmd:e(bias)}}maximum bias of {cmd:e(est)}{p_end}
{synopt:{cmd:e(HLCi)}}half-length of a two-sided CI based on {cmd:e(est)}, such that the CI is given by ({cmd:e(est)}-{cmd:e(HLCi)}, {cmd:e(est)}+{cmd:e(HLCi)}){p_end}
{synopt:{cmd:e(TCiL)}}= {cmd:e(est)}+{cmd:e(HLCi)}, lower bound of a two-sided CI{p_end}
{synopt:{cmd:e(TCiU)}}= {cmd:e(est)}-{cmd:e(HLCi)}, upper bound of a two-sided CI{p_end}
{synopt:{cmd:e(OCiL)}}lower end-point of a one-sided CI based on {cmd:e(est)}{p_end}
{synopt:{cmd:e(OCiU)}}upper end-point of a one-sided CI based on {cmd:e(est)}{p_end}
{synopt:{cmd:e(effObs)}}effective number of observations used for estimation{p_end}
{synopt:{cmd:e(Leverage)}}maximum leverage for estimation {cmd:e(est)}{p_end}
{synopt:{cmd:e(cutoff)}}RD cutoff of the running variable{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rdhonest}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(title)}}title of output{p_end}
{synopt:{cmd:e(se_method)}}method for estimating standard error of estimate{p_end}
{synopt:{cmd:e(opt_crit)}}optimality criterion that bandwidth is designed to optimize{p_end}
{synopt:{cmd:e(kernel)}}kernel function used in the local regression{p_end}
{synopt:{cmd:e(rd)}}regression discontinuity class (sharp or fuzzy){p_end}
{synopt:{cmd:e(savewgtest)}}name of the new variable that stores estimated weights{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(runvar)}}name of running variable{p_end}
{synopt:{cmd:e(treat)}}name of treatment variable, for fuzzy RD{p_end}
{synopt:{cmd:e(covar)}}list of covariates{p_end}
{synopt:{cmd:e(cluster)}}list of clustering variables{p_end}
{synopt:{cmd:e(weight)}}name of weighting variable{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker AK2018a}{...}
{phang}
Armstrong, T. B., & Kolesár, M. (2018). Optimal inference in a class of regression models. {it:Econometrica}, 86, 655-683. {browse "https://doi.org/10.3982/ECTA14434":doi:10.3982/ECTA14434}{p_end}

{marker AK2020}{...}
{phang}
Armstrong, T. B., & Kolesár, M. (2020). Simple and honest confidence intervals in nonparametric regression. {it:Quantitative Economics}, 11(1), 1-39. {browse "https://doi.org/10.3982/QE1199":doi:10.3982/QE1199}{p_end}

{marker IK2012}{...}
{phang}
Imbens, G., & Kalyanaraman, K. (2012). Optimal bandwidth choice for the regression discontinuity estimator. {it:The Review of economic studies}, 79(3), 933-959. {browse "https://doi.org/10.1093/restud/rdr043":doi:10.1093/restud/rdr043}{p_end}

{marker KR2018}{...}
{phang}
Kolesár, M., & Rothe, C. (2018). Inference in regression discontinuity designs with a discrete running variable. {it:American Economic Review}, 108, 2277–2304. {browse "https://doi.org/10.1257/aer.20160945":doi:10.1257/aer.20160945}{p_end}


{marker citation}
{title:Cite rdhonest as follow}

{phang}
Citation for rdhonest, the Stata package.{p_end}


{marker author}{...}
{title:Authors}

{pstd}Timothy Armstrong{p_end}
{pstd}University of Southern California{p_end}

{pstd}Michal Kolesár{p_end}
{pstd}Princeton University{p_end}

{pstd}Yugen Chen{p_end}
{pstd}University of Southern California{p_end}

{pstd}Sai Zhang{p_end}
{pstd}University of Southern California{p_end}

{pstd}Kwok-Hao Lee{p_end}
{pstd}Princeton University{p_end}


{marker also_see}{...}
{title:Also see}

{p 4 14 2}
Development version: net install rdhonest, from("https://raw.githubusercontent.com/SaiChrisZHANG/RDHonest-vStata/master/current/"){p_end}

{p 4 14 2}
Article: potential stata journal {p_end}

{p 4 14 2}
{browse "https://github.com/kolesarm/RDHonest":Github page} of RDHonest {p_end}

{title:End}