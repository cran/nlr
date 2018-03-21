\name{ntp}
\alias{ntp}
\docType{data}
\title{ntp data
%%   ~~ data name/kind ... ~~
}
\description{Cromium concentration in blood and kidney of Mouse, Rat, and Guinea pig.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(ntp)}
\format{
  The format is:data.frame
 chr "ntp"
 \itemize{
  \item dm.k, cm.k, Dose (x) and Chromium Concentration (y) in Mouse kidney. 
  \item  dr.b, cr.b, Dose (x) and Chromium Concentration (y) Rat blood.
  \item dr.k, cr.k, Dose (x) and Chromium Concentration (y) Rat kidney.
  \item dp.b, cp.b, Dose (x) and Chromium Concentration (y) Guinea pig blood.
  \item dp.k, cp.k, Dose (x) and Chromium Concentration (y) in Guinea pig kidney.
 }
}
\details{Lim et al. (2010) illustrated the Weighted M-Estimate (WME) methodology with real data from National Toxicology study Program (NTP 2007. pp. 11-12). NTP (2007) conducted 3-month and 2-year studies, where rodents were exposed to CrVI administered in drinking water as sodium dichromate dihydrate.The dose concentrations were 0, 2.87, 8.62, 28.7, 86.2, 287, and 862mg sodium dichromate dihyrate/L (to yield 0, 1, 3, 10, 30, 100, and 300 mg chromium/L). When animals were sacrificed, total chromium concentrations in blood, kidneys, and femurs were measured.
Lim et al. (2010) proposed the hill model for the data 
\deqn{
\theta_0 + \frac{\theta_1 x_i^{\theta_2}}{\theta_3^{\theta_2}+x_i^{\theta_2}}
}
with heterogeneous standard deviation 
\deqn{
\tau_0+\frac{\tau_1}{1+e^{-\tau_2 x_i}} 
}
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{National Toxicology Program, 2007. NTP Toxicity Studies of Sodium Dichromate Dihydrate (CAS No. 7789-12-0) Administered in Drinking Water to Male and Female F344/N Rats and B6C3F1 Mice and Male BALB/c and am3-C57BL/6 Mice. Toxicity Report Series 72, 1-G4, U.S. Department of Health and Human Services, Public Health Service, National Institutes of Health, RTP, North Carolina, U.S.A.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Lim, C., Sen, P. K., Peddada, S. D. (2010). Statistical inference in nonlinear regression under heteroscedasticity. Sankhya B 72:202-218.
}
\examples{
data(ntp)
## maybe strnt ntp}
\keyword{datasets}
