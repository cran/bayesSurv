\name{tandmobRoos}
\alias{tandmobRoos}
\title{Signal Tandmobiel data, version Roos}
\description{
  This is the dataset resulting from a longitudinal prospective dental
  study performed in Flanders (North of Belgium) in 1996 -- 2001. The
  cohort of 4\,468 randomly sampled children who attended the first year
  of the basic school at the beginning of the study was annualy dental
  examined by one of 16 trained dentists.
  The original dataset consists thus of at most 6 dental observations
  for each child.
  
  The dataset presented here contains mainly the
  information on the emergence and caries times summarized in the
  interval-censored observations. Some baseline covariates are also
  included here.

  For more detail on the design of the study see Vanobbergen et
  al. (2000).

  This is the version of the dataset used first by Leroy et al. (2005)
  and contains a subset of the \code{\link{tandmob2}}. Some children
  were removed to satisfy inclusion criteria given in Leroy et
  al. (2005). Additionally, left-censored emergence times of the
  permanent first molars are adjusted according to the eruption stage
  (see Leroy et al., 2005).
  
  This data set was then used in the analyses presented in
  Komárek and Lesaffre (2006, 2008).

  \bold{IMPORTANT NOTICE:} It is possible to use these data for your
  research work under the condition that each manuscript is first
  approved by \cr
  Prof. Emmanuel Lesaffre \cr
  Biostatistical Centre Katholieke Universiteit Leuven \cr
  Kapucijnenvoer 35 \cr
  B-3000 Leuven \cr
  Belgium \cr
  <\code{emmanuel.lesaffre@med.kuleuven.be}>  
}
\usage{data(tandmobRoos)}
\format{a~data frame with 4\,394 rows and the following
  variables
  \describe{
  \item{IDNR}{identification number of a child}
  \item{GENDER}{character \emph{boy} or \emph{girl}}
  \item{DOB}{character, date of birth in the format DDmmmYY}
  \item{PROVINCE}{factor, code of the province with
    \describe{
      \item{0 = }{Antwerpen}
      \item{1 = }{Vlaams Brabant}
      \item{2 = }{Limburg}
      \item{3 = }{Oost Vlaanderen}
      \item{4 = }{West Vlaanderen}
    }  
  }
  \item{EDUC}{factor, code of the educational system with
    \describe{
      \item{0 = }{Free}
      \item{1 = }{Community school}
      \item{2 = }{Province/council school}
    }  
  }
  \item{GIRL}{numeric, 0 = \emph{boy}, 1 = \emph{girl}}
  \item{EBEG.xx}{lower limit of the emergence (in years of age) of the
    permanent tooth xx. In contrast to \code{\link{tandmob2}}, the lower
    emergence limit for the permanent first molars that were originally
    left-censored, are adjusted according to the eruption stage (see
    Leroy, 2005 for more details).

    xx takes values 16, 26, 36, 46 (permanent first molars).
  }
  \item{EEND.xx}{upper limit of the emergence (in years of age) of the
    permanent tooth xx. \code{NA} if the emergence was right-censored.

    xx takes values as for the variable \code{EBEG.xx}.
  }
  \item{FBEG.xx}{lower limit for the caries time (in years of age, `F'
    stands for `failure') of the permanent tooth xx. \code{NA} if the
    caries time was left-censored.

    xx takes values as for the variable \code{EBEG.xx}.    
  }
  \item{FEND.xx}{upper limit for the caries time (in years of age, `F'
    stands for `failure') of the permanent tooth xx. \code{NA} if the
    caries time was right-censored.

    xx takes values as for the variable \code{EBEG.xx}.

    Unfortunately, for all teeth except 16, 26, 36 and 46 almost all
    the caries times are right-censored. For teeth 16, 26, 36, 46, the
    amount of right-censoring is only about 25\%.
  }
  \item{TOOTH.xx}{numeric, 0 or 1. Equal to 1 if the information
    concerning the permanent tooth was available, 0 if the permanent
    tooth xx was removed from the dataset by Kris.

    xx takes values 16, 26, 36, 46.

    These variables are almost useless for ordinary users.
  }
  \item{Txxd}{numeric, 0 or 1. It is equal to 1 if the deciduous tooth
    xx was decayed, 0 otherwise.

    xx takes values 54, 64, 74, 84 (deciduous first molars),
    55, 65, 75, 85 (deciduous second molars).
  }
  \item{Txxm}{numeric, 0 or 1. It is equal to 1 if the deciduous tooth
    xx was missing due to caries, 0 otherwise.

    xx takes values 54, 64, 74, 84 (deciduous first molars),
    55, 65, 75, 85 (deciduous second molars).
  }
  \item{Txxf}{numeric, 0 or 1. It is equal to 1 if the deciduous tooth
    xx was filled, 0 otherwise.

    xx takes values 54, 64, 74, 84 (deciduous first molars),
    55, 65, 75, 85 (deciduous second molars).
  }
  \item{Txxs}{numeric, 0 or 1. It is equal to 1 if the deciduous tooth
    xx was sound, 0 otherwise.

    xx takes values 54, 64, 74, 84 (deciduous first molars),
    55, 65, 75, 85 (deciduous second molars).
  }
  \item{SEAL.xx}{numeric, 0 or 1. It is equal to 1 if the permanent
    first molar xx was sealed in pits and fissures (a form of
    protection), 0 otherwise.

    xx takes values 16, 26, 36, 46 (permanent first molars).    
  }
  \item{FREQ.BR}{numeric, 0 or 1. It is equal to 1 if the child brushes
    daily the teeth, equal to 0 if he/she brushes less than once a day.
  }
  \item{PLAQUE.xx.1}{numeric, 0 or 1. It is equal to 1 if there was
    occlusal plaque in pits and fissures of the permanent tooth xx. It
    is equal to 0 if there was either no plaque present or the plaque
    was present on the total occlusal surface.

    xx takes values 16, 26, 36, 46 (permanent first molars).    
  }  
  \item{PLAQUE.xx.2}{numeric, 0 or 1. It is equal to 1 if there was
    occlusal plaque on the total occlusal surface of the permanent tooth xx. It
    is equal to 0 if there was either no plaque present or the plaque
    was present only in pits and fissures.

    xx takes values 16, 26, 36, 46 (permanent first molars).    
  }
  }
}
\source{Biostatistical Centre, Katholieke Universiteit Leuven,
  Kapucijnenvoer 35, 3000 Leuven, Belgium

  \code{URL:} \code{http://med.kuleuven.be/biostat/}

  Data collection was supported by Unilever, Belgium.
  The Signal Tandmobiel project comprises the following partners:
  D. Declerck (Dental School, Catholic University Leuven),
  L. Martens (Dental School, University Ghent),
  J. Vanobbergen (Oral Health Promotion and Prevention, Flemish Dental Association),
  P. Bottenberg (Dental School, University Brussels),
  E. Lesaffre (Biostatistical Centre, Catholic University Leuven),
  K. Hoppenbrouwers (Youth Health Department, Catholic University Leuven; Flemish Association for Youth Health Care).
}
\references{
  Komárek, A. and Lesaffre, E. (2008).
  Bayesian accelerated failure time model with multivariate doubly-interval-censored data
  and flexible distributional assumptions.
  \emph{Journal of the American Statistical Association}, \bold{103}, 523--533.

  Komárek, A. and Lesaffre, E. (2006).
  Bayesian semi-parametric accelerated failurew time model for paired
  doubly interval-censored data.
  \emph{Statistical Modelling}, \bold{6}, 3--22.

  Leroy, R., Bogaerts, K., Lesaffre, E., and Declerck, D. (2005).
  Effect of caries experience in primary molars on cavity formation
  in the adjacent permanent first molar.
  \emph{Caries Research}, \bold{39}, 342--349.
  
  Vanobbergen, J., Martens, L., Lesaffre, E., and Declerck, D. (2000).
  The Signal-Tandmobiel project -- a longitudinal intervention health promotion
  study in Flanders (Belgium): baseline and first year results.
  \emph{European Journal of Paediatric Dentistry}, \bold{2}, 87--96.
}  
\keyword{datasets}
