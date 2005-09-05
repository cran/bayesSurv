\name{cgd}
\alias{cgd}
\title{Chronic Granulomatous Disease data}
\description{
  Dataset from Fleming and Harrington (1991).

  Data from a multicenter placebo-controlled randomized trial of gamma
  inferon in patients with chronic granulomatous disease (CGD). 128
  patients were randomized to either gamma inferon (\eqn{n=63}{n=63}) or
  placebo (\eqn{n=65}{n=65}). For each patient the times from study
  entry to initial and any recurrent serious infections are
  available. There is a minimum of one and a maximum of eight
  (recurrent) infection times per patient, with a total of 203 records.
}  
\usage{data(cgd)}
\format{a~data frame with 203 rows and 17 columns. There are the
  following variables contained in the data frame.
  \describe{
    \item{hospit}{code of the hospital}
    \item{ID}{case identification number}
    \item{RDT}{date randomization onto study (mmddyy)}
    \item{IDT}{date of onset of serious infection, or date the patient was taken off the study (mmddyy)}
    \item{trtmt}{treatment code, 1 = rIFN-gamma, 2 = placebo}
    \item{inherit}{pattern of inheritance, 1 = X-linked, 2 = autosomal recessive}
    \item{age}{age in years}
    \item{height}{height of the patient in cm}
    \item{weight}{weight of the patient in kg}
    \item{cortico}{using corticosteroids at time of study entry, 1 =
      yes, 2 = no}
    \item{prophy}{using prophylatic antibiotics at time of study entry,
      1 = yes, 2 = no}
    \item{gender}{1 = male, 2 = female}
    \item{hcat}{hospital category, 1 = US-NIH,  2 = US-other, 2 = Europe
    - Amsterdam, 4 = Europe - other}
    \item{T1}{elapsed time (in days) from randomization (from sequence = 1 record) to diagnosis of a serious infection
               or, if a censored observation, elapsed time from randomization to censoring date;
               computed as IDT - RDT (from sequence = 1 record)}
    \item{T2}{0, for sequence = 1 record, if sequence > 1, T2 = T1(from previous record) + 1}
    \item{event}{censoring indicator, 1 = non-censored observation, 2 = censored observation}
    \item{sequence}{sequence number, for each patient, the infection recods are in sequence number order}
  }  
}  
\source{Appendix D.2 of Fleming and Harrington (1991).
}  
\references{
  Fleming, T. R. and Harrington, D. P. (1991).
  \emph{Counting Processes and Survival Analysis}.
  New York: John Wiley and Sons.
}  
\keyword{datasets}
