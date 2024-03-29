% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matching_criterions.R
\name{no_additional_matching}
\alias{no_additional_matching}
\alias{additional_wcvp_matching}
\title{Matching criterions}
\usage{
no_additional_matching(enrich_database_extract, message = "")

additional_wcvp_matching(enrich_database_extract, message = "")
}
\arguments{
\item{enrich_database_extract}{extract of enrichment database (often records with identical taxonomic names).}

\item{message}{Matching message.}
}
\value{
Always returns a list of 2 with:
\itemize{
\item \verb{$row} The rows of \code{enrich_database_extract} corresponding to the best match/s according to the criterion.
\item \verb{$message} A combination of the input message the method used to determine a match.
}
}
\description{
Custom matching criterions used within \code{match_all_issue()} and \code{get_match_from_multiple()}.
}
\details{
Within BGSmartR we have two in-built custom matching criterions.
\itemize{
\item \code{no_additional_matching()} is used when no custom matching wants to be performed.
\item \code{additional_wcvp_matching()} is the default for additional matching to the World Checklist of Vascular Plants (WCVP) database. This method includes looking ahead to see if potential matches lead to the same accepted plant name, and matching on the taxon_status of the records (Accepted or Synonym).
}

These can be used within \code{\link[=match_all_issue]{match_all_issue()}} and \code{\link[=get_match_from_multiple]{get_match_from_multiple()}} via the input \code{matching_criterion}. See vignette \verb{Method of Matching taxonomic records} for further details on custom matching and how to create new matching criterions.
}
