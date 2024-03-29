% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AssessMe.R
\name{accuracy_hpc}
\alias{accuracy_hpc}
\title{Test the robustness of clusters}
\usage{
accuracy_hpc(generate, crossvali = iterat, ntree = 200, loreg = F)
}
\arguments{
\item{generate}{output of generate_input_accuracy_hpc function.}

\item{crossvali}{numeric index to iterate through indixes of n-fold sampled cells for training data.}

\item{ntree}{number of trees to grow for random forest based reclassification.}

\item{loreg}{logical. If \code{T}, perform n-fold crossvalidation based on multionomial logistic regression. Default = \code{F}.}
}
\description{
accuracy function to be executed on a high performing computing cluster system.
}
\examples{
give_accuracy <- accuracy(generate = tes_generate, crossvali = iterat, ntree = 100)
}
