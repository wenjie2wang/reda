##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2025
##
## This file is part of the R package reda.
##
## The R package reda is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reda is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


##' Simulated Sample Dataset for Demonstration
##'
##' A simulated data frame with covariates named
##' \code{ID}, \code{time}, \code{event}, \code{group}, \code{x1},
##' and \code{gender}, where
##' \itemize{
##'     \item \code{ID}: Subjects identification;
##'     \item \code{time}: Event or censoring time;
##'     \item \code{event}: Event indicator, 1 = event, 0 = censored;
##'     \item \code{group}: Treatment group indicator;
##'     \item \code{x1}: Continuous variable.
##'     \item \code{gender}: Gender of subjects.
##' }
##'
##' @details
##' The sample dataset is originally simulated by the thinning
##' method developed by Lewis and Shedler (1979) and
##' further processed for a better demonstration purpose.
##' See Fu et al. (2016) for details also.
##'
##' @docType data
##' @name simuDat
##' @format A data frame with 500 rows and 6 variables.
##' @references
##' Lewis, P. A., & Shedler, G. S. (1979).
##' Simulation of nonhomogeneous Poisson processes by thinning.
##' \emph{Naval Research Logistics Quarterly}, 26(3), 403--413.
##'
##' Fu, H., Luo, J., & Qu, Y. (2016).
##' Hypoglycemic events analysis via recurrent time-to-event (HEART) models.
##' \emph{Journal Of Biopharmaceutical Statistics}, 26(2), 280--298.
NULL


##' Valve Seats Dataset
##'
##' Valve seats wear out in certain diesel engines, each with 16 valve seats.
##' The dataset served as an example of recurrence data in Nelson (1995),
##' which consists of valve-seat replacements on 41 engines in a fleet.
##' The covariates are named
##' \code{ID}, \code{Days}, and \code{No.}, where
##' \itemize{
##'     \item \code{ID}: The engine number;
##'     \item \code{Days}: Engine age in days;
##'     \item \code{No.}: Event indicator, '1' for a valve-seat replacement
##'         and, '0' for the censoring age of an engine.
##' }
##' @docType data
##' @name valveSeats
##' @format A data frame with 89 rows and 3 variables.
##' @references
##' Nelson, W. (1995), Confidence Limits for Recurrence
##' Data-Applied to Cost or Number of Product Repairs, \emph{Technometrics},
##' 37, 147--157.
NULL
