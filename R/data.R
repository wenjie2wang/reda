################################################################################
##
##   R package heart by Haoda Fu, Jun Yan, and Wenjie Wang
##   Copyright (C) 2015
##
##   This file is part of the R package heart.
##
##   The R package heart is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package heart is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package heart. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#' Simulated Dataset for Demonstration
#'
#' A data frame with covariates named 
#' 'ID', 'time', 'event', 'group' and 'X1'.
#' 
#' The dataset is simulated by the thinning method developed 
#' by \emph{Lewis and Shedler (1979)}. 
#' See \emph{Fu et al. (2014)} for more details.
#' @docType data
#' @name simuDat
#' @format data frame
#' @references 
#' Lewis, P. and Shedler, G. (1979), 
#' “Simulation of nonhomogeneous Poisson processes by thinning,” 
#' \emph{Naval Research Logistics Quarterly}, 26, 403–413.
#' 
#' Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
#' "Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models," 
#' \emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
#' @importFrom utils head 
#' @examples 
#' data(simuDat)
#' head(simuDat)
#' str(simuDat)
NULL
