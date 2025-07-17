#' oneInt
#'
#' @name oneInt
#' @docType data
#' @description A simulated dataset based on the elementary interaction model.
#' @examples
#'
#' m1 <- "
#' # Outer Model
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#' # Inner Model
#'   Y ~ X + Z + X:Z
#' "
#'
#' est <- modsem(m1, data = oneInt)
#' summary(est)
NULL


#' TPB
#'
#' @name TPB
#' @docType data
#' @description A simulated dataset based on the Theory of Planned Behaviour
#' @examples
#'
#' tpb <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC + INT:PBC
#' "
#'
#' est <- modsem(tpb, data = TPB)
#' summary(est)
NULL


#' TPB_1SO
#'
#' @name TPB_1SO
#' @docType data
#' @description A simulated dataset based on the Theory of Planned Behaviour,
#' where INT is a higher order construct of ATT, SN, and PBC.
#' @examples
#' tpb <- '
#'   # First order constructs
#'   ATT =~ att1 + att2 + att3
#'   SN  =~ sn1 + sn2 + sn3
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   BEH =~ b1 + b2
#' 
#'   # Higher order constructs
#'   INT =~ ATT + PBC + SN
#' 
#'   # Higher order interaction
#'   INTxPBC =~ ATT:PBC + SN:PBC + PBC:PBC
#'   
#'   # Structural model
#'   BEH ~ PBC + INT + INTxPBC
#' '
#' 
#' \dontrun{
#' est_ca <- modsem(tpb, data = TPB_1SO, method = "ca")
#' summary(est_ca)
#' 
#' est_dblcent  <- modsem(tpb, data = TPB_1SO, method = "dblcent")
#' summary(est_dblcent)
#' }
NULL


#' TPB_2SO
#'
#' @name TPB_2SO
#' @docType data
#' @description A simulated dataset based on the Theory of Planned Behaviour,
#' where INT is a higher order construct of ATT and SN, and PBC is a higher order
#' construct of PC and PB.
#' @examples
#' tpb <- '
#'   # First order constructs
#'   ATT =~ att1 + att2 + att3
#'   SN  =~ sn1 + sn2 + sn3
#'   PB =~ pb1 + pb2 + pb3
#'   PC =~ pc1 + pc2 + pc3
#'   BEH =~ b1 + b2
#' 
#'   # Higher order constructs
#'   INT =~ ATT + SN
#'   PBC =~ PC + PB
#' 
#'   # Higher order interaction
#'   INTxPBC =~ ATT:PC + ATT:PB + SN:PC + SN:PB
#'   
#'   # Structural model
#'   BEH ~ PBC + INT + INTxPBC
#' '
#' 
#' \dontrun{
#' est <- modsem(tpb, data = TPB_2SO, method = "ca")
#' summary(est)
#' }
NULL


#' TPB_UK
#'
#' @name TPB_UK
#' @docType data
#' @description A dataset based on the Theory of Planned Behaviour from a
#' UK sample. 4 variables with high communality were selected for each
#' latent variable (ATT, SN, PBC, INT, BEH), from two time points (t1 and t2).
#'
#' @source
#' Gathered from a replication study by Hagger et al. (2023).
#' 
#' Obtained from \url{https://doi.org/10.23668/psycharchives.12187}
#'
#' @examples
#'
#' tpb_uk <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'  ATT =~ att3 + att2 + att1 + att4
#'  SN =~ sn4 + sn2 + sn3 + sn1
#'  PBC =~ pbc2 + pbc1 + pbc3 + pbc4
#'  INT =~ int2 + int1 + int3 + int4
#'  BEH =~ beh3 + beh2 + beh1 + beh4
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'  # Causal Relationsships
#'  INT ~ ATT + SN + PBC
#'  BEH ~ INT + PBC
#'  BEH ~ INT:PBC
#' "
#'
#' est <- modsem(tpb_uk, data = TPB_UK)
#' summary(est)
NULL


#' Jordan subset of PISA 2006 data
#'
#' @name jordan
#' @docType data
#' @description The data stem from the large-scale assessment study PISA 2006
#' (Organisation for Economic Co-Operation and Development, 2009) where
#' competencies of 15-year-old students in reading, mathematics, and science
#' are assessed using nationally representative samples in 3-year cycles.
#' In this example, data from the student background questionnaire from the
#' Jordan sample of PISA 2006 were used. Only data of students with complete
#' responses to all 15 items (N = 6,038) were considered.
#'
#' @format
#' A data frame of fifteen variables and 6,038 observations:
#'
#' \describe{
#' \item{enjoy1}{indicator for enjoyment of science, item ST16Q01: I generally have fun when I am learning <broad science> topics.}
#' \item{enjoy2}{indicator for enjoyment of science, item ST16Q02: I like reading about <broad science>.}
#' \item{enjoy3}{indicator for enjoyment of science, item ST16Q03: I am happy doing <broad science> problems.}
#' \item{enjoy4}{indicator for enjoyment of science, item ST16Q04: I enjoy acquiring new knowledge in <broad science>.}
#' \item{enjoy5}{indicator for enjoyment of science, item ST16Q05: I am interested in learning about <broad science>.}
#' \item{academic1}{indicator for academic self-concept in science, item ST37Q01: I can easily understand new ideas in <school science>.}
#' \item{academic2}{indicator for academic self-concept in science, item ST37Q02: Learning advanced <school science> topics would be easy for me.}
#' \item{academic3}{indicator for academic self-concept in science, item ST37Q03: I can usually give good answers to <test questions> on <school science> topics.}
#' \item{academic4}{indicator for academic self-concept in science, item ST37Q04: I learn <school science> topics quickly.}
#' \item{academic5}{indicator for academic self-concept in science, item ST37Q05: <School science> topics are easy for me.}
#' \item{academic6}{indicator for academic self-concept in science, item ST37Q06: When I am being taught <school science>, I can understand the concepts very well.}
#' \item{career1}{indicator for career aspirations in science, item ST29Q01: I would like to work in a career involving <broad science>.}
#' \item{career2}{indicator for career aspirations in science, item ST29Q02: I would like to study <broad science> after <secondary school>.}
#' \item{career3}{indicator for career aspirations in science, item ST29Q03: I would like to spend my life doing advanced <broad science>.}
#' \item{career4}{indicator for career aspirations in science, item ST29Q04: I would like to work on <broad science> projects as an adult.}
#' }
#'
#' @source
#' This version of the dataset, as well as the description was gathered from the
#' documentation of the 'nlsem' package (https://cran.r-project.org/package=nlsem),
#' where the only difference is that the names of the variables were changed
#'
#' Originally the dataset was gathered by the Organisation for Economic Co-Operation and Development (2009).
#' Pisa 2006: Science competencies for tomorrow's world (Tech. Rep.).
#' Paris, France. Obtained from: https://www.oecd.org/pisa/pisaproducts/database-pisa2006.htm
#'
#' @examples
#' \dontrun{
#' m1 <- "
#'   ENJ =~ enjoy1 + enjoy2 + enjoy3 + enjoy4 + enjoy5
#'   CAREER =~ career1 + career2 + career3 + career4
#'   SC =~ academic1 + academic2 + academic3 + academic4 + academic5 + academic6
#'   CAREER ~ ENJ + SC + ENJ:ENJ + SC:SC + ENJ:SC
#' "
#'
#' est <- modsem(m1, data = jordan, method = "qml")
#' summary(est)
#' }
NULL
