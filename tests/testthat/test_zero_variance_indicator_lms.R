devtools::load_all()

set.seed(123)
n <- 100  # number of observations

# Step 1: Simulate latent variables
# We define latent variables: psc, SWE, decisions, work_engagement, irritation, stress, intention
# Means and variances can be adjusted as needed

latent_means <- rep(0, 6)  # [psc, SWE, decisions_not_influenced_latent, work_engagement, irritation, stress]
latent_cov <- matrix(c(
    1, 0.5,  0.4,   0,     0,   0,
  0.5,   1,  0.4,   0,     0,   0,
  0.4, 0.4,    1,   0,     0,   0,
    0,   0,    0,   1,   0.3, 0.3,
    0,   0,    0, 0.3,     1, 0.4,
    0,   0,    0, 0.3,   0.4,   1
), nrow = 6, byrow = TRUE)
lvs <- c("psc", "SWE", "dec", "eng", "irr", "strs")
dimnames(latent_cov) <- list(lvs, lvs)

latent_data <- as.data.frame(mvtnorm::rmvnorm(n = n, mean = latent_means, sigma = latent_cov))
names(latent_data) <- lvs

# Step 2: Simulate indicators with loadings and residual variances
# We assume loadings of ~0.8 and residuals to account for remaining variance

make_indicators <- function(latent, prefix, items, loading = 0.8) {
  residual_sd <- sqrt(1 - loading^2)
  indicators <- matrix(NA, nrow = n, ncol = items)
  for (i in 1:items) {
    indicators[, i] <- loading * latent + rnorm(n, 0, residual_sd)
  }
  indicators <- as.data.frame(indicators)
  colnames(indicators) <- paste0(prefix, "_", 1:items)
  return(indicators)
}

# Measurement model
PSC <- make_indicators(latent_data$psc, "PSC", 4)
SWE <- make_indicators(latent_data$SWE, "Innovation", 4)
ENG <- make_indicators(latent_data$eng, "engagement", 3)
STRS <- make_indicators(latent_data$strs, "stress", 3)
IRR <- make_indicators(latent_data$irr, "irritation", 3)

# Single-indicator latent: decisions and intention (use dummy for decisions)
decisions_not_influenced <- 0.8 * latent_data$dec + rnorm(n, 0, sqrt(1 - 0.8^2))
intention <- ifelse(latent_data$eng + rnorm(n) < 0, 0, 1)  # dummy outcome ~ work_engagement

# Step 3: Create dummy variable names matching your SEM model
colnames(PSC) <- c("PSC_address_problems", "PSC_nobody_intentional", "PSC_dare_risk", "PSC_skills_appreciated")
colnames(SWE) <- c("Innovation_push_through", "Innovation_inspire", "Innovation_innovation_potential", "Innovation_trust")
colnames(ENG) <- c("fit_energetic", "inspiring_work", "completely_absorbed")
colnames(STRS) <- c("too_little_time_tasks", "tasks_inadequately_fulfilled", "pointless_tasks")
colnames(IRR) <- c("thinking_about_work_issues_at_home", "irritable", "feeling_like_bundle_of_nerves")

# Final dataset
dummy_data <- data.frame(
  PSC,
  SWE,
  ENG,
  STRS,
  IRR,
  decisions_not_influenced = decisions_not_influenced,
  intention_to_change_leadership_binary = intention
)


# Job Demands-Resources (JDR) model for SEM
JDR_Model_modsem <- '
# Measurement model
  psc =~ PSC_address_problems + PSC_nobody_intentional + PSC_dare_risk + PSC_skills_appreciated  # social resource
  SWE =~ Innovation_push_through + Innovation_inspire + Innovation_innovation_potential + Innovation_trust  # personal resource
  work_engagement =~ fit_energetic + inspiring_work + completely_absorbed
  stress =~ too_little_time_tasks + tasks_inadequately_fulfilled + pointless_tasks
  irritation =~ thinking_about_work_issues_at_home + irritable + feeling_like_bundle_of_nerves  # Strain
  decisions_not_influenced_latent =~ decisions_not_influenced  # Job Demands
  intention_to_change_leadership_binary_latent =~ intention_to_change_leadership_binary  # Outcome

# Covariances observed
  Innovation_push_through ~~ Innovation_inspire
  inspiring_work ~~ completely_absorbed
  tasks_inadequately_fulfilled ~~ pointless_tasks

# Covariances
  psc ~~ SWE + decisions_not_influenced_latent
  SWE ~~ decisions_not_influenced_latent
  irritation ~~ work_engagement

# Structural model
  work_engagement ~ psc + SWE + decisions_not_influenced_latent
  irritation ~ psc + SWE + decisions_not_influenced_latent
  stress ~ work_engagement + irritation
  intention_to_change_leadership_binary_latent ~ work_engagement

# Interactions
  work_engagement ~ psc:decisions_not_influenced_latent
  work_engagement ~ SWE:decisions_not_influenced_latent
'


testthat::expect_error(
  modsem(model = JDR_Model_modsem, data = dummy_data, method = "lms"),
  regex = ".*zero residual variance.*"
)

est_qml <- modsem(model = JDR_Model_modsem, data = dummy_data, method = "qml")
print(summary(est_qml))
print(diag(vcov(est_qml)))


JDR_Model_linear <- '
# Measurement model
  PSC =~ PSC_address_problems + PSC_nobody_intentional + PSC_dare_risk + PSC_skills_appreciated  # social resource
  SWE =~ Innovation_push_through + Innovation_inspire + Innovation_innovation_potential + Innovation_trust  # personal resource
  WE =~ fit_energetic + inspiring_work + completely_absorbed
  STRESS =~ too_little_time_tasks + tasks_inadequately_fulfilled + pointless_tasks
  IRRITATE =~ thinking_about_work_issues_at_home + irritable + feeling_like_bundle_of_nerves  # Strain
  NONINF =~ decisions_not_influenced  # Job Demands
  INT =~ intention_to_change_leadership_binary  # Outcome

# Structural model
  WE ~ PSC + SWE + NONINF
  IRRITATE ~ PSC + SWE + NONINF
  STRESS ~ WE + IRRITATE
  INT ~ WE
'

est_lin_lms <- modsem(JDR_Model_linear, dummy_data, "lms")
testthat::expect_equal(est_lin_lms$iterations, 2)
print(summary(est_lin_lms, H0 = FALSE))
