devtools::load_all()

set.seed(42)
n <- 100
# Generate dummy data
dummy_data <- data.frame(
  PSC_address_problems = sample(1:5, n, replace = TRUE),
  PSC_nobody_intentional = sample(1:5, n, replace = TRUE),
  PSC_dare_risk = sample(1:5, n, replace = TRUE),
  PSC_skills_appreciated = sample(1:5, n, replace = TRUE),
  Innovation_push_through = sample(1:5, n, replace = TRUE),
  Innovation_inspire = sample(1:5, n, replace = TRUE),
  Innovation_innovation_potential = sample(1:5, n, replace = TRUE),
  Innovation_trust = sample(1:5, n, replace = TRUE),
  fit_energetic = sample(1:5, n, replace = TRUE),
  inspiring_work = sample(1:5, n, replace = TRUE),
  completely_absorbed = sample(1:5, n, replace = TRUE),
  too_little_time_tasks = sample(1:5, n, replace = TRUE),
  tasks_inadequately_fulfilled = sample(1:5, n, replace = TRUE),
  pointless_tasks = sample(1:5, n, replace = TRUE),
  thinking_about_work_issues_at_home = sample(1:5, n, replace = TRUE),
  irritable = sample(1:5, n, replace = TRUE),
  feeling_like_bundle_of_nerves = sample(1:5, n, replace = TRUE),
  decisions_not_influenced = sample(1:5, n, replace = TRUE),
  intention_to_change_leadership_binary = sample(0:1, n, replace = TRUE)
)

# Save to .sav file
library(haven)

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

  # Optional (commented out):
  # irritation ~ psc:decisions_not_influenced_latent
  # irritation ~ SWE:decisions_not_influenced_latent
'


testthat::expect_error(
  modsem(model = JDR_Model_modsem, data = dummy_data, method = "lms"),
  regex = ".*zero residual variance.*"
)

testthat::expect_warning(
  modsem(model = JDR_Model_modsem, data = dummy_data, method = "qml"),
  regex = ".*Hessian matrix.*"
)
