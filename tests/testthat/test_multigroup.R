devtools::load_all()
model <- ' 
# Outer Model
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6

# Inner Model
  visual ~ c(v_t1, v_t2) * textual

# Intercepts
  visual ~ c(b01, b02) * 1
  x1 ~ 0 * 1 # fix, to avoid identification issues

# Interaction effects, and main effects of school
  textual_x_school := v_t2 - v_t1 # interaction
  main_school := b02 - b01
'

fit_lavaan <- lavaan::sem(model, 
                          data = lavaan::HolzingerSwineford1939, 
                          group = "school")

# using modsem 
fit_modsem <- modsem(model, 
                     data = lavaan::HolzingerSwineford1939, 
                     group = "school")

# the group after the | is the reference group
plot_interaction(x = "textual|Pasteur", z = "main_school",,
                 y = "visual|Pasteur", xz = "textual_x_school", vals_z = c(0, 1),
                 rescale = FALSE, model = fit_modsem)

lavaanEst <- lavaan::parameterEstimates(fit_lavaan)
lavaanEst[is.na(lavaanEst)] <- -999
modsemEst <- lavaan::parameterEstimates(fit_modsem$lavaan)
modsemEst[is.na(modsemEst)] <- -999
testthat::expect_equal(lavaanEst, modsemEst)

model <- ' 
# Outer Model
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9

# Inner Model
  visual ~ textual + speed + speed:textual
'

# using modsem 
fit_modsem <- modsem(model, 
                     data = lavaan::HolzingerSwineford1939, 
                     group = "school")

plot_interaction(x = "textual|Pasteur", z = "speed|Pasteur",,
                 y = "visual|Pasteur", xz = "speed:textual|Pasteur", vals_z = c(-1, 1),
                 model = fit_modsem)
