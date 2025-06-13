devtools::load_all()

tpb_uk <- "
# Outer Model (Based on Hagger et al., 2007)
 ATT =~ att3 + att2 + att1 + att4
 SN =~ sn4 + sn2 + sn3 + sn1
 PBC =~ pbc2 + pbc1 + pbc3 + pbc4
 INT =~ int2 + int1 + int3 + int4
 BEH =~ beh3 + beh2 + beh1 + beh4

# Inner Model (Based on Steinmetz et al., 2011)
 INT ~ ATT + SN + PBC
 BEH ~ INT + PBC
 BEH ~ INT:PBC
"

# Chronbach's Alpha ------------------------------------------------------------
corrected <- relcorr_single_item(syntax = tpb_uk, data = TPB_UK)
print(corrected)
syntax <- corrected$syntax
data   <- corrected$data

est_dca <- modsem(syntax, data = data, method = "dblcent")
est_lms <- modsem(syntax, data = data, method="lms", nodes=32)


if (FALSE) {
Data_C1 <- read.csv(file = '~/Downloads/Example_C1.csv')

# Specify the measurement model (Example.C1.Model.Measure) 
Example.C1.Model.Measure <- '
  JDemand =~ JobD1 + JobD2 + JobD3
  JResource =~ JobRes1+ JobRes2 + JobRes3 + JobRes4 + JobRes5 + JobRes6
  HomeSick =~ HomeS1 + HomeS2 + HomeS3 + HomeS4 + HomeS5 + HomeS6 + HomeS7 + HomeS8 + HomeS9 + HomeS10 +
              HomeS11 + HomeS12 + HomeS13 + HomeS14 + HomeS15 + HomeS16 + HomeS17 + HomeS18 + HomeS19 + HomeS20
  EStability =~ EmoStab1 + EmoStab2 + EmoStab3 + EmoStab4 + EmoStab5 + EmoStab6
  Openness =~ Open1 + Open2 + Open3 + Open4 + Open5 + Open6
  Performance =~ TaskP1 + TaskP2 + TaskP3

  # Structural model #
  Performance ~ a1*JResource + HomeSick + EStability + Openness
  Performance ~ z1*JResource:HomeSick + w1a*JResource:EStability + HomeSick:EStability
  Performance ~ w1b*JResource:Openness + HomeSick:Openness

  # Control variable #
  Performance ~ JDemand
'
corrected <- relcorr_single_item(Example.C1.Model.Measure, Data_C1)
print(corrected)
#> Average Variance Extracted:
#>   JDemand:        0.663
#>   JResource:      0.781
#>   HomeSick:       0.919
#>   EStability:     0.781
#>   Openness:       0.778
#>   Performance:    0.658
#> 
#> Construct Reliability:
#>   JDemand:        0.732
#>   JResource:      0.735
#>   HomeSick:       0.861
#>   EStability:     0.737
#>   Openness:       0.705
#>   Performance:    0.697
#> 
#> Generated Syntax:
#>   Performance ~ a1*JResource
#>   Performance ~ HomeSick
#>   Performance ~ EStability
#>   Performance ~ Openness
#>   Performance ~ z1*JResource:HomeSick
#>   Performance ~ w1a*JResource:EStability
#>   Performance ~ HomeSick:EStability
#>   Performance ~ w1b*JResource:Openness
#>   Performance ~ HomeSick:Openness
#>   Performance ~ JDemand
#>   JDemand =~ 1*composite_JDemand_
#>   JResource =~ 1*composite_JResource_
#>   HomeSick =~ 1*composite_HomeSick_
#>   EStability =~ 1*composite_EStability_
#>   Openness =~ 1*composite_Openness_
#>   Performance =~ 1*composite_Performance_
#>   composite_JDemand_ ~~ 0.124747077738605*composite_JDemand_
#>   composite_JResource_ ~~ 0.0769133338087597*composite_JResource_
#>   composite_HomeSick_ ~~ 0.0238661172257718*composite_HomeSick_
#>   composite_EStability_ ~~ 0.254605663714866*composite_EStability_
#>   composite_Openness_ ~~ 0.287668887959575*composite_Openness_
#>   composite_Performance_ ~~ 0.0705913680595234*composite_Performance_
#>   
#> Generated Items:
#>   'data.frame': 422 obs. of  6 variables:
#>     $ composite_JDemand_    : num  2.67 3.33 2.33 2.67 4 ...
#>     $ composite_JResource_  : num  2.67 3.5 3.17 3.5 3.33 ...
#>     $ composite_HomeSick_   : num  2.55 3 3.1 3 3.1 3.05 3.35 3.35 2.75 3.55 ...
#>     $ composite_EStability_ : num  5 5.5 4.5 5.17 4 ...
#>     $ composite_Openness_   : num  3.17 4.83 4.83 4.33 5.33 ...
#>     $ composite_Performance_: num  3 3.33 3 4 3.33 ...

est_dca <- modsem(corrected$syntax, data = corrected$data, method="dblcent")  # method = dblcent, rca, uca, ca, pind
est_lms <- modsem(corrected$syntax, data=corrected$data, method="lms")
}
