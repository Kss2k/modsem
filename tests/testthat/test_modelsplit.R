devtools::load_all()

tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + PBC + SN
  BEH ~ INT + PBC 
  BEH ~ INT:PBC  
'

parTable <- modsemify(tpb)
splitParTable(parTable)
est_qml <- modsem(tpb, TPB, method = "qml")
summary(est_qml)


model <- '
  # behavior
  behavior =~ var82 + 0.991*var83
  
  # intention
  intent =~ var77 + var78 + var79 + var80
  
  # motivation
  motivat =~ var73 + var74 + var75 + var76
  var73 ~~ var74
  
  # attitude
  attitude =~ var33 + var34 + var35 + var36 + var37
  var34 ~~ var37
  
  # procurement goal beliefs (e*w)
  EWprocur =~ v7_v12 + v8_v13 + v9_v14 + v10_v15 + v11_v16
  v7_v12 ~~ v8_v13
  v10_v15 ~~ v11_v16
  
  # other attitudional beliefs (e*w) - environment
  EWEoAttBe =~ v38_v46 + v40_v48 + v42_v50 + v43_v51
  
  # other attitudional beliefs (e*w) - shoping
  EWSoAttBe =~ v39_v47 + v41_v49 + v44_v52 + v45_v53
  
  # subjective norm
  sub_norm =~ var54 + var55 + var56 + var57
  var54 ~~ var55
  
  # approval goal beliefs (e*w) - social environment
  EWEapprov =~ v17_v25 + v18_v26 + v19_v27
  
  # approval goal beliefs (e*w) - others
  EWOapprov =~ v22_v30 + v23_v31 + v24_v32
  
  # perceived behavioral controls
  pbc =~ var58 + var59 + var60 + var61 + var62
  var61 ~~ var62
  var58 ~~ var61
  
  # control beliefs (e*w)
  EWcb =~ v63_v68 + v64_v69 + v65_v70 + v66_v71 + v67_v72
  v63_v68 ~~ v64_v69

  # structural part
  behavior  ~ intent + pbc                       + intent:pbc
  intent    ~ motivat + pbc                      + motivat:pbc
  motivat   ~ attitude + sub_norm   + EWprocur   + attitude:EWprocur +
    EWEapprov  + sub_norm:EWEapprov + EWOapprov  + sub_norm:EWOapprov
  attitude  ~ EWEoAttBe + EWSoAttBe + EWprocur   + EWEoAttBe:EWprocur +
                                                  EWSoAttBe:EWprocur
  sub_norm  ~ EWEapprov + EWOapprov
  pbc       ~ EWcb
  EWEoAttBe ~ EWprocur
  EWSoAttBe ~ EWprocur
  
  EWEoAttBe ~~ EWSoAttBe + EWEapprov + EWOapprov + EWcb
  EWSoAttBe ~~ EWEapprov + EWOapprov + EWcb
    
  pbc ~~ sub_norm
'

splitParTable(modsemify(model))

df <- read.table("sketches/sam_interaction_paper/inter_lavaan.dat")
names(df) <- c("var33", "var34", "var35", "var36", "var37", "var54", "var55",
                 "var56", "var57", "var58", "var59", "var60", "var61", "var62",
                 "var73", "var74", "var75", "var76", "var77", "var78", "var79",
                 "var80", "var82", "var83", "v7_v12", "v8_v13", "v9_v14",
                 "v10_v15", "v11_v16", "v38_v46", "v40_v48", "v42_v50",
                 "v43_v51", "v39_v47", "v41_v49", "v44_v52", "v45_v53",
                 "v17_v25", "v18_v26", "v19_v27", "v22_v30", "v23_v31",
                 "v24_v32", "v63_v68", "v64_v69", "v65_v70", "v66_v71",
                 "v67_v72")
df[df == -999] <- as.numeric(NA)
dim(df)
names(df)

est_qml <- modsem(model, 
 
