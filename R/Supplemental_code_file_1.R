####Vector definitions####

#Kd.app = c("L-Aspartic acid-Mg",
#           "Glutathione-Mg",
#           "Glucose 1-P-Mg",
#           "ATP-Mg",
#           "AMP-Mg",
#           "L-Alanine-Mg",
#           "L-Asparagine-Mg",
#           "Pyruvic acid-Mg",
#           "L-Aspartic acid-Mn",
#           "Glutathione-Mn",
#           "Glucose 1-P-Mn",
#           "ATP-Mn",
#           "AMP-Mn",
#           "L-Alanine-Mn",
#           "L-Asparagine-Mn",
#           "Pyruvic acid-Mn",
#           "L-Aspartic acid-Zn",
#           "Glutathione-Zn",
#           "Glucose 1-P-Zn",
#           "ATP-Zn",
#           "AMP-Zn"
#           "L-Alanine-Zn",
#           "L-Asparagine-Zn",
#           "Pyruvic acid-Zn"
#           "L-Aspartic acid-Ca",
#           "Glutathione-Ca",
#           "Glucose 1-P-Ca"
#           "ATP-Ca",
#           "AMP-Ca",
#           "L-Alanine-Ca",
#           "L-Asparagine-Ca",
#           "Pyruvic acid-Ca")

#x =       c("Mg",                   #1
#            "Mn",                   #2
#            "Zn",                   #3
#            "Ca",                   #4
#            "L-Aspartic acid",      #5
#            "Glutathione",          #6
#            "Glucose 1-P",          #7
#            "ATP",                  #8
#            "AMP",                  #9
#            "L-Alanine",            #10
#            "L-Asparagine",         #11
#            "Pyruvic acid",         #12
#            "L-Aspartic acid_Mg",   #13
#            "Glutathione_Mg",       #14
#            "Glucose 1-P_Mg",       #15
#            "ATP_Mg",               #16
#            "AMP_Mg",               #17
#            "L-Alanine_Mg",         #18
#            "L-Asparagine_Mg",      #19
#            "Pyruvic acid_Mg",      #20
#            "L-Aspartic acid_Mn",   #21
#            "Glutathione_Mn",       #22
#            "Glucose 1-P_Mn",       #23
#            "ATP_Mn",               #24
#            "AMP_Mn",               #25
#            "L-Alanine_Mn",         #26
#            "L-Asparagine_Mn",      #27
#            "Pyruvic acid_Mn",      #28
#            "L-Aspartic acid_Zn",   #29
#            "Glutathione_Zn",       #30
#            "Glucose 1-P_Zn",       #31
#            "ATP_Zn",               #32
#            "AMP_Zn",               #33
#            "L-Alanine_Zn",         #34
#            "L-Asparagine_Zn",      #35
#            "Pyruvic acid_Zn",      #36
#            "L-Aspartic acid_Ca",   #37
#            "Glutathione_Ca",       #38
#            "Glucose 1-P_Ca",       #39
#            "ATP_Ca",               #40
#            "AMP_Ca",               #41
#            "L-Alanine_Ca",         #42
#            "L-Asparagine_Ca",      #43
#            "Pyruvic acid_Ca")      #44

####omega(x)####

omega = function(z) {
  x = z^2
  y = c()

  #Kd Mg
  y[1] = log10(x[1]*x[5]/(Kd.app[1]*x[13]))
  y[2] = log10(x[1]*x[6]/(Kd.app[2]*x[14]))
  y[3] = log10(x[1]*x[7]/(Kd.app[3]*x[15]))
  y[4] = log10(x[1]*x[8]/(Kd.app[4]*x[16]))
  y[5] = log10(x[1]*x[9]/(Kd.app[5]*x[17]))
  y[6] = log10(x[1]*x[10]/(Kd.app[6]*x[18]))
  y[7] = log10(x[1]*x[11]/(Kd.app[7]*x[19]))
  y[8] = log10(x[1]*x[12]/(Kd.app[8]*x[20]))

  #Kd Mn
  y[9]  = log10(x[2]*x[5]/(Kd.app[9]*x[21]))
  y[10] = log10(x[2]*x[6]/(Kd.app[10]*x[22]))
  y[11] = log10(x[2]*x[7]/(Kd.app[11]*x[23]))
  y[12] = log10(x[2]*x[8]/(Kd.app[12]*x[24]))
  y[13] = log10(x[2]*x[9]/(Kd.app[13]*x[25]))
  y[14] = log10(x[2]*x[10]/(Kd.app[14]*x[26]))
  y[15] = log10(x[2]*x[11]/(Kd.app[15]*x[27]))
  y[16] = log10(x[2]*x[12]/(Kd.app[16]*x[28]))


  #Kd Zn
  y[17] = log10(x[3]*x[5]/(Kd.app[17]*x[29]))
  y[18] = log10(x[3]*x[6]/(Kd.app[18]*x[30]))
  y[19] = log10(x[3]*x[7]/(Kd.app[19]*x[31]))
  y[20] = log10(x[3]*x[8]/(Kd.app[20]*x[32]))
  y[21] = log10(x[3]*x[9]/(Kd.app[21]*x[33]))
  y[22] = log10(x[3]*x[10]/(Kd.app[22]*x[34]))
  y[23] = log10(x[3]*x[11]/(Kd.app[23]*x[35]))
  y[24] = log10(x[3]*x[12]/(Kd.app[24]*x[36]))

  #Kd Ca
  y[25] = log10(x[4]*x[5]/(Kd.app[25]*x[37]))
  y[26] = log10(x[4]*x[6]/(Kd.app[26]*x[38]))
  y[27] = log10(x[4]*x[7]/(Kd.app[27]*x[39]))
  y[28] = log10(x[4]*x[8]/(Kd.app[28]*x[40]))
  y[29] = log10(x[4]*x[9]/(Kd.app[29]*x[41]))
  y[30] = log10(x[4]*x[10]/(Kd.app[30]*x[42]))
  y[31] = log10(x[4]*x[11]/(Kd.app[31]*x[43]))
  y[32] = log10(x[4]*x[12]/(Kd.app[32]*x[44]))

  #Total M2+
  y[33] = (x[1] + x[13] + x[14] + x[15] + x[16] + x[17] + x[18] + x[19] + x[20] - v.c[9])        #Mg
  y[34] = 10*(x[2] + x[21] + x[22] + x[23] + x[24] + x[25] + x[26] + x[27] + x[28] - v.c[10])    #Mn
  y[35] = 1000*(x[3] + x[29] + x[30] + x[31] + x[32] + x[33] + x[34] + x[35] + x[36] - v.c[11])  #Zn
  y[36] = 1000*(x[4] + x[37] + x[38] + x[39] + x[40] + x[41] + x[42] + x[43] + x[44] - v.c[12])  #Ca


  #Total metabolites
  y[37] = x[5]  + x[13] + x[21] + x[29] + x[37] - v.c[1]
  y[38] = x[6]  + x[14] + x[22] + x[30] + x[38] - v.c[2]
  y[39] = x[7]  + x[15] + x[23] + x[31] + x[39] - v.c[3]
  y[40] = x[8]  + x[16] + x[24] + x[32] + x[40] - v.c[4]
  y[41] = x[9]  + x[17] + x[25] + x[33] + x[41] - v.c[5]
  y[42] = x[10] + x[18] + x[26] + x[34] + x[42] - v.c[6]
  y[43] = x[11] + x[19] + x[27] + x[35] + x[43] - v.c[7]
  y[44] = x[12] + x[20] + x[28] + x[36] + x[44] - v.c[8]

  return(y)
}
