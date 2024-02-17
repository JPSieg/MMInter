#'Function that calculates the log10 of a metal ion affinity constant for Ionic strength = 0 using the Debye-Huckel equation
#'
#'@param log10K Molar affinity constant of a metabolite binding to a metal ion
#'@param I ionic strength
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param MX.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.0.calculator = function(log10K, I, M.charge, X.charge, XM.charge, A = 0.524){
  y.XM = 10^(0.1*(XM.charge^2)*I - (A*(XM.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.M = 10^(0.1*(M.charge^2)*I - (A*(M.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  log10K.0 = log10K - log10(y.XM/(y.M*y.X))
}

#'Function that calculates the log10 of a metal ion affinity constant for any Ionic strength  using the Debye-Huckel equation
#'
#'@param log10K.ref Molar affinity constant of a metabolite binding to a metal ion that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param XM.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.IS.calculator = function(log10K.ref, I.ref, I, M.charge, X.charge, XM.charge, A = 0.524){
  y.XM = 10^(0.1*(XM.charge^2)*I.ref - (A*(XM.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.M = 10^(0.1*(M.charge^2)*I.ref - (A*(M.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.X = 10^(0.1*(X.charge^2)*I.ref - (A*(X.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  log10K.0 = log10K.ref - log10(y.XM/(y.M*y.X))
  y.XM = 10^(0.1*(XM.charge^2)*I - (A*(XM.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.M = 10^(0.1*(M.charge^2)*I - (A*(M.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  log10K = log10K.0 + log10(y.XM/(y.M*y.X))
}

#'Function that calculates the pKa of any Ionic strength using the Debye-Huckel equation
#'
#'@param pKa.ref pKa for the reference that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param X.charge charge of the base
#'@param HX.charge Charge of acid
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
pKa.IS.calculator = function(pKa.ref, I.ref, I, X.charge, HX.charge, A = 0.524){
  y.HX = 10^(0.1*(HX.charge^2)*I.ref - (A*(HX.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.H = 10^(0.1*I.ref - (A*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.X = 10^(0.1*(X.charge^2)*I.ref - (A*(X.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  pKa.0 = pKa.ref - log10(y.HX/(y.H*y.X))
  print("pKa.0")
  print(pKa.0)
  y.HX = 10^(0.1*(HX.charge^2)*I - (A*(HX.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.H = 10^(0.1*I - (A*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  pKa = pKa.0 + log10(y.HX/(y.H*y.X))
  print("pKa")
  print(pKa)
}

