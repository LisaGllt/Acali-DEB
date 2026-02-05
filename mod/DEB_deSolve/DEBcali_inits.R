initParms <- function(newParms = NULL) {
  parms <- c(
    muOM = 0,
    rOM_ClxHorse = 0,
    Fm = 0,
    kapX = 0,
    pAm = 0,
    v = 0,
    kap = 0,
    pM = 0,
    pT = 0,
    Eg = 0,
    Shape = 0,
    w = 0,
    kJ = 0,
    Ehb = 0,
    Ehp = 0,
    L_coc = 0,
    E_coc = 0,
    kapR = 0,
    TAH = 0,
    TH = 0,
    TA = 0,
    Tref = 0,
    Winit = 0,
    Texp = 0,
    WeightSoilCosm = 0
  )

  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
    parms[names(newParms)] <- newParms
  }

  parms <- within(as.list(parms), {
  })
  out <- .C("getParms",  as.double(parms),
            out=double(length(parms)),
            as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

Outputs <- c(
    "Energy",
    "Maturity",
    "Reproduction",
    "Weight",
    "Organic_matter"
)

initStates <- function(parms, newStates = NULL) {
  Y <- c(
    E = 0.0,
    L = 0.0,
    Eh = 0.0,
    R = 0.0,
    OM_soil = 0.0,
    OM_horse = 0.0,
    W = 0.0,
    OM = 0.0
  )

  Y <- within(c(as.list(parms),as.list(Y)), {
    Y["E"] <- Em 
    Y["L"] <- pow ( ( Winit / ( 1 +1 * w ) ) , ( 1.0 / 3.0 ) ) 
    Y["Eh"] <- Ehb 
    Y["R"] <- 1e-6 

  })$Y

  if (!is.null(newStates)) {
    if (!all(names(newStates) %in% c(names(Y)))) {
      stop("illegal state variable name in newStates")
    }
    Y[names(newStates)] <- newStates
  }

.C("initState", as.double(Y));
Y
}
