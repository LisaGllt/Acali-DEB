#source("SA_Sobol_Trout_ChronicScenario.R")

setwd("/ccc/cont004/home/ineris/beaudoir/Trout_PBPK/ChronicScenario/Global")
source('/ccc/cont004/home/ineris/beaudoir/Trout_PBPK/RainbowTrout_PBPK.r')
source('/ccc/cont004/dsku/zinc/home/user/ineris/beaudoir/Trout_PBPK/lowry_plots.R')

cheminWork = "/ccc/work/cont004/ineris/beaudoir/Trout_PBPK/ChronicScenario/Global/"


# package
library("sensitivity")
library("parallel")
library("gplots")
library("reshape2")
library("ggplot2")

#=========================================================================================================================
# (1)                                               INPUTS
#=========================================================================================================================
mean.f = function(x){ mean(x,na.rm=T) }

# Experiment characteristics
nsim       = 30000
processors = 25
nboot      = 1000 # nb replicat
NOut = 3          # number of outputs

# Parameter names to take into account

  Names  = c(
    # Physiological parameters
    "Bw_i",   "Bw_ref",  "v",  "g",  "KM",  "EHm",  "EHb",  
    "shape",  "a_BW_L",  "b_BW_L" ,
    
    # Environmental condition
    "Temperature", "TA",  "f",  "OEE",  "S",  
  
    "Qb_ref",    "Vo2_ref",
    "sc_blood", # volume scaling factor : fraction of BW (%)
    "sc_gonads",  "sc_brain",   "sc_liver",   "sc_fat",   "sc_skin",
    "sc_GIT",     "sc_kidney",  "sc_rp",  
    
    "frac_gonads",  # Fraction of arterial blood flow
    "frac_brain",  "frac_liver", "frac_fat", "frac_skin", "frac_GIT" , 
    "frac_kidney", "frac_rp",    "a_Fpp",    "a_Fs" , 
    
    "log_Kow", "a_PC" , "b_PC" , 
    
    "water_liver", # relative value (%)
    "water_brain" , "water_gonads", "water_fat", "water_skin", "water_GIT",
    "water_kidney", "water_rp",     "water_pp",
    
    "lipids_liver" , # relative value (%)
    "lipids_brain",    "lipids_gonads" ,    "lipids_fat" ,    "lipids_skin" ,
    "lipids_GIT" ,    "lipids_kidney" ,    "lipids_rp",     "lipids_pp" , 
    
    "Cl", "Ke_bile",  "Ke_urine" )

  NP    = length(Names)            # nbr param

# Model information
  TimeExpo = 1000
  TimeMesure = c(0.99*TimeExpo)
  ntemp = length( TimeMesure )

# Storage Objects
  SARes = vector("list", 4) 
  names(SARes) = c( "Plan experience","Sortie","First order indices", "Total indices") 
  fichierSim = paste0(cheminWork, "PlanExp.txt")
  
#==========================================================================================================================
# (2)                                    Paramter distributions
#==========================================================================================================================       

# Distribution des parametres
  Distribution = rep("Unif", NP)
  
 Medians = c( 
  
    Bw_i   = 500, Bw_ref = 500, v = 0.38,   g = 0.196,  KM = 0.027, 
    EHm    = 452, EHb    = 52,  shape  = 0.113, a_BW_L = 0.0051, b_BW_L = 3.3513  , 

    Temperature  = 8,   TA= 6930,   f = 0.9, 
  
    Qb_ref  = 28.8,  Vo2_ref = 1.512,  OEE  = 0.71,   S = 0.9,  
    
    sc_blood  = 0.0527, sc_gonads = 0.199 ,  sc_brain  = 0.0049,  sc_liver  = 0.0165,   
    sc_fat    = 0.095,  sc_skin   = 0.10  ,  sc_GIT    = 0.099 ,  sc_kidney = 0.0076,
    sc_rp     = 0.0189,
    
    frac_gonads = 0.0138,  frac_brain  = 0.0002, frac_liver  = 0.0158,  frac_fat = 0.07, 
    frac_skin   = 0.0570  , frac_GIT    = 0.176, frac_kidney = 0.0817,   frac_rp = 0.0952   , 
    a_Fpp = 0.4,    a_Fs  = 0.1,

    log_Kow   = 4.96 , a_PC      = 0.78 ,  b_PC      = 0.82 ,
    
    water_liver = 0.746   ,  water_brain = 0.76    , water_gonads= 0.65,
    water_fat   = 0.05    , water_skin  = 0.667   ,  water_GIT   = 0.60,  
    water_kidney= 0.789   ,  water_rp    = 0.68    , water_pp    = 0.769, 
    
    lipids_liver = 0.045 , lipids_brain = 0.076 , lipids_gonads= 0.076 ,
    lipids_fat   = 0.942 , lipids_skin  = 0.029 , lipids_GIT   = 0.497 ,
    lipids_kidney= 0.052 , lipids_rp    = 0.032 , lipids_pp    = 0.03  , 

    Cl = 1, Ke_bile  = 1,  Ke_urine = 1   )  

  Lower=  Medians - Medians*0.1
  Upper=  Medians + Medians*0.1
 
  #Special cases
  Lower["log_Kow"]=0
  Upper["log_Kow"]=6
  
  Lower["Cl"]=0
  Upper["Cl"]=6192
  
  Lower["Ke_bile"]=0
  Upper["Ke_bile"]=5
  
  Lower["Ke_urine"]=0
  Upper["Ke_urine"]=50
  
#==========================================================================================================================
# (3)                                         Plan d'experience
#==========================================================================================================================

	# the first random sample.
	X1 = matrix(NA, nrow=nsim, ncol=NP) 

	for(i in 1:NP){
	  if( Distribution[i] == "Unif" ) { X1[,i] = runif(nsim, min = Lower[i], max= Upper[i]) }
	} 

	# the second random sample.
	X2 = matrix(NA, nrow=nsim, ncol=NP)
	for(i in 1:NP){
	  if( Distribution[i] == "Unif" ) { X2[,i] = runif(nsim, min= Lower[i], max= Upper[i]) }
	}

	colnames(X1)=colnames(X2)=Names

	#simulation experimental design
	sa = soboljansen(model = NULL, X1, X2, nboot =  nboot, conf = 0.95)
	SARes[[1]] = sa$X

	#write.table(sa$X ,file=fichierSim,sep="\t",quote=F) 

	# simulation experiement splitting (parallel computing)
	nexp = nrow(sa$X)
	Group.Sim =split(1:nexp, seq(1,nexp,l=10) )


  save(Group.Sim  , file = paste0(cheminWork,"GroupSim.RData") ) 

#=========================================================================================================================
# (4)                                               parallelisation
#=========================================================================================================================

  # (2.1) the initialization function
  #===============================================  
  prepro <- function(dummy) {
    source('/ccc/cont004/home/ineris/beaudoir/Trout_PBPK/RainbowTrout_PBPK.r')

	  library("sensitivity") }

  # (2.3) Simulation function
  #===============================================
  Simulation =function(NS, Parameter, TimeM, Texpo ){
    
	Res=PBPK.RT.F(

	  Bw_i   = as.numeric( Parameter[NS,"Bw_i"]),  
	  Bw_ref = as.numeric( Parameter[NS,"Bw_ref"])  , 
	  v      = as.numeric( Parameter[NS, "v"]),  
	  g      = as.numeric( Parameter[NS, "g"]), 
	  KM     = as.numeric( Parameter[NS, "KM"]),  
	  EHm    = as.numeric( Parameter[NS, "EHm"]), 
	  EHb    = as.numeric( Parameter[NS, "EHb"]),  
	  shape  = as.numeric( Parameter[NS, "shape"]),
	  a_BW_L = as.numeric( Parameter[NS, "a_BW_L"]),  
	  b_BW_L = as.numeric( Parameter[NS, "b_BW_L"]),  

	  # Environmental condition
	  Temperature  = as.numeric( Parameter[NS,"Temperature"])   ,   
	  TA           = as.numeric( Parameter[NS,"TA"]) ,   
	  f            = as.numeric( Parameter[NS,"f"])   ,  
	  
	  Qb_ref  = as.numeric( Parameter[NS,"Qb_ref"]),
	  Vo2_ref = as.numeric( Parameter[NS,"Vo2_ref"]),  
	  OEE     = as.numeric( Parameter[NS,"OEE"]),  
	  S       = as.numeric( Parameter[NS,"S"]), 
	  
	  sc_blood  = as.numeric( Parameter[NS,"sc_blood"]),  
	  sc_gonads = as.numeric( Parameter[NS,"sc_gonads"]),  
	  sc_brain  = as.numeric( Parameter[NS,"sc_brain"]),   
	  sc_liver  = as.numeric( Parameter[NS,"sc_liver"]),    
	  sc_fat    = as.numeric( Parameter[NS,"sc_fat"]),  
	  sc_skin   = as.numeric( Parameter[NS,"sc_skin"]),  
	  sc_GIT    = as.numeric( Parameter[NS,"sc_GIT"]),   
	  sc_kidney = as.numeric( Parameter[NS,"sc_kidney"]), 
	  sc_rp     = as.numeric( Parameter[NS,"sc_rp"]), 
	  
	  frac_gonads = as.numeric( Parameter[NS,"frac_gonads"]),   
	  frac_brain  = as.numeric( Parameter[NS,"frac_brain"]), 
	  frac_liver  = as.numeric( Parameter[NS,"frac_liver"]),   
	  frac_fat    = as.numeric( Parameter[NS,"frac_fat"]),  
	  frac_skin   = as.numeric( Parameter[NS,"frac_skin"]), 
	  frac_GIT    = as.numeric( Parameter[NS,"frac_GIT"]), 
	  frac_kidney = as.numeric( Parameter[NS,"frac_kidney"]), 
	  frac_rp     = as.numeric( Parameter[NS,"frac_rp"]), 
	  a_Fpp       = as.numeric( Parameter[NS,"a_Fpp"]), 
	  a_Fs        = as.numeric( Parameter[NS,"a_Fs"]), 
	  
	  # Chemical parameters
	  log_Kow   = as.numeric( Parameter[NS,"log_Kow"]),
	  a_PC      = as.numeric( Parameter[NS,"a_PC"]), 
	  b_PC      = as.numeric( Parameter[NS,"b_PC"]),
	  
	  water_liver = as.numeric( Parameter[NS,"water_liver"]),  
	  water_brain = as.numeric( Parameter[NS,"water_brain"]),
	  water_gonads= as.numeric( Parameter[NS,"water_gonads"]),
	  water_fat   = as.numeric( Parameter[NS,"water_fat"]), 
	  water_skin  = as.numeric( Parameter[NS,"water_skin"]),  
	  water_GIT   = as.numeric( Parameter[NS,"water_GIT"]),  
	  water_kidney= as.numeric( Parameter[NS,"water_kidney"]),  
	  water_rp    = as.numeric( Parameter[NS,"water_rp"]), 
	  water_pp    = as.numeric( Parameter[NS,"water_pp"]),
	  
	  lipids_liver = as.numeric( Parameter[NS,"lipids_liver"]),
	  lipids_brain = as.numeric( Parameter[NS,"lipids_brain"]),
	  lipids_gonads= as.numeric( Parameter[NS,"lipids_gonads"]),
	  lipids_fat   = as.numeric( Parameter[NS,"lipids_fat"]), 
	  lipids_skin  = as.numeric( Parameter[NS,"lipids_skin"]), 
	  lipids_GIT   = as.numeric( Parameter[NS,"lipids_GIT"]),
	  lipids_kidney= as.numeric( Parameter[NS,"lipids_kidney"]), 
	  lipids_rp    = as.numeric( Parameter[NS,"lipids_rp"]),
	  lipids_pp    = as.numeric( Parameter[NS,"lipids_pp"]),

	  #Metabolism
	  Cl       = as.numeric( Parameter[NS,"Cl"]),
	  Ke_bile  = as.numeric( Parameter[NS,"Ke_bile"]), 
	  Ke_urine = as.numeric( Parameter[NS,"Ke_urine"]),  

	  # Exposure parameter
	  WaterQuantity  =  1.21e+09, # microg/mL
	  IngestQuantity = 0,
	  ivQuantity     = 0,
	  
	  Vout = 1E12,
	  Texpo  = Texpo,# days
	  Times  = c(0,TimeM) # days
	  )

	tmpS = NULL
	tmpS = c(tmpS, Res$Sim[ 2:nrow(Res$Sim), "C_tot"] )                  # Total concentration
  tmpS = c(tmpS, Res$Sim[ 2:nrow(Res$Sim), "Concentration.C_art"])   # arterial conc
  tmpS = c(tmpS, Res$Sim[ 2:nrow(Res$Sim), "Concentration.C_liver"])     # liver conc
  return(c(NS, tmpS))
}

  # (2.3) the quit function
  #===============================================
  postpro <- function(x) { }

#=========================================================================================================================
# (5)                                 Calcul des sorties pour chaque vecteur de parametres
#========================================================================================================================= 
	#Test
   # NS = 14
   # Parameter=sa$X
   # TimeM   = TimeMesure
   # Texpo   = TimeExpo
   # parameters = as.numeric(Parameter[NS,])
   # Simulation(NS = 14, Parameter=sa$X, TimeM   = TimeMesure, Texpo   = TimeExpo )

  # PlanExp = read.table("PlanExp153.txt",sep="\t",dec=".",h=T)
	start_t = proc.time()

	# Creation d'un cluster
	cl <- makeCluster(processors)  
  
	#the initialization function
	invisible(parLapply(cl, 1:processors, prepro))

	for ( i in 1:length(Group.Sim) ){

			print(paste(i,"/",length(Group.Sim),sep=""))

			# Distribution des simulations
			result <- parSapply(cl = cl,  
						X = Group.Sim[[i]], 
						FUN = Simulation, 
						simplify = F,
						Parameter=sa$X, 
						TimeM   = TimeMesure,
						Texpo   = TimeExpo)

	# ecrit les resultats et re-initialise les objets
	SortieBrut=NULL
	for (z in 1:length(result)){ 
	if( length( result[[z]] ) < ntemp+1 ){ 
	print( result[[z]] )
	tmp = c(  result[[z]], rep(NA, ntemp + 1 - length( result[[z]] ) ))   
	}else{ tmp = result[[z]]  }
	  
	SortieBrut= rbind( SortieBrut, tmp)}
	
	colnames(SortieBrut) = c( "NumSim", 
	                          paste("Ci", TimeMesure, sep="."),
	                          paste("Cart", TimeMesure, sep="."),
	                          paste("Cliver", TimeMesure, sep=".") ) 
	
	rownames(SortieBrut) = 1:nrow(SortieBrut)

			save(SortieBrut, file = paste0(cheminWork,"SortieF",i,".RData")  )
			SortieBrut=NULL
			rm(result)

			}

	end_t = proc.time() - start_t
	end_t[1:5]
	print(end_t[1:5])

    # stop cluster
    stopCluster(cl)

save.image( paste0(cheminWork,"ASChronicPBPK.RData") )

# load( paste0(cheminWork,"ASChronicPBPK.RData") )
#=========================================================================================================================
# (5)                                 Calcul indices de sensibilite
#========================================================================================================================= 


Sortie.brut = NULL
for ( i in 1:length(Group.Sim) ){ # Group.Sim = 1:15
  load( paste0(cheminWork,"SortieF",i,".RData")  )
  summary(SortieBrut)
  Sortie.brut = rbind(Sortie.brut, SortieBrut) 
  rm(SortieBrut)
  }

NumSim = Sortie.brut[,1]
SimOrder = order(NumSim)  
Sortie.brut = Sortie.brut[SimOrder, ]
Sortie = Sortie.brut[,-1]

postscript("AS_Jansen_PBPK_Chronic.ps", paper = "a4", horizontal=F,pointsize=12)

#jpeg(file = "./res/SA_Chronic_Plot%d.jpeg",  width = 1344, height = 840, pointsize = 14, quality=60)

  par(mfrow=c(2,2), las=2, cex=0.7)

  FOI          = TI          = TI.borninf           = TI.bornsup          = matrix(NA, nrow = NP, ncol = ntemp*NOut)      
  rownames(FOI)= rownames(TI)= rownames(TI.borninf) = rownames(TI.bornsup)= Names

  
 for (i in 1: ncol(Sortie)) {
   print(i)
    tell(x=sa,y=Sortie[,i], nboot =  nboot, conf = 0.95)               
    FOI[,i]       = sa$S[,1]                                          # First order indices
    TI[,i]        = sa$T[,1]                                          # Total indices
    TI.borninf[,i] = sa$T[,4]
    TI.bornsup[,i] = sa$T[,5]                                         # Conf interval FOI
      
   plot(sa, main=colnames(Sortie)[i])                                 # 3 graphes (1 par sortie : un pour la croissance, un pour la maturit?
                                                                      # et un pour la reproduction) --> 3 r?sultats diff?rents 
 }

#=========================================================================================================================
# (6)                                         Graphiques
#========================================================================================================================= 

# Index moyenne de tous les temps
par(mfrow=c(1,1), las=2, cex=1, mai=c(1,0.7,0.35,0.1))

    for (i in 1:ncol(FOI)){
    FOI[ FOI[,i] < 0 , i ] =0 
    FOI[ FOI[,i] > TI[,i] , i ] =  TI[FOI[,i] > TI[,i] ,i] 
    }
    
FOI.t = apply(FOI, 1, mean.f)                                         # str(FOI.t) => conversion en matrice
TI.t = apply(TI, 1, mean.f)
sorting = order(TI.t, decreasing = F)
TI.t = TI.t[sorting]
FOI.t = FOI.t[sorting]

temp = t(cbind(FOI.t, TI.t))
barplot(temp, col=c("lightgrey","darkgrey"),horiz= F, beside=T , main="mean three times") 


# tous les temps 

  # Ci total
  FOI.L = as.matrix(FOI[,1:length(TimeMesure)])                       # as.matrix sinon pb de dimension avec apply
  TI.L  = as.matrix(TI[,1:length(TimeMesure)])
  
  FOI.L.t = apply(FOI.L, 1, mean.f)
  TI.L.t  = apply(TI.L,  1, mean.f)
  
  sorting = order(TI.L.t, decreasing = F)
  TI.L.t  = TI.L.t[sorting]
  FOI.L.t = FOI.L.t[sorting]

  FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
  tempC    = t(cbind(FOI.L.t, TI.L.t))
  #barplot(tempC, col=t(c("lightgrey","darkgrey")),horiz = T, beside =T , main="TI and FOI, Ci total") 

  # Ci arterial
  FOI.M = as.matrix(FOI[,2:2]) #as.matrix(FOI[,2])
  TI.M  = as.matrix(TI[,2:2])
  
  FOI.M.t = apply(FOI.M, 1, mean.f)
  TI.M.t  = apply(TI.M,  1, mean.f)
  
  sorting = order(TI.M.t, decreasing = F)
  TI.M.t  = TI.M.t[sorting]
  FOI.M.t = FOI.M.t[sorting]

  FOI.M.t = ifelse(FOI.M.t <= 0, 0, FOI.M.t)
  #TI.M.t = ifelse(FOI.t <= 0, 0, FOI.M.t)
  tempM = t(cbind(FOI.M.t, TI.M.t))
  #barplot(tempM, col=t(c("lightgrey","darkgrey")),horiz = T, beside =T , main="TI and FOI, C arterial") 

  # Ci liver
  FOI.R = as.matrix(FOI[,3 : ncol(FOI)])
  TI.R  = as.matrix(TI[,3 : ncol(TI)])
  
  FOI.R.t = apply(FOI.R, 1, mean.f)
  TI.R.t  = apply(TI.R,  1, mean.f)
  
  sorting =order( TI.R.t, decreasing = F)
  TI.R.t  = TI.R.t[sorting]
  FOI.R.t = FOI.R.t[sorting]

  FOI.R.t = ifelse(FOI.R.t <= 0, 0, FOI.R.t)
  tempR = t(cbind(FOI.R.t, TI.R.t))
  #barplot(tempR, col=t(c("lightgrey","darkgrey")),horiz = T, beside =T , main="TI and FOI, C liver") 
  
############################### Figure article ##################################
par(mfrow=c(1,3), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))

barplot(tempC, col=t(c("lightgrey","darkgrey")),horiz = T, beside =T , main="", 
        ylab="Ci total", cex.lab=1.5 , xlim=c(0,0.6) )  
mtext("[A]", line=-15,adj=0.5, cex=2)

barplot(tempM, col=t(c("lightgrey","darkgrey")),horiz = T, beside =T , main="", 
        ylab="C arterial", cex.lab=1.5, xlim=c(0,0.6)  ) 
mtext("[B]", line=-15,adj=0.5, cex=2)

barplot(tempR, col=t(c("lightgrey","darkgrey")),horiz = T, beside =T , main="", 
        ylab="C liver", cex.lab=1.5, xlim=c(0,0.6) ) 
mtext("[C]", line=-15,adj=0.5, cex=2)


############################### Lowry plots ####################################

par(mfrow=c(1,1))
	
  data <- data.frame(Parameter = Names, 
				Interaction = TI.t- FOI.t, 
				Main.Effect = FOI.t)
 
  #Fortify data and dump contents into plot function environment
  data_list <- fortify_lowry_data(data=data, 
  			          param_var = "Parameter",
			          main_var  = "Main.Effect",
			          inter_var = "Interaction")

  data_list$mdata = data_list$mdata[order(data_list$mdata$variable,decreasing=T),]
  
  sel = (length(Names)-9):length(Names)
  #sel = 1:length(Names)
  data_list$mdata = data_list$mdata[ data_list$mdata[,1] %in% Names[sel], ]
  data_list$data = data_list$data[ data_list$data$Parameter %in% Names[sel], ]
  
  p <- ggplot(data_list$data) +
    geom_bar(
    aes_string(x = data_list$mdata$Parameter, y = "value", fill = "variable"),
    data = data_list$mdata, 
    stat = "identity") +
    #scale_fill_manual(values = c("darkgrey","lightgrey")) +
    geom_ribbon( 
     aes(x = .numeric.param, ymin = .cumulative.main.effect, ymax = .valid.ymax),
      data = data_list$data,
      alpha = 0.5) +
    xlab("Parameters") +
    ylab("Total Effects (= Main Effects + Interactions)") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    scale_fill_grey(start = 0.5,end=0.2) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.direction = "horizontal") +
    guides(fill = guide_legend(reverse=TRUE))
  plot(p)
################################################################################

par(mfrow=c(2,3), las=1, mai=c(0.35,1,0.35,0.1))

for (i in 1:(ntemp*NOut)){
  sorting  = order(TI[,i], decreasing = F)
  TI.temp  = TI[sorting,i]
  FOI.temp = FOI[sorting,i]
  TI.bornsup.temp = TI.bornsup[sorting,i]
  TI.borninf.temp = TI.borninf[sorting,i]
  
  BorneInf = t(cbind( rep(0,nrow(TI.borninf)), TI.borninf.temp) )
  BorneSup = t(cbind( rep(0,nrow(TI.borninf)), TI.bornsup.temp) )
  
  temp = t(cbind(  FOI.temp,TI.temp))
  barplot2(temp,  col=t(c("lightgrey","darkgrey")),
  horiz = T,beside =T ,cex.names=0.7,	 
  main=paste("TI and FOI", colnames(Sortie)[i] ),
           plot.ci = TRUE,   ci.l = BorneInf,  ci.u = BorneSup )  
}


dev.off()
