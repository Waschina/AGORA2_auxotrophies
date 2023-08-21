library(sybilSBML)
library(cplexAPI)
library(MicrobiomeGS2)
library(parallel)

SYBIL_SETTINGS("SOLVER", "cplexAPI")
cl <- makeCluster(10)


AAs <- c("ala","val","leu","ile","pro","cys","met","gly","trp","phe","lys","arg",
         "his","tyr","thr","glu","gln","asp","asn","ser")
aa_exIDs <- lapply(AAs, function(x) {
  # single free AA
  if(x != "gly")
    out <- c(paste0("EX_",x,"_L__40__e__41__"),
             paste0("EX_",x,"_D__40__e__41__"))
  if(x == "gly")
    out <- "EX_gly__40__e__41__"
  
  # dipeptides
  dipep1 <- paste0("EX_",x,AAs,"__40__e__41__")
  dipep2 <- paste0("EX_",AAs,x,"__40__e__41__")
  
  return(c(out, dipep1, dipep2))
})
names(aa_exIDs) <- AAs

model_files <- dir("/mnt/nuuk/2023/AGORA2/RDS/", full.names = TRUE)
models <- lapply(model_files, readRDS)
names(models) <- gsub("\\.RDS$","",basename(model_files))

clusterExport(cl, c("AAs","aa_exIDs"))

auxout <- parLapply(cl, models, fun = function(mod) {
  library(sybil)
  SYBIL_SETTINGS("SOLVER", "cplexAPI")
  auxpred <- rep(1, length(AAs))
  names(auxpred) <- AAs
  grWT <- optimizeProb(mod)
  for(aai in AAs) {
    rel_ex <- intersect(mod@react_id, aa_exIDs[[aai]])
    lbs <- mod@lowbnd[match(rel_ex,mod@react_id)]
    
    # keep other part of dipeptides
    keep_ex <- gsub(aai,"",rel_ex)
    keep_ex <- gsub("^EX_|__40__e__41__","",keep_ex)
    
    mod_mod <- mod
    
    for(keepi in keep_ex) {
      if(keepi != "gly")
        keepi <- paste0(keepi, "_L")
      
      mod_mod <- addReact(mod_mod,
                          id = paste0("IF_adhoc_", keepi),
                          met = paste0(keepi,"__91__c__93__"),
                          Scoef = 1)
      
    }
    
    mod_mod <- rmReact(mod_mod, rel_ex)
    
    # actual auxotrophy prediction via growth estimate
    grRed <- optimizeProb(mod_mod)
    if(grRed@lp_obj < grWT@lp_obj*1e-12)
      auxpred[aai] <- 0
  }
  return(auxpred)
})

stopCluster(cl)


auxout <- do.call("rbind", auxout)

zut <- as.data.table(auxout)
zut$AGORA2 <- rownames(auxout)

fwrite(zut, "AGORA2_predicted_auxotrophies.csv")
