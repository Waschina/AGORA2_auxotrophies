library(sybilSBML)
errlog <- file.create("errlog.txt")
files <- dir("/mnt/nuuk/2023/AGORA2/sbml/", full.names = TRUE)

for(file in files)  {
  cat(basename(file),"\n")
  # (1) Read sbml
  mod <- tryCatch(readSBMLmod(file, validateSBML = FALSE), error = function(cond) { return(NA) })
  
  if(class(mod) == "modelorg") {
    # (2) export RDS
    exfile <- paste0("/mnt/nuuk/2023/AGORA2/RDS/", gsub("\\.xml\\.gz$",".RDS",basename(file)))
    saveRDS(mod, exfile)
  } else {
    errm <- validateSBMLdocument(file)
    errm <- paste0("\n\nFile: ",file,"\n",
                   "line: ", errm@sbmlErrors[[1]]$line,", column: ", errm@sbmlErrors[[1]]$column,"\n",
                   "message: ",errm@sbmlErrors[[1]]$message)
    cat(errm, file = "errlog.txt", append = TRUE)
  }
  gc()
}