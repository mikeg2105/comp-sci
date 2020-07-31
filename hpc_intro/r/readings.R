main <- function() {
  
  dat <- read.csv(file = "../data/inflammation-01.csv", header = FALSE)
  mean_per_patient <- apply(dat, 1, mean)
  cat(mean_per_patient, sep = "\n")
}

main()