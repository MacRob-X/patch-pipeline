# Checking if the 2016 and 2019 taxonomies are identical
# Robert MacDonald
# 04 April 2024

# read in taxonomies 
taxo2016 <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2016_02_29.csv", strings = F)
taxo2019 <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv", strings = F)

# check if each column is identical in both taxonomies
for (i in 1:ncol(taxo2016)){
  if(identical(taxo2016[, i], taxo2019[, i])){
    print(paste0("Column ", i, " identical in both taxonomies."))
  }
}


# all columns are identical