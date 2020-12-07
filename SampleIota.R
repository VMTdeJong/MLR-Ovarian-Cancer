# Gives a random split of the iota data set.
# iota is assumed to be the iota data set
# n.small is the desired size of the smallest outcome category.
# Returns a list of 2.
#   "inc" is the data included for development
#   "ex"  is the data excluded for development, possibly to be used for validation.

SampleIota <- function(iota, n.small) {
  # Sizes of categories of full data set
  nt.small  <- min(table(iota$outcome))
  nt.middle <- median(table(iota$outcome))
  nt.large  <- max(table(iota$outcome))
  
  # Sampling fraction
  fraction  <- n.small / nt.small
  n.middle  <- round(fraction * nt.middle)
  n.large   <- round(fraction * nt.large)
  
  # Names of outcomes
  which.small  <- which.min(table(iota$outcome))
  which.large  <- which.max(table(iota$outcome))
  which.middle <- table(iota$outcome)[-c(which.small, which.large)]
  
  name.small  <- names(which.small)
  name.middle <- names(which.middle)
  name.large  <- names(which.large)
  
  # Indexes of observations per category
  full.i.small  <- which(iota$outcome == name.small)
  full.i.middle <- which(iota$outcome == name.middle)
  full.i.large  <- which(iota$outcome == name.large)
  
  # Indexes of sampled observations
  i.small   <- sample(full.i.small, n.small)
  i.middle  <- sample(full.i.middle, n.middle)
  i.large   <- sample(full.i.large, n.large)
  
  list("inc" = iota[ c(i.small, i.middle, i.large), ],
       "ex"  = iota[-c(i.small, i.middle, i.large), ])
}