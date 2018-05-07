sumup2binom <- function(pos, ntot) {
  if(pos >= ntot) {
    warning("Number of positives is larger or equal to total number of individuals!")
    return(-1)
  }

  res <- c(rep(1, pos), rep(0, ntot-pos))
  res
}

sumup2binom_ages <- function(pos, tot, ages) {
  if(!(length(pos) ==  length(ages)) && (length(ages) == length(tot))) {
    warning("Length of provides vectors are not equal to each other!")
    return(-1)
  }

  total <- sum(tot)
  index <- 1

  bin_res <- rep(0, total)
  age_res <- rep(0, total)

  for(i in 1:length(pos)) {
    bin_res[index:(index+tot[i]-1)] <- c(rep(1, pos[i]), rep(0, tot[i]-pos[i]))
    age_res[index:(index+tot[i]-1)] <- rep(ages[i], tot[i])
    index = index + tot[i]
  }

  return(list(bin=bin_res, age=age_res))
}

binom2sumup <- function(data) {
  if(length(data) == 0) {
    warning("Length of data must be larger than zero!")
    return(-1)
  }

  res <- list(pos=sum(data), ntot=length(data))
  res
}

binom2sumup_ages <- function(bin, ages) {
  if(!(length(bin) == length(ages))) {
    warning("Length of provides vectors are not equal to each other!")
    return(-1)
  }

  unique_ages <- unique(ages)
  num_unique_ages <- length(unique_ages)

  pos_res <- rep(0, num_unique_ages)
  tot_res <- rep(0, num_unique_ages)
  age_res <- rep(0, num_unique_ages)

  for(i in 1:num_unique_ages) {
    ind_ages <- which(unique_ages[i] == ages)

    pos_res[i] <- sum(bin[ind_ages])
    tot_res[i] <- length(ind_ages)
    age_res[i] <- unique_ages[i]
  }

  return(list(pos=pos_res, tot=tot_res, age=age_res))
}
