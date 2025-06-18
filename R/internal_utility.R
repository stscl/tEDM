.check_character = \(x,...){
  xstrs = c(x,...)
  for (i in xstrs){
    if (!inherits(i,"character")) {
      stop("Please check the characters in the function input.")
    }
  }
  return(xstrs)
}

.check_inputelementnum = \(x,n,condsnum = NULL){
  if (is.null(condsnum) || length(x) == 1){
    res = rep(x,length.out = n)
  } else if (length(x) == 2) {
    res = c(rep(x[1],2),rep(x[2],condsnum))
  } else {
    res = c(x[1:2],rep(x[c(-1,-2)],length.out = condsnum))
  }
  return(res)
}

.check_parallellevel = \(parallel.level){
  pl = 0
  if (parallel.level != "low"){
    pl = 1
  }
  return(pl)
}

.internal_library = \(data){
  nnaindice = which(apply(is.na(data),1,\(.x) !any(.x)))
  return(nnaindice)
}

.uni_ts = \(data,target){
  if (is.null(target)) return(rep(0,nrow(data)))
  target = .check_character(target)
  res = data[,"target",drop = TRUE]
  return(res)
}

.multivar_ts = \(data,columns){
  columns = .check_character(columns)
  res = as.matrix(data[,columns,drop = FALSE])
  return(res)
}
