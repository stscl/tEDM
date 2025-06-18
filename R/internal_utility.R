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

.internal_varname = \(conds = NULL){
  .varname = c("cause","effect")
  return(.varname)
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

.internal_xmapdf_binding = \(varname,x_xmap_y,y_xmap_x,bidirectional,keyname = "libsizes"){
  colnames(y_xmap_x) = c(keyname,"y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = as.data.frame(y_xmap_x)

  if (bidirectional){
    colnames(x_xmap_y) = c(keyname,"x_xmap_y_mean","x_xmap_y_sig",
                           "x_xmap_y_upper","x_xmap_y_lower")
    x_xmap_y = as.data.frame(x_xmap_y)
    resdf = x_xmap_y |>
      dplyr::full_join(y_xmap_x, by = keyname) |>
      dplyr::arrange({{keyname}})
  } else {
    resdf = dplyr::arrange(y_xmap_x,{{keyname}})
  }

  return(resdf)
}

.bind_xmapdf = \(varname,x_xmap_y,y_xmap_x,bidirectional){
  resdf = .internal_xmapdf_binding(varname,x_xmap_y,y_xmap_x,bidirectional)
  res = list("xmap" = resdf, "varname" = varname, "bidirectional" = bidirectional)
  class(res) = 'ccm_res'
  return(res)
}

.bind_xmapdf2 = \(varname,x_xmap_y,y_xmap_x,bidirectional){

  tyxmapx = y_xmap_x[,c(1,2,4:6),drop = FALSE]
  dyxmapx = y_xmap_x[,c(1,3,7:9),drop = FALSE]
  txxmapy = NULL
  dxxmapy = NULL
  if(bidirectional){
    txxmapy = x_xmap_y[,c(1,2,4:6),drop = FALSE]
    dxxmapy = x_xmap_y[,c(1,3,7:9),drop = FALSE]
  }

  txmap = .internal_xmapdf_binding(varname[1:2],txxmapy,tyxmapx,bidirectional)
  dxmap = .internal_xmapdf_binding(varname[1:2],dxxmapy,dyxmapx,bidirectional)

  res = list("pxmap" = dxmap, "xmap" = txmap,
             "varname" = varname[1:2],
             "bidirectional" = bidirectional)
  class(res) = 'pcm_res'
  return(res)
}

.bind_intersectdf = \(varname,x_xmap_y,y_xmap_x,bidirectional){
  resdf = .internal_xmapdf_binding(varname,x_xmap_y,y_xmap_x,bidirectional,keyname = "neighbors")
  res = list("xmap" = resdf, "varname" = varname, "bidirectional" = bidirectional)
  class(res) = 'cmc_res'
  return(res)
}

.bind_xmapself = \(x,varname,...){
  res = list("xmap" = x,"varname" = varname)
  class(res) = "xmap_self"
  return(res)
}
