.internal_xmapdf_print = \(x,keyname = "libsizes",significant = FALSE){
  resdf = x[[1]]
  bidirectional = x$bidirectional
  if (bidirectional){
    if (significant) {
      resdf = resdf[resdf$x_xmap_y_sig < 0.05 & resdf$y_xmap_x_sig < 0.05,
                    c(keyname, "y_xmap_x_mean", "x_xmap_y_mean"),drop = FALSE]
    } else {
      resdf = resdf[,c(keyname, "y_xmap_x_mean", "x_xmap_y_mean"),drop = FALSE]
    }
    names(resdf) = c(keyname,
                     paste0(x$varname[1], "->", x$varname[2]),
                     paste0(x$varname[2], "->", x$varname[1]))
  } else {
    if (significant) {
      resdf = resdf[resdf$y_xmap_x_sig < 0.05,c(keyname, "y_xmap_x_mean"),drop = FALSE]
    } else {
      resdf = resdf[,c(keyname, "y_xmap_x_mean"),drop = FALSE]
    }
    names(resdf) = c(keyname,
                     paste0(x$varname[1], "->", x$varname[2]))
  }
  return(resdf)
}

#' print ccm result
#' @noRd
#' @export
print.ccm_res = \(x,significant = FALSE,...){
  print(.internal_xmapdf_print(x,significant = significant))
}

#' print cmc result
#' @noRd
#' @export
print.cmc_res = \(x,significant = FALSE,...){
  print(.internal_xmapdf_print(x,"neighbors",significant = significant))
}

#' print pcm result
#' @noRd
#' @export
print.pcm_res = \(x,significant = FALSE,...){
  pxmap = x[-2]
  xmap = x[-1]

  cat('-------------------------------------- \n')
  cat("***partial cross mapping prediction*** \n")
  cat('-------------------------------------- \n')
  print(.internal_xmapdf_print(pxmap,significant = significant))
  cat("\n------------------------------ \n")
  cat("***cross mapping prediction*** \n")
  cat('------------------------------ \n')
  print(.internal_xmapdf_print(xmap,significant = significant))
}

#' print xmap_self result
#' @noRd
#' @export
print.xmap_self = \(x,...){
  res = as.matrix(x$xmap)
  if (x$method == "smap"){
    cat(paste0("The suggested theta for variable ", x$varname, " is ", OptThetaParm(res)), "\n")
  } else {
    if (x$method == "simplex"){
      res = OptEmbedDim(res)
    } else {
      res = OptICparm(res)
    }
    cat(paste0("The suggested E and k for variable ", x$varname, " is ", res[1], " and ", res[2]), "\n")
    if (res[1] == 1 && x$tau == 0) warning("When tau = 0, E should not be 1")
  }
}

#' plot ccm result
#' @noRd
#' @export
plot.ccm_res = \(x, family = "serif", legend_texts = NULL,
                 legend_cols = c("#ed795b","#608dbe"),
                 draw_ci = FALSE, ci_alpha = 0.25,
                 xbreaks = NULL, xlimits = NULL,
                 ybreaks = seq(0, 1, by = 0.1),
                 ylimits = c(-0.05, 1),
                 ylabel = expression(rho), ...){
  resdf = x[[1]]
  bidirectional = x$bidirectional

  if(is.null(xbreaks)) xbreaks = resdf$libsizes
  if(is.null(xlimits)) xlimits = c(min(xbreaks)-1,max(xbreaks)+1)
  if (is.null(legend_texts)){
    pval = resdf |>
      dplyr::slice_tail(n = 1) |>
      dplyr::select(x_xmap_y_sig,y_xmap_x_sig) |>
      unlist() |>
      round(3)
    legend_texts = c(paste0(x$varname[2], " xmap ", x$varname[1], ", p = ", pval[2]),
                     paste0(x$varname[1], " xmap ", x$varname[2], ", p = ", pval[1]))
  }
  legend_texts = .check_inputelementnum(legend_texts,2)
  legend_cols = .check_inputelementnum(legend_cols,2)
  names(legend_cols) = c("x - y","y - x")

  ci_alpha = .check_inputelementnum(ci_alpha,2)

  fig1 = ggplot2::ggplot(data = resdf,
                         ggplot2::aes(x = libsizes)) +
    ggplot2::geom_line(ggplot2::aes(y = y_xmap_x_mean,
                                    color = "x - y"),
                       lwd = 1.25)

  if (draw_ci) {
    fig1 = fig1 +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = y_xmap_x_lower,
                                        ymax = y_xmap_x_upper),
                           alpha = ci_alpha[1], fill = legend_cols[1])
  }

  if (bidirectional){
    fig1 = fig1 + ggplot2::geom_line(ggplot2::aes(y = x_xmap_y_mean,
                                                  color = "y - x"),
                                     lwd = 1.25)

    if (draw_ci) {
      fig1 = fig1 +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = x_xmap_y_lower,
                                          ymax = x_xmap_y_upper),
                             alpha = ci_alpha[2], fill = legend_cols[2])
    }
  }

  fig1 = fig1 +
    ggplot2::scale_x_continuous(breaks = xbreaks, limits = xlimits,
                                expand = c(0, 0), name = "Library size") +
    ggplot2::scale_y_continuous(breaks = ybreaks, limits = ylimits,
                                expand = c(0, 0), name = ylabel) +
    ggplot2::scale_color_manual(values = legend_cols,
                                labels = legend_texts,
                                name = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(family = family),
                   axis.text.x = ggplot2::element_text(angle = 30),
                   axis.title = ggplot2::element_text(family = family),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "inside",
                   legend.justification = c(0.05,1),
                   legend.background = ggplot2::element_rect(fill = 'transparent'),
                   legend.text = ggplot2::element_text(family = family))
  return(fig1)
}

#' plot cmc result
#' @noRd
#' @export
plot.cmc_res = \(x, ...){
  xmap = x[-1]
  class(xmap) = "ccm"
  draw_ci = FALSE
  fig1 = plot.ccm_res(xmap,draw_ci = draw_ci,ylabel = "Causal Score",...)
  return(fig1)
}

#' plot pcm result
#' @noRd
#' @export
plot.pcm_res = \(x, partial = TRUE, ...){
  indice = ifelse(partial,-2,-1)
  xmap = x[indice]
  class(xmap) = "ccm"
  fig1 = plot.ccm_res(xmap,...)
  return(fig1)
}
