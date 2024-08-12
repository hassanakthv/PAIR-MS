#
# This file is knitted from ../notebooks/lab_notebook.Rmd

run_model <- function(df, val_func, excl_func=NULL) {
  peaks <- unique(df$peak)
  files <- unique(df$file)

  coll_df <- data.frame("file"=character(),
                        "peak"=character(),
                        "value"=numeric(),
                        "group"=character())

  for (u_peak in peaks) {
    for (u_file in files) {
      input_df <- subset(df, peak==u_peak & file==u_file)
      if (!is.null(excl_func)) {
        input_df <- excl_func(input_df)
      }
      group <- unique(input_df$group)[1]
      output_val <- val_func(input_df$gamma)
      coll_entry <- data.frame("file"=u_file, "peak"=u_peak, "value"=output_val, "group"=group)
      coll_df <- rbind(coll_df, coll_entry)
    }
  }
  names(coll_df) <- c("file", "peak", "value", "group")
  return(coll_df)
}


get_cov <- function(coll_df) {

  peaks <- unique(coll_df$peak)
  groups <- unique(coll_df$group)

  cov_df <- data.frame("group"=character(),
                      "peak"=character(),
                      "cov"=numeric())

  for (u_peak in peaks) {
    for (u_group in groups) {
      group_data <- subset(coll_df, peak==u_peak & group==u_group)
      group_mean <- mean(group_data$value)
      group_sd <- sd(group_data$value)
      cov <- group_sd / group_mean * 100
      cov_entry <- data.frame("group"=u_group, "peak"=u_peak, "cov"=cov)
      cov_df <- rbind(cov_df, cov_entry)
    }
  }

  return(cov_df)

}

plot_dens_values <- function(df, coll_df) {
  n_groups = length(unique(df$group))

  ggplot(df, aes(x=gamma, col=peak)) +
    scale_y_sqrt() + scale_x_sqrt() + # make visualisation more obvious
    geom_density() +
    facet_wrap(c("file", "group"), as.table=TRUE, ncol=n_groups) +
    theme_minimal() +
    geom_vline(data=coll_df, aes(xintercept=value, col=peak))
}

run_analysis <- function(df, val_func_in, excl_func_in, val_func_name=NULL, title=NULL) {
  coll_df <- run_model(df, val_func=val_func_in, excl_func=excl_func_in)
  p <- plot_dens_values(df, coll_df)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  print(p)

  cov_entry <- get_cov(coll_df)
  ifelse (is.null(val_func_name),
          cov_entry$func <- deparse(substitute(val_func_in)),
          cov_entry$func <- val_func_name)
  cov_entry$excl <- deparse(substitute(excl_func_in))

  print(kable(coll_df))
  print(kable(cov_entry))

  return(cov_entry)
}


fit_gauss <- function(df) {
  res <- fitdist(df, dnorm, method="mle", lower=0)$estimate['mean']
  return(res)
}

fit_weib <- function(df) {
  fit <- fitdist(df, "weibull")
  # return the weibull peak (i.e. mode)
  shape <- fit$estimate[["shape"]]
  scale <- fit$estimate[["scale"]]
  res <- scale*(1-1/shape)^(1/shape)
  return(res)
}

fit_gamma <- function(df) {
  fit <- fitdistr(df, "gamma")
  # return gamma peak (i.e. mode)
  shape <- fit$estimate[["shape"]]
  rate <- fit$estimate[["rate"]]
  res <- (shape - 1) / rate
  return(res)
}

max_dens <- function(df) {
  df_dens <- density(df)
  xs <- df_dens[["x"]]
  ys <- df_dens[["y"]]
  maxval <- xs[ys==max(ys)]
  return(maxval)
}

strong_grad <- function (df) {
  dens <- density(df)
  xs <- dens[["x"]]
  ys <- dens[["y"]]

  grad <- diff(ys)
  grad_max_y <- ys[(grad==max(grad))]
  grad_min_y <- ys[(grad==min(grad))]

  # take the more exclusive gradient
  thresh <- max(grad_min_y, grad_max_y)

  ys <- ys - thresh
  ys[ys<0] <- 0
  ys <- ys/sum(ys)
  res <- sum(xs*ys)
  return(res)
}

hz_func <- function(df, prop) {
  if (prop == 0) {return(mean(df))}
  else {
    dens <- density(df)
    xs <- dens[["x"]]
    ys <- dens[["y"]]

    thresh <- (prop/100) * max(ys)

    ys <- ys - thresh
    ys[ys<0] <- 0
    ys <- ys/sum(ys)
    res <- sum(xs*ys)
    return(res)
  }
}

hz_mean_50 <- function(df) {
  return(hz_func(df, 50))
}

hz_mean_25 <- function(df) {
  return(hz_func(df, 25))
}

hz_mean_75 <- function(df) {
  return(hz_func(df, 75))
}

hz_one_peak <- function(df) {
  dens <- density(df)
  xs <- dens[["x"]]
  ys <- dens[["y"]]

  local_maxs_y <- ys[diff(sign(diff(ys)))==-2]
  if (length(local_maxs_y)==1) {
    ys <- ys/sum(ys)
    value <- sum(xs*ys)
  } else {
    sorted_maxs_y <- sort(local_maxs_y, decreasing=TRUE)
    scnd_peak <- sorted_maxs_y[[2]]
    ys_corr <- ys - scnd_peak
    ys_corr[ys_corr<0] <- 0
    ys_corr <- ys_corr/sum(ys_corr)
    value <- sum(xs*ys_corr)
  }
  return(value)
}

excl_iqm <- function(df) {
  qs <- quantile(df$gamma)
  out_df <- subset(df, gamma>qs['25%'] & gamma<qs['75%'])
  return (out_df)
}


