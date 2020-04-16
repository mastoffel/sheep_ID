get_raneff <- function (inla_mod, scales = c("var", "sd")) {
        if (length(scales) == 2) scales <- scales[[1]]
        
        irp <- inla_mod$internal.marginals.hyperpar
        hrp <- inla_mod$marginals.hyperpar
        hypnames <- names(irp)
        iip <- grep("precision", hypnames)
        for (i in 1:length(irp)) {
                if (i %in% iip) {
                        irp[[i]] <- hyper_sd(irp[[i]], scales, internal = TRUE)
                }
                else {
                        irp[[i]] <- hrp[[i]]
                        hypnames[i] <- names(hrp)[i]
                }
        }
        #ts <- t(sapply(irp, density_summary))
        ts <- purrr::map_df(irp, density_summary)
        hypnames <- sub("Log precision", scales, hypnames)
        ts <- cbind(term = hypnames, ts) %>% 
                tidyr::separate(col = term, into = c("scale", "term"), sep = " for ") 
        tibble::as_tibble(ts)[c(2,1,3:8)]
}

hyper_sd <- function (prec, scales, internal = FALSE) {
        if (scales == "var"){
                trans_fun <- function(x) 1/x
        } else if (scales == "sd") {
                trans_fun <- function(x) 1/sqrt(x)
        } else {
                stop("scales has to be var or sd")
        }
        
        if (internal) {
                inla.tmarginal(function(x)  trans_fun(exp(x)), prec)
        }
        else {
                inla.tmarginal(function(x)  trans_fun(x), prec)
        }
}

density_summary <- function(dens) {
        m <- INLA::inla.emarginal(function(xx) c(xx, xx^2), dens)
        q <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), dens)
        s <- sqrt(max(0, m[2] - m[1]^2))
        md <- INLA::inla.mmarginal(dens)
        c(mean = m[1], std_err = s, ci_lower = q[1], median = q[2], ci_upper = q[3], 
          mode = md)
}
