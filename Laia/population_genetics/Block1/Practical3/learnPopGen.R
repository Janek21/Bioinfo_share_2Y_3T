selection <- function (p0 = 0.01, w = c(1, 0.9, 0.8), time = 100, show = "p", 
    pause = 0, ...) 
{
    if (hasArg(add)) 
        add <- list(...)$add
    else add <- FALSE
    if (hasArg(color)) 
        color <- list(...)$color
    else color <- "black"
    if (hasArg(equil)) 
        equil <- list(...)$equil
    else equil <- FALSE
    if (hasArg(lwd)) 
        lwd <- list(...)$lwd
    else lwd <- 2
    if (hasArg(lty)) 
        lty <- list(...)$lty
    else lty <- "solid"
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(xlim)) 
        xlim <- list(...)$xlim
    else xlim <- NULL
    if (w[2] > w[1] && w[2] > w[3]) {
        W <- w/w[2]
        eq <- (1 - W[3])/((1 - W[1]) + (1 - W[3]))
    }
    else if (w[2] < w[1] && w[2] < w[3]) {
        W <- w/w[2]
        eq <- (1 - W[3])/((1 - W[1]) + (1 - W[3]))
    }
    else if (1 %in% which(w == max(w))) {
        eq <- 1
    }
    else eq <- 0
    if (show == "surface") {
        p <- 0:100/100
        wbar <- p^2 * w[1] + 2 * p * (1 - p) * w[2] + (1 - p)^2 * 
            w[3]
        plot(p, wbar, type = "l", xlim = xlim, ylim = c(0, 
            1.1 * max(w)), main = if (is.null(main)) 
            expression(paste("mean fitness (", bar(w), 
                ")", sep = ""))
        else main, ylab = expression(bar(w)), col = color)
        if (equil) 
            abline(v = eq, lty = "dotted")
    }
    else if (show == "deltap") {
        p <- 0:100/100
        wbar <- p^2 * w[1] + 2 * p * (1 - p) * w[2] + (1 - p)^2 * 
            w[3]
        deltap <- (p/wbar) * (p * w[1] + (1 - p) * w[2] - wbar)
        plot(p, deltap, type = "l", xlim = xlim, main = if (is.null(main)) 
            expression(paste(Delta, "p as a function of p", 
                sep = ""))
        else main, ylab = expression(paste(Delta, "p", 
            sep = "")), col = color)
        lines(c(0, 1), c(0, 0), lty = 2)
        if (equil) 
            abline(v = eq, lty = "dotted")
    }
    else {
        if (show == "cobweb") {
            p <- 0:100/100
            wbar <- p^2 * w[1] + 2 * p * (1 - p) * w[2] + (1 - 
                p)^2 * w[3]
            p2 <- (p/wbar) * (p * w[1] + (1 - p) * w[2] - wbar) + 
                p
            plot(p, p2, type = "l", xlim = xlim, xlab = expression(p[t]), 
                ylab = expression(p[t + 1]), main = if (is.null(main)) 
                  expression(paste(p[t + 1], " as a function of ", 
                    p[t], sep = ""))
                else main, col = color)
            lines(c(0, 1), c(0, 1), lty = 2)
            if (equil) {
                abline(v = eq, lty = "dotted")
                abline(h = eq, lty = "dotted")
            }
            dev.flush()
        }
        p <- wbar <- vector()
        p[1] <- p0
        wbar[1] <- p[1]^2 * w[1] + 2 * p[1] * (1 - p[1]) * w[2] + 
            (1 - p[1])^2 * w[3]
        for (i in 2:time) {
            p[i] <- p[i - 1]
            p[i] <- (p[i]^2 * w[1] + p[i] * (1 - p[i]) * w[2])/wbar[i - 
                1]
            wbar[i] <- p[i]^2 * w[1] + 2 * p[i] * (1 - p[i]) * 
                w[2] + (1 - p[i])^2 * w[3]
            ii <- (i - 1):i
            if (show == "p") {
                if (i == 2 && !add) 
                  plot(1:i, p, type = "l", xlim = if (is.null(xlim)) 
                    c(0, time)
                  else xlim, ylim = c(0, 1), xlab = "time", 
                    main = if (is.null(main)) 
                      "frequency of A"
                    else main, col = color, lwd = lwd, lty = lty)
                else lines(ii, p[ii], type = "l", col = color, 
                  lwd = lwd, lty = lty)
            }
            else if (show == "q" && !add) {
                if (i == 2) 
                  plot(1:i, 1 - p, type = "l", xlim = if (is.null(xlim)) 
                    c(0, time)
                  else xlim, ylim = c(0, 1), xlab = "time", 
                    ylab = "q", main = if (is.null(main)) 
                      "frequency of a"
                    else main, col = color, lwd = lwd, lty = lty)
                else lines(ii, 1 - p[ii], type = "l", col = color, 
                  lwd = lwd, lty = lty)
            }
            else if (show == "fitness" && !add) {
                if (i == 2) 
                  plot(1:i, wbar, type = "l", xlim = if (is.null(xlim)) 
                    c(0, time)
                  else xlim, ylim = c(0, 1.1 * max(w)), xlab = "time", 
                    main = if (is.null(main)) 
                      expression(paste("mean fitness (", 
                        bar(w), ")", sep = ""))
                    else main, ylab = expression(bar(w)), col = color)
                else lines(ii, wbar[ii], type = "l", col = color)
            }
            else if (show == "cobweb") {
                lines(c(p[i - 1], p[i - 1]), c(p[i - 1], p[i]), 
                  col = color)
                lines(c(p[i - 1], p[i]), c(p[i], p[i]), col = color)
            }
            else {
                message("not a recognized option")
                break
            }
            if (equil) {
                if (show == "p") 
                  abline(h = eq, lty = "dotted")
                if (show == "q") 
                  abline(h = 1 - eq, lty = "dotted")
                if (show == "fitness") 
                  abline(h = eq^2 * w[1] + 2 * eq * (1 - eq) * 
                    w[2] + (1 - eq)^2 * w[3], lty = "dotted")
            }
            dev.flush()
            Sys.sleep(pause)
        }
    }
    object <- list(p0 = p0, w = w, time = time, p = if (show %in% 
        c("surface", "deltap")) NULL else p, equilibrium = eq)
    class(object) <- "selection"
    invisible(object)
}

drift.selection <- function (p0 = 0.5, Ne = 100, w = c(1, 1, 1), ngen = 400, nrep = 10, 
    colors = NULL, ...) 
{
    if (is.null(colors)) 
        colors <- rainbow(nrep)
    w <- (w/max(w))[3:1]
    gametes <- rep(0, 2 * Ne)
    if (p0 > 0) 
        gametes[1:round(p0 * 2 * Ne)] <- 1
    gametes <- replicate(nrep, gametes, simplify = FALSE)
    p <- lapply(gametes, mean)
    for (i in 1:ngen) {
        genotypes <- lapply(gametes, function(x) matrix(sample(x), 
            length(x)/2, 2))
        fitness <- lapply(genotypes, function(x, w) w[rowSums(x) + 
            1], w = w)
        selected <- lapply(fitness, function(prob, x) cbind(sample(x, 
            prob = prob, replace = TRUE), sample(x, prob = prob, 
            replace = TRUE)), x = Ne)
        copy <- replicate(nrep, matrix(sample(1:2, 2 * Ne, replace = TRUE), 
            Ne, 2), simplify = FALSE)
        gametes <- mapply(function(g, s, c) c(diag(g[s[, 1], 
            ][, c[, 1]]), diag(g[s[, 2], ][, c[, 2]])), g = genotypes, 
            s = selected, c = copy, SIMPLIFY = FALSE)
        for (j in 1:nrep) p[[j]][i + 1] <- mean(gametes[[j]])
    }
    plot(0:ngen, p[[1]], type = "l", col = colors[1], lwd = 2, 
        ylim = c(0, 1), xlab = "time (generations)", ylab = "f(A)")
    if (nrep > 1) 
        nulo <- mapply(lines, x = replicate(nrep - 1, 0:ngen, 
            simplify = FALSE), y = p[2:nrep], col = colors[2:nrep], 
            lwd = 2)
    class(p) <- "drift.selection"
    attr(Ne, "Ne") <- Ne
    attr(p, "w") <- w[3:1]
    invisible(p)
}