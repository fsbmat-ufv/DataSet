# Versão corrigida da função HeckmanBS com proteção para rho e start
HeckmanBS_safe <- function(selection, outcome, data = sys.frame(sys.parent()), start = NULL) {
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("selection", "data", "subset"), names(mf), 0)
        mfS <- mf[c(1, m)]
        mfS$drop.unused.levels <- TRUE
        mfS$na.action <- na.pass
        mfS[[1]] <- as.name("model.frame")
        names(mfS)[2] <- "formula"
        mfS <- eval(mfS, parent.frame())
        mtS <- terms(mfS)
        XS <- model.matrix(mtS, mfS)
        NXS <- ncol(XS)
        YS <- model.response(mfS)
        YSLevels <- levels(as.factor(YS))
        
        m <- match(c("outcome", "data", "subset", "weights", "offset"), names(mf), 0)
        mfO <- mf[c(1, m)]
        mfO$na.action <- na.pass
        mfO$drop.unused.levels <- TRUE
        mfO[[1]] <- as.name("model.frame")
        names(mfO)[2] <- "formula"
        mfO <- eval(mfO, parent.frame())
        mtO <- attr(mfO, "terms")
        XO <- model.matrix(mtO, mfO)
        NXO <- ncol(XO)
        YO <- model.response(mfO)
        
        if (is.null(start)) {
                start <- rep(0.1, NXS + NXO + 2)
        }
        if (any(!is.finite(start))) stop("❌ ERRO: 'start' contém valores não finitos.")
        
        loglik_BS <- function(par) {
                print(par)
                NXS <- ncol(XS)
                NXO <- ncol(XO)
                gamma <- par[1:NXS]
                beta <- par[(NXS + 1):(NXS + NXO)]
                phi1 <- par[NXS + NXO + 1]
                rho  <- par[NXS + NXO + 2]
                
                if (!is.finite(rho) || rho <= -1 || rho >= 1) {
                        cat("⚠️ rho inválido:", rho, "\n")
                        return(NA)
                }
                if (phi1 < 0 || !is.finite(phi1)) {
                        cat("⚠️ phi1 inválido:", phi1, "\n")
                        return(NA)
                }
                
                if (any(!is.finite(c(XS, XO, YO)))) {
                        cat("⚠️ Algum valor em XS, XO ou YO não é finito\n")
                        return(NA)
                }
                
                XS0 <- XS[YS == 0, , drop = FALSE]
                XS1 <- XS[YS == 1, , drop = FALSE]
                YO[is.na(YO)] <- 0
                YO1 <- YO[YS == 1]
                XO1 <- XO[YS == 1, , drop = FALSE]
                
                XS0.g <- exp(as.numeric((XS0) %*% gamma))
                XS1.g <- exp(as.numeric((XS1) %*% gamma))
                XO1.b <- exp(as.numeric((XO1) %*% beta))
                
                term0 <- ((YO1 * (phi1 + 1)/(phi1 * XO1.b))^(1/2) - ((phi1 * XO1.b)/(YO1 * (phi1 + 1)))^(1/2))
                term1 <- exp((-phi1/4) * (term0)^2)
                term2 <- (((phi1 + 1)/(phi1 * XO1.b * YO1))^(1/2) + ((phi1 * XO1.b)/((phi1 + 1) * (YO1^3)))^(1/2))
                term3 <- (1/(2 * sqrt(2 * pi))) * ((phi1/2)^(1/2))
                term4 <- (((2)/(2 * XS1.g * (1 - rho^2)))^(1/2))
                term5 <- ((2 * XS1.g)/2) - 1
                term6 <- rho * (phi1/(2 * (1 - rho^2)))^(1/2)
                integrand <- term4 * term5 + term6 * term0
                term7 <- pnorm(integrand, log.p = TRUE)
                term8 <- ((1)^(1/2)) * (((2)/(1 * XS0.g))^(1/2) - ((1 * XS0.g)/2)^(1/2))
                FT2 <- pnorm(term8, log.p = TRUE)
                loglik <- sum(log(term2) + log(term1) + log(term3) + term7) + sum(FT2)
                
                return(-loglik)
        }
        
        gradlik_BS <- function(par) {
                rep(0, length(par))  # placeholder simples
        }
        
        fit <- optim(start, loglik_BS, gradlik_BS, method = "BFGS", hessian = TRUE, control = list(fnscale = 1))
        
        names(fit$par) <- c(colnames(XS), colnames(XO), "phi1", "rho")
        
        result <- list(
                Coefficients   = fit$par,
                loglik         = -fit$value,
                aic            = 2 * length(fit$par) - 2 * (-fit$value),
                bic            = -2 * (-fit$value) + length(fit$par) * log(nrow(XS)),
                counts         = fit$counts,
                hessian        = fit$hessian,
                prop_sigmaBS   = sqrt(diag(solve(fit$hessian))),
                nObs           = length(YS),
                N1             = sum(YS == 1),
                N0             = sum(YS == 0),
                message        = if (fit$convergence == 0) "Convergência bem-sucedida" else "⚠️ Não convergiu"
        )
        class(result) <- "HeckmanBS"
        return(result)
}


start_bs <- c(-1, 1, -20, 20, 0.5)  # ou ajuste conforme desejar

modelo_bs <- HeckmanBS_safe(
        selection = sel ~ log_dbh,
        outcome   = H_considering_damage ~ log_dbh,
        data      = df_clean,
        start     = start_bs
)
start_bs <- c(-1, 1, -20, 20, 0.5)
is.finite(start_bs)  # deve dar TRUE para todos


length(start_bs)       # deve ser 5
ncol(XS)               # quantos coeficientes na seleção
ncol(XO)               # quantos coeficientes na regressão
