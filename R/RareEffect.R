# Read Group file and split variants by functional annotations
# groups: character vector of annotation labels to look for in the group file
#         e.g. c("lof", "missense", "synonymous")
read_groupfile <- function(groupfile_name, gene_name, groups) {
    groupfile <- file(groupfile_name, "r")
    line <- 0
    var <- NULL
    anno <- NULL

    while (TRUE) {
        line <- line + 1
        marker_group_line <- readLines(groupfile, n = 1)

        if (length(marker_group_line) == 0) {
            break
        }

        marker_group_line_list <- strsplit(marker_group_line, split = "[ \t]+")[[1]]

        if (marker_group_line_list[1] == gene_name) {
            if (marker_group_line_list[2] == "var") {
                var <- marker_group_line_list
            } else {
                anno <- marker_group_line_list
            }
        }
    }

    # Build a named list: one element per group with the matching variants
    out <- list()
    for (g in groups) {
        idx <- which(anno == g)
        out[[g]] <- var[idx]
    }

    close(groupfile)
    return(out)
}

read_matrix_by_one_marker <- function(objGeno, var_list, sampleID) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        print("No variants in the list")
        return(mat)
    }
    is_flipped <- rep(FALSE, length(var_list))
    for (i in 1:length(var_list)) {
        t_GVec <- rep(0, n_samples)
        idx <- which(var_list[i] == objGeno$markerInfo$ID)
        SAIGE::Unified_getOneMarker(t_genoType = objGeno$genoType,
            t_gIndex_prev = objGeno$markerInfo$genoIndex_prev[idx],
            t_gIndex = objGeno$markerInfo$genoIndex[idx],
            t_ref = "2",
            t_alt = "1",
            t_marker = objGeno$markerInfo$ID[idx],
            t_pd = objGeno$markerInfo$POS[idx],
            t_chr = toString(objGeno$markerInfo$CHROM[idx]),
            t_altFreq = 0,
            t_altCounts = 0,
            t_missingRate = 0,
            t_imputeInfo = 0,
            t_isOutputIndexForMissing = TRUE,
            t_indexForMissing = 0,
            t_isOnlyOutputNonZero = FALSE,
            t_indexForNonZero = 0,
            t_GVec = t_GVec,
            t_isImputation = FALSE
        )

	t_GVec[t_GVec < 0] <- 0     # Convert missing to zero
	t_GVec[t_GVec > 2] <- 0     # Sometimes missing value is a very large positive number
        # If MAF > 0.5, flip allele
        if (sum(t_GVec) > n_samples) {
            t_GVec <- 2 - t_GVec
            is_flipped[i] <- TRUE
        }
        t_GVec_sp <- as(t_GVec, "sparseVector")
        t_GVec_sp_mat <- as(t_GVec_sp, "Matrix")

        mat <- cbind(mat, t_GVec_sp_mat)
    }

    colnames(mat) <- var_list
    rownames(mat) <- sampleID

    return(list(mat, is_flipped))
}

collapse_matrix <- function(objGeno, var_list, sampleID, modglmm, macThreshold = 10, rare_maf_threshold = 0.01) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        print("No variants in the list")
        return(list(mat, NULL))
    }
    tmp <- read_matrix_by_one_marker(objGeno, var_list, sampleID)
    mat <- tmp[[1]]
    is_flipped <- tmp[[2]]
    flipped_var <- as.data.table(cbind(var_list, is_flipped))
    mat <- mat[which(rownames(mat) %in% modglmm$sampleID), , drop = FALSE]
    MAF <- colSums(mat) / (2 * nrow(mat))
    idx_rare <- which((MAF < rare_maf_threshold) & (colSums(mat) >= macThreshold))
    idx_UR <- which((MAF > 0) & (colSums(mat) < macThreshold))
    # mat <- mat[, idx, drop = FALSE]
    if (length(idx_rare) == 0 & length(idx_UR) == 0) {
        print(paste0("No variants with MAF < ", rare_maf_threshold, " or MAF > 0.99"))
        # mat <- mat[, idx, drop = FALSE]
        mat <- Matrix::Matrix(ncol=0, nrow=0)
        return(list(mat, NULL))
    }
    mat_rare <- mat[, idx_rare, drop = FALSE]
    mat_UR <- mat[, idx_UR, drop = FALSE]
    UR_rowsum <- rowSums(mat_UR)
    UR_rowsum[which(UR_rowsum > 1)] <- 1
    mat_UR_collapsed <- cbind(
	if (!is.null(mat_rare) && ncol(mat_rare) > 0) mat_rare,
	if (!is.null(UR_rowsum) && length(UR_rowsum) > 0) UR_rowsum
	)
    # mat_UR_collapsed <- cbind(mat_rare, UR_rowsum)
    return(list(mat_UR_collapsed, flipped_var))
}

get_range <- function(v) {
    # input is a named list of variant vectors by functional annotation
    n_groups <- length(v)
    start_pos <- Inf
    end_pos <- 0
    for (i in 1:n_groups) {
        if (length(v[[i]]) > 0) {
            pos <- as.numeric(str_split_fixed(v[[i]], ":", 4)[,2]) # Modified in 0.3
            sub_start_pos <- min(pos)
            sub_end_pos <- max(pos)
            print(sub_start_pos)
            print(sub_end_pos)
            if (start_pos > sub_start_pos) {
                start_pos <- sub_start_pos
            }
            if (end_pos < sub_end_pos) {
                end_pos <- sub_end_pos
            }
        }
    }
    return (list(start_pos, end_pos))
}


# Read phenotype file
read_pheno <- function(pheno_file, pheno_code, iid_col = "f.eid") {
    pheno <- data.table::fread(pheno_file, quote = "")
    out <- subset(pheno, select = c(iid_col, pheno_code))
    return(out)
}


calc_log_lik <- function(delta, S, UtY, Y_UUtY) {
    k <- length(S)
    n <- nrow(Y_UUtY)

    log_lik1 <- 0
    for (i in 1:k) {
        log_lik1 <- log_lik1 + (UtY[i, ])^2 / (S[i] + delta)
    }

    log_lik2 <- 1 / delta * sum((Y_UUtY)^2)

    out <- -0.5 * (n * log(2 * pi) + sum(log(S + delta)) + (n - k) * log(delta)
                   + n + n * log(1 / n * (log_lik1 + log_lik2)))

    return(as.numeric(out))
}

calc_post_beta <- function(K, G, delta, S, UtY, U) {
    K_sparse <- as(K, "dgCMatrix")
    if (length(S) == 1) {
        S <- as.matrix(S)
    }

    out <- K_sparse %*% t(G) %*% U %*% diag(1 / (S + delta)) %*% (UtY)
    return(out)
}

# Run FaST-LMM to obtain posterior beta
fast_lmm <- function(G, Y) {
    # if sum(G) == 0, let effect size = 0
    print("Estimating beta using FaST-LMM")
    G[is.na(G)] <- 0
    GtG <- t(G) %*% G
    if (sum(G, na.rm = T) == 0) {
        return (list(as.matrix(0), 0, 1e6, GtG))
    }
    Y <- as.matrix(Y)
    K <- diag(1, nrow = ncol(G))
    L <- chol(K)
    W <- G %*% L
    W_sparse <- as(W, "dgCMatrix")
    svd_mat <- sparsesvd::sparsesvd(W_sparse)

    U <- svd_mat$u
    S <- (svd_mat$d)^2

    UtY <- t(U) %*% Y

    # Y_UUtY = Y - UUtY
    Y_UUtY <- Y - U %*% UtY

    opt <- optim(par = 1, fn = calc_log_lik, S = S, UtY = UtY, Y_UUtY = Y_UUtY,
                 method = c("Brent"), lower = 0, upper = 1e6, control = list(fnscale = -1))
    opt_delta <- opt$par

    tr_GtG <- sum(diag(GtG %*% K))

    post_beta <- calc_post_beta(K, G, opt_delta, S, UtY, U)

    return (list(post_beta, tr_GtG, opt_delta, GtG))
}


calc_gene_effect_size <- function(G, lof_ncol, post_beta, Y) {
    beta_lof <- post_beta[1:lof_ncol]
    beta_lof <- abs(beta_lof)
    lof_prs <- G[, 1:lof_ncol, drop = F] %*% beta_lof
    lof_prs_norm <- (lof_prs - mean(lof_prs)) / sd(lof_prs)
    m1 <- lm(Y ~ lof_prs_norm[, 1])

    return(m1$coefficients[2])
}

mom_estimator_marginal <- function(G, y) {
    n <- length(y)
    G <- as(G, "dgCMatrix")
    G[is.na(G)] <- 0
    # tr(G Sigma G^T G Sigma G^T) = sum((G Sigma G^T )^2) = sum(G^T G Sigma)^2
    Sigma <- diag(1, nrow = ncol(G))

    system.time({
        t1 <- sum((crossprod(G) %*% Sigma)^2)
    })
    t2 <- sum(diag(crossprod(G) %*% Sigma))
    A <- matrix(c(t1, t2, t2, n), ncol = 2)
    c1 <- as.numeric(t(y) %*% G %*% Sigma %*% t(G) %*% y)
    c2 <- sum(y^2)
    b <- matrix(c(c1, c2), ncol = 1)
    var_comp <- solve(A) %*% b
    return (var_comp)
}

# Generalized MoM estimator for N groups
# G_list: list of genotype matrices (one per group)
# y: phenotype residual vector
mom_estimator_joint <- function(G_list, y) {
    n <- length(y)
    n_groups <- length(G_list)

    # Prepare sparse matrices and covariance components
    G_sparse <- list()
    Sigma <- list()
    L_mat <- list()
    for (i in 1:n_groups) {
        G_sparse[[i]] <- as(G_list[[i]], "dgCMatrix")
        G_sparse[[i]][is.na(G_sparse[[i]])] <- 0
        Sigma[[i]] <- diag(1, nrow = ncol(G_sparse[[i]]))
        L_mat[[i]] <- chol(Sigma[[i]])
    }

    # Build the (n_groups+1) x (n_groups+1) system of equations
    dim_A <- n_groups + 1
    A <- matrix(0, nrow = dim_A, ncol = dim_A)

    # Diagonal: t_ii = sum((G_i^T G_i Sigma_i)^2)
    for (i in 1:n_groups) {
        A[i, i] <- sum((crossprod(G_sparse[[i]]) %*% Sigma[[i]])^2)
    }

    # Off-diagonal: t_ij = sum((t(G_i L_i) %*% (G_j L_j))^2)
    if (n_groups > 1) {
        for (i in 1:(n_groups - 1)) {
            for (j in (i + 1):n_groups) {
                val <- sum((t(G_sparse[[i]] %*% L_mat[[i]]) %*% (G_sparse[[j]] %*% L_mat[[j]]))^2)
                A[i, j] <- val
                A[j, i] <- val
            }
        }
    }

    # Last row/column: t_i,(n_groups+1) = sum(diag(G_i^T G_i Sigma_i))
    for (i in 1:n_groups) {
        val <- sum(diag(t(G_sparse[[i]]) %*% G_sparse[[i]] %*% Sigma[[i]]))
        A[i, dim_A] <- val
        A[dim_A, i] <- val
    }
    A[dim_A, dim_A] <- n

    # Right-hand side
    b_vec <- numeric(dim_A)
    for (i in 1:n_groups) {
        b_vec[i] <- as.numeric(t(y) %*% G_sparse[[i]] %*% Sigma[[i]] %*% t(G_sparse[[i]]) %*% y)
    }
    b_vec[dim_A] <- sum(y^2)
    b <- matrix(b_vec, ncol = 1)

    var_comp <- solve(A) %*% b
    return (var_comp)
}

# Generalized joint BLUP for N groups
# G_list: list of genotype matrices
# tau_list: list of tau values (one per group)
# psi: residual variance (sigma_sq)
# y: phenotype residual vector
calculate_joint_blup <- function(G_list, tau_list, psi, y) {
    n_groups <- length(G_list)
    G_sparse_list <- list()
    for (i in 1:n_groups) {
        G_sparse_list[[i]] <- as(G_list[[i]], "dgCMatrix")
    }
    G <- do.call(cbind, G_sparse_list)
    G[is.na(G)] <- 0

    # Build block-diagonal Sigma = bdiag(tau1*I1, tau2*I2, ...)
    sigma_blocks <- list()
    for (i in 1:n_groups) {
        sigma_blocks[[i]] <- tau_list[[i]] * diag(1, ncol(G_list[[i]]))
    }
    Sigma <- bdiag(sigma_blocks)
    beta <- solve(t(G) %*% G / psi + solve(Sigma)) %*% t(G) %*% y / psi
    return (beta)
}

weight_cal <- function(beta_k, delta = 10^(-5), gamma = 2, q = 0, factor2 = 1, n_per_group){
	# delta=10^(-5); gamma=2; q=0
	p <- length(beta_k)
	w_out = rep(0, p)
	idx_null = NULL
	idx<-which(abs(beta_k) <= delta)
	if(length(idx)> 0){
		beta_k1<-beta_k[idx]
		b_delta = abs(beta_k1/delta)^{gamma}
		logw1 = (q-2) * log(delta) + (q-2)/gamma* log(1+ b_delta)
		w_out[idx]<-exp(logw1 * factor2)
		idx_null = idx
	}
	idx<-which(abs(beta_k) > delta)
	if(length(idx)> 0){
		beta_k1<-beta_k[idx]
		b_delta_inv = abs(delta/beta_k1)^{gamma}
		logw2 = (q-2) * log(abs(beta_k1)) + (q-2)/gamma* log(1+ b_delta_inv)
		w_out[idx]<-exp(logw2)
	}

	# make sum of all to be p
	a1 = 1/w_out
	a1_out = a1/sum(a1)*p
	w_out = 1/a1_out

	return(list(w_out=w_out))
}

adaptive_ridge <- function(X, y, lambda, q = 0, delta = 1e-5, gamma = 2, max_iter = 100, tol = 0.01, sigma_sq = 1, n_per_group) {
    w <- rep(1, ncol(X))     # Initialize weights
    beta <- rep(0, ncol(X))  # Initialize beta
    beta_all<-NULL
    w_all<-NULL
    idx_set<-1:length(w)

    # save to reduce redundant computation
    Xy = t(X) %*% y
    XX = t(X) %*% X
    for (k in 1:max_iter) {
        W <- diag(w)
        # using solve is better than calculating inverse matrix first...
        beta_new <- solve(XX + lambda * W, Xy)  # Ridge regression for initial estimate

        # Update weights
        w_new_re <- weight_cal(beta_new, delta=delta, gamma=gamma, q=q, n_per_group=n_per_group)
        w_new = w_new_re$w_out

        # print(sqrt(sum((beta_new - beta)^2)))
        if (sqrt(sum((beta_new - beta)^2)) < sum(beta^2)*tol) {
            # print(paste0("Converged at iteration: ", k))
            break
        }

        beta <- beta_new
        w <- w_new
        w_all<-cbind(w_all, w)
        beta_all<-cbind(beta_all, beta)
    }

    A <- XX
    diag(A) <- diag(A) + lambda * w
    L <- chol(A)
    A_inv <- chol2inv(L)

    cov_beta <- as.numeric(sigma_sq) * (A_inv %*% XX %*% A_inv)
    se_beta <- sqrt(diag(cov_beta))

    list(beta = beta, weights = w, iterations = k, beta_all=beta_all, w_all=w_all, se_beta = se_beta)
}

# groups: a comma-separated string of annotation group labels, e.g. "lof,missense,synonymous"
#         These must match the annotation labels in the group file.
# collapseGroups: a named logical vector indicating which groups to collapse,
#                 e.g. c(lof=TRUE, missense=FALSE, synonymous=FALSE)
#                 If NULL, no groups are collapsed.
run_RareEffect <- function(rdaFile, chrom, geneName, groupFile, traitType, bedFile, bimFile, famFile, macThreshold, collapseGroups = NULL, apply_AR, outputPrefix, groups = "lof,missense,synonymous") {
    # Parse groups from comma-separated string
    group_names <- trimws(strsplit(groups, ",")[[1]])
    n_groups <- length(group_names)

    # Load SAIGE step 1 results
    load(rdaFile)

    if (traitType == "binary") {
        v <- modglmm$fitted.values * (1 - modglmm$fitted.values)    # v_i = mu_i * (1 - mu_i)
        modglmm$residuals <- sqrt(v) * (1 / (modglmm$fitted.values * (1 - modglmm$fitted.values)) * (modglmm$y - modglmm$fitted.values))
    }
    sigma_sq <- var(modglmm$residuals)

    # Set PLINK object
    bim <- fread(bimFile)
    fam <- fread(famFile)
    sampleID <- as.character(fam$V2)
    n_samples <- length(sampleID)

    objGeno <- SAIGE::setGenoInput(bgenFile = "",
            bgenFileIndex = "",
            vcfFile = "",
            vcfFileIndex = "",
            vcfField = "",
            savFile = "",
            savFileIndex = "",
            sampleFile = "",
            bedFile=bedFile,
            bimFile=bimFile,
            famFile=famFile,
            idstoIncludeFile = "",
            rangestoIncludeFile = "",
            chrom = "",
            AlleleOrder = "alt-first",
            sampleInModel = sampleID
        )

    print("Analysis started")

    var_by_func_anno <- read_groupfile(groupFile, geneName, group_names)

    # Print which groups have no variants
    for (g in group_names) {
        if (length(var_by_func_anno[[g]]) == 0) {
            print(paste0(g, " variant does not exist."))
        }
    }

    # Remove variants not in plink file
    for (g in group_names) {
        var_by_func_anno[[g]] <- var_by_func_anno[[g]][which(var_by_func_anno[[g]] %in% objGeno$markerInfo$ID)]
    }

    # Read and collapse genotype matrices for each group
    mat_collapsed <- list()
    mat_flipped <- list()
    for (g in group_names) {
        should_collapse <- FALSE
        if (!is.null(collapseGroups) && g %in% names(collapseGroups)) {
            should_collapse <- collapseGroups[[g]]
        }

        if (should_collapse) {
            tmp <- collapse_matrix(objGeno, var_by_func_anno[[g]], sampleID, modglmm, macThreshold = 0)
            mat_collapsed[[g]] <- tmp[[1]]
            mat_flipped[[g]] <- tmp[[2]]
            mat_collapsed[[g]] <- Matrix::Matrix(rowSums(mat_collapsed[[g]]), ncol = 1, sparse = TRUE, dimnames = list(rownames(mat_collapsed[[g]]), NULL))
        } else {
            tmp <- collapse_matrix(objGeno, var_by_func_anno[[g]], sampleID, modglmm, macThreshold)
            mat_collapsed[[g]] <- tmp[[1]]
            mat_flipped[[g]] <- tmp[[2]]
        }
    }

    # Set column name for UR column
    for (g in group_names) {
        if (ncol(mat_collapsed[[g]]) > 0) {
            colnames(mat_collapsed[[g]])[ncol(mat_collapsed[[g]])] <- paste0(g, "_UR")
        }
    }

    # Make combined genotype matrix
    nonempty_mat <- list()
    for (g in group_names) {
        if (ncol(mat_collapsed[[g]]) > 0) {
            nonempty_mat <- c(nonempty_mat, list(mat_collapsed[[g]]))
        }
    }

    G <- do.call(cbind, nonempty_mat)

    # Track number of columns per group
    group_ncol <- list()
    for (g in group_names) {
        group_ncol[[g]] <- ncol(mat_collapsed[[g]])
    }

    # Obtain residual vector and genotype matrix with the same order
    y_tilde <- cbind(modglmm$sampleID, modglmm$residuals)
    y_tilde <- y_tilde[which(y_tilde[,1] %in% sampleID),]
    G_reordered <- G[match(y_tilde[,1], rownames(G)),]
    if (traitType == "binary") {
        vG_reordered <- as.vector(sqrt(v)) * G_reordered    # Sigma_e^(-1/2) G
    }
    n_samples <- nrow(G_reordered)

    # Compute column offsets for each group
    col_start <- list()
    col_end <- list()
    offset <- 0
    for (g in group_names) {
        col_start[[g]] <- offset + 1
        col_end[[g]] <- offset + group_ncol[[g]]
        offset <- offset + group_ncol[[g]]
    }

    # Define matrices by functional annotation
    post_beta_marginal <- list()
    group_mat <- list()
    for (g in group_names) {
        if (group_ncol[[g]] == 0) {
            group_mat[[g]] <- NULL
            post_beta_marginal[[g]] <- NULL
        } else {
            cols <- col_start[[g]]:col_end[[g]]
            if (traitType == "binary") {
                group_mat[[g]] <- vG_reordered[, cols, drop = F]
            } else {
                group_mat[[g]] <- G_reordered[, cols, drop = F]
            }
        }
    }

    # Run FaST-LMM for each group (estimate marginal variance component tau)
    tau <- list()
    tau_mom_marginal <- list()
    tr_GtG <- list()
    delta_group <- list()
    GtG_group <- list()
    fast_lmm_results <- list()

    for (g in group_names) {
        if (group_ncol[[g]] > 0) {
            fl <- fast_lmm(G = group_mat[[g]], Y = as.numeric(y_tilde[,2]))
            post_beta_marginal[[g]] <- fl[[1]]
            tr_GtG[[g]] <- fl[[2]]
            delta_group[[g]] <- fl[[3]]
            GtG_group[[g]] <- fl[[4]]
            tau[[g]] <- as.numeric(sigma_sq / delta_group[[g]])
            tau_mom_marginal[[g]] <- mom_estimator_marginal(G = group_mat[[g]], y = as.numeric(y_tilde[,2]))
        } else {
            tau[[g]] <- 0
        }
    }

    # Estimate variance component tau jointly (MoM estimator), only if all groups are non-empty
    all_nonempty <- all(sapply(group_names, function(g) group_ncol[[g]] > 0))

    tau_adj <- list()
    if (all_nonempty) {
        G_list_for_mom <- lapply(group_names, function(g) group_mat[[g]])
        tau_mom_joint <- mom_estimator_joint(G_list = G_list_for_mom, y = as.numeric(y_tilde[,2]))

        # Adjust variance component tau
        for (i in 1:n_groups) {
            g <- group_names[i]
            tau_adj[[g]] <- tau[[g]] * tau_mom_joint[i] / tau_mom_marginal[[g]][1]
        }
    } else {
        # Use marginal variance component tau if any group is empty
        for (g in group_names) {
            tau_adj[[g]] <- tau[[g]]
        }
    }

    # Estimate gene-level heritability
    h2_adj <- list()
    for (g in group_names) {
        if (group_ncol[[g]] > 0) {
            if (traitType == "binary") {
                h2_adj[[g]] <- max(tau_adj[[g]] * tr_GtG[[g]] / (tau_adj[[g]] * tr_GtG[[g]] + sigma_sq * sum(1/v)), 0)
            } else {
                h2_adj[[g]] <- max(tau_adj[[g]] * tr_GtG[[g]] / (tau_adj[[g]] * tr_GtG[[g]] + sigma_sq * n_samples), 0)
            }
        } else {
            h2_adj[[g]] <- 0
        }
    }

    # Obtain effect size jointly if all groups have tau_adj > 0
    all_tau_positive <- all(sapply(group_names, function(g) tau_adj[[g]] > 0))

    if (all_tau_positive) {
        G_list_for_blup <- lapply(group_names, function(g) group_mat[[g]])
        tau_list_for_blup <- lapply(group_names, function(g) tau_adj[[g]])
        post_beta <- calculate_joint_blup(
            G_list = G_list_for_blup,
            tau_list = tau_list_for_blup,
            psi = as.numeric(sigma_sq),
            y = as.numeric(y_tilde[,2])
        )
    } else {
        # If not, calculate beta marginally
        post_beta <- do.call(rbind, lapply(group_names, function(g) post_beta_marginal[[g]]))
        post_beta <- as.vector(post_beta)
    }

    post_beta <- as.vector(post_beta)

    if (apply_AR == TRUE) {
        # Build lambda vector from per-group tau_adj and group sizes
        lambda_parts <- list()
        n_per_group <- list()
        for (g in group_names) {
            sigma_g <- rep(1, ncol(group_mat[[g]]))
            lambda_parts <- c(lambda_parts, list(as.numeric(sigma_sq) / (tau_adj[[g]] * sigma_g)))
            n_per_group[[g]] <- group_ncol[[g]]
        }
        lambda <- do.call(c, lambda_parts)
        result_AR <- adaptive_ridge(G_reordered, as.numeric(y_tilde[,2]), lambda, q = 0, delta = 1e-5, gamma = 2, max_iter = 5, tol = 0.01, sigma_sq = sigma_sq, n_per_group = n_per_group)
        post_beta <- as.vector(result_AR$beta)
        se_beta <- (as.vector(result_AR$se_beta))
    }

    if (apply_AR == FALSE) {
        # Obtain prediction error variance (PEV)
        PEV_list <- list()
        for (g in group_names) {
            if (group_ncol[[g]] > 0) {
                GtG_g <- GtG_group[[g]] / as.numeric(sigma_sq)
                diag(GtG_g) <- diag(GtG_g) + 1 / tau_adj[[g]]
                PEV_list[[g]] <- diag(solve(GtG_g))
            } else {
                PEV_list[[g]] <- NULL
            }
        }
        PEV <- abs(do.call(c, PEV_list))
    } else {
        PEV <- se_beta^2
    }

    # Apply Firth bias correction for binary phenotype
    if (traitType == "binary") {
        y_binary <- cbind(modglmm$sampleID, modglmm$y)
        offset1 <- modglmm$linear.predictors - modglmm$coefficients[1]
        l2.var = 1
        maxit = 50

        beta_firth_parts <- list()
        for (g in group_names) {
            if (group_ncol[[g]] > 0) {
                cols <- col_start[[g]]:col_end[[g]]
                G_g_sp <- as(G_reordered[, cols, drop = F], "sparseMatrix")
                nMarker_g <- ncol(G_g_sp)
                firth_result <- Run_Firth_MultiVar_Single(G_g_sp, modglmm$obj.noK, as.numeric(y_binary[,2]), offset1, nMarker_g, l2.var=1/(2*tau_adj[[g]]), Is.Fast=FALSE, Is.Sparse=TRUE)[,2]
                beta_firth_parts <- c(beta_firth_parts, list(firth_result))
            }
        }

        beta_firth <- do.call(c, beta_firth_parts)
        effect <- ifelse(abs(post_beta) < log(2), post_beta, beta_firth)
    } else {
        effect <- post_beta
    }

    effect <- as.vector(effect)
    # output related to single-variant effect size
    variant <- colnames(G)

    # Sign of the effect: use the first group if it has exactly 1 variant,
    # otherwise use a weighted sum approach
    first_group <- group_names[1]
    first_ncol <- group_ncol[[first_group]]
    if (first_ncol == 1) {
        sgn <- sign(post_beta[1])
    } else if (first_ncol == 0) {
        print(paste0(first_group, " variant does not exist, so the sign of the effect size is calculated by the sign of other groups."))
        MAC <- colSums(G_reordered)
        sgn <- sign(sum(MAC * post_beta))
    } else {
        MAC <- colSums(G_reordered[, c(1:first_ncol), drop = F])
        sgn <- sign(sum(MAC * post_beta[c(1:first_ncol)]))
    }

    h2_values <- sapply(group_names, function(g) h2_adj[[g]])
    h2_all <- sum(h2_values) * sgn
    h2 <- c(h2_values, h2_all)
    group_labels <- c(group_names, "all")

    effect_out <- as.data.frame(cbind(variant, effect, PEV))
    effect_out$effect <- as.numeric(effect_out$effect)
    effect_out$PEV <- as.numeric(effect_out$PEV)
    h2_out <- rbind(group_labels, h2)

    # Find flipped var and change the sign of the effect size
    flipped_var_all <- do.call(rbind, lapply(group_names, function(g) mat_flipped[[g]]))
    flipped_var <- flipped_var_all[is_flipped == TRUE,]$var_list
    effect_out[which(effect_out$variant %in% flipped_var),]$effect <- -effect_out[which(effect_out$variant %in% flipped_var),]$effect

    effect_outname <- paste0(outputPrefix, "_effect.txt")
    h2_outname <- paste0(outputPrefix, "_h2.txt")

    write.table(effect_out, effect_outname, row.names=F, quote=F)
    write.table(h2_out, h2_outname, row.names=F, col.names=F, quote=F)

    print("Analysis completed.")
}
