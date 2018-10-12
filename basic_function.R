# fungsi custom yang bisa mendapatkan mean dari data.frame
# mean() dalam R belum tersedia fitur tersebut
mean_custom <- function(column) {
    if(is.data.frame(column)) {
        column <- column[, sapply(column, class) == 'numeric']
    }
    mean_vec <- sapply(
        column,
        function(x) Reduce('+', x) / length(x)
    )
    mean_vec
}

# mencari matriks kovarian, format output sama dengan cov()
cov_custom <- function(table) {
    if (!is.data.frame(table)) {
        stop('data class must be data.frame')
    }
    
    col_names <- colnames(table)
    nvar <- length(col_names)
    result <- matrix(data = NA, nrow = nvar, ncol = nvar, dimnames = list(col_names, col_names))
    for(i in col_names) {
        for (j in col_names) {
            var = var(table[i], table[j])
            result[i,j] = var
            if (i != j) result[j,i] = var
        }
    }
    result
}

onesample_z2 <- function(table, to_compare, mu0, sigma) {
    p <- length(to_compare)
    data <- table[,to_compare]
    n <- nrow(table)
    
    mean_vec <- mean_custom(data)
    diff <- mean_vec - mu0
    z2 <- n * t(diff) %*% solve(sigma) %*% diff
    p_val <- pchisq(z2, p, lower.tail = FALSE)
    cat('Z2 value: ', z2, ' ')
    cat('p-val Chisq: ', p_val, ' ')
    cat(if(p_val < 0.05) 'Reject Ho' else 'Fail to reject Ho', 'with alpha 0.05')
}

onesample_t2 <- function(table, to_compare, mu0) {
    p <- length(to_compare)
    data <- table[,to_compare]
    n <- nrow(data)
    df <- n - 1
    
    mean_vec <- mean_custom(data)
    cov_mtx <- cov_custom(data)
    diff <- mean_vec - mu0
    t2 <- n * t(diff) %*% solve(cov_mtx) %*% diff
    
    f_const <- (df-p+1) / (df*p)
    f_stat <- f_const * t2
    p_val <-  pf(f_stat, p, (df-p+1), lower.tail = FALSE)
    cat('T2 value is: ', t2, ' with df: ', df)
    cat('\n')
    cat('F p-val: ', p_val, ' ')
    cat(if(p_val < 0.05) 'Reject Ho' else 'Fail to reject Ho', 'with alpha 0.05')
}

twosample_t2 <- function(table, index, to_index, to_compare) {
    p <- length(to_compare)
    data_sub <- split(table, index)
    data_sub <- data_sub[to_index]
    simple_data <- lapply(data_sub, function(x) x[,to_compare])
    
    mean_vec <- lapply(simple_data, mean_custom)
    cov_mtx <- lapply(simple_data, cov_custom)
    n <- lapply(simple_data, nrow)
    N <- Reduce('+', n)
    df <- N - 2
    w_mtx <- mapply(function(n, cov) (n-1) * cov, n, cov_mtx, SIMPLIFY = FALSE)
    s_pool <- Reduce('+', w_mtx) / df
    y_diff <- Reduce('-', mean_vec)
    constant <- Reduce('*', n) / N
    t2 <- constant * t(y_diff) %*% solve(s_pool) %*% y_diff
    f_const <- (N-p-1) / ((N-2)*p)
    f_stat <- f_const * t2
    p_val <-  pf(f_stat, p, (N-p-1), lower.tail = FALSE)
    cat('T2 value: ', t2, ' with df: ', df)
    cat('\n')
    cat('F p-val: ', p_val, ' ')
    cat(if(p_val < 0.05) 'Reject Ho' else 'Fail to reject Ho', 'with alpha 0.05')
}

twosample_paired <- function(data1, data2, compare.col) {
    n <- nrow(data1)
    p <- length(compare.col)
    if (n != nrow(data2)) {
        stop('data1 and data2 must have same number of row')
    }
    
    diff_data <- data.frame(data1[compare.col] - data2[compare.col])
    mean_vec <- mean_custom(diff_data)
    cov_mtx <- cov_custom(diff_data)
    t2 <- n * t(mean_vec) %*% solve(cov_mtx) %*% mean_vec
    f_const <- (n-p) / ((n-1)*p)
    f_stat <- f_const * t2
    p_val <-  pf(f_stat, p, (n-p), lower.tail = FALSE)
    cat('T2 value: ', t2)
    cat('\n')
    cat('F p-val: ', p_val, ' ')
    cat(if(p_val < 0.05) 'Reject Ho' else 'Fail to reject Ho', 'with alpha 0.05')
}

#mode bisa t2 atau bonferroni
confint_smlt <- function(table, column, mode = 't2') {
    p <- length(column)
    n <- nrow(table)
    sub_data <- table[column]
    mean_vec <- mean_custom(sub_data)
    cov_mtx <- cov_custom(sub_data)
    
    param <- list()
    switch (mode,
        t2 = {
            f_stat <- qf(0.05, p, n, lower.tail = FALSE)
            t2_const <- sqrt(p * (n-1) * f_stat / (n - p))
            param['const'] <- t2_const
        },
        bonferroni = {
            t_const <- qt(0.05 / (2*p), n - 1, lower.tail = FALSE)
            param['const'] <- t_const
        }
    )
    
    cat('mu with', mode, 'method', sep = ' ')
    cat('\n')
    cat('variable','lower_bound','upper_bound', sep = '\t\t')
    cat('\n')
    for (i in 1:p) {
        a <- numeric(p)
        a[i] <- 1
        variability <- param$const * sqrt(t(a) %*% cov_mtx %*% a) * sqrt(1/n)
        mean_val <- t(a) %*% mean_vec
        cat(column[i], mean_val - variability, mean_val + variability, sep = '\t')
        cat('\n')
    }
}

oneway_manova <- function(data, index, column) {
    p <- nlevels(index)
    if (p != 3) {
        stop('number of groups must be three')
    }
    
    xbar <- mean_custom(data[,column])
    n <- nrow(data)
    grouped_data <- split(data, index)
    grouped_data <- lapply(grouped_data, function(x) x[,column])
    
    ng <- sapply(grouped_data, nrow)
    xg_bar <- lapply(grouped_data, mean_custom)
    cov.g <- lapply(grouped_data, cov_custom)
    
    b_list <- mapply(
        function(n, xi.bar) n * (xi.bar-xbar) %*% t(xi.bar-xbar),
        ng, xg_bar, SIMPLIFY = FALSE
    )
    b_mtx <- Reduce('+', b_list)
    w_list <- mapply(function(n, s) (n-1) * s, ng, cov.g, SIMPLIFY = FALSE)
    w_mtx <- Reduce('+', w_list)
    
    wilks_lmd <- det(w_mtx) / det(w_mtx + b_mtx)
    f_stat <- ((n-p-2) / p) * ((1-sqrt(wilks_lmd)) / sqrt(wilks_lmd))
    p_val <- pf(f_stat, 2*p, 2*(n-p-2), lower.tail = FALSE)
    cat('F-stat value:',f_stat,'with p-val:',p_val,sep = ' ')
    cat('\n')
    cat('set alpha 0.05 then',if(p_val < 0.05) 'reject Ho' else 'fail to reject Ho')
}