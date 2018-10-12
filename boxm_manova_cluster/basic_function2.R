# data("iris")
# head(iris)
# iris$new <- ifelse(iris$Sepal.Length > 5, 'a', 'b')
# iris$new <- as.factor(iris$new)

# my_data <- read.csv('/storage/Code/data-r/apg/manova2_test.csv')
# my_data$pend_krt <- as.factor(my_data$pend_krt)
# my_data$jumlah_art <- as.factor(my_data$jumlah_art)

twoway_anova <- function (respon,faktora,faktorb, alfa) {
    
    n.a <- tapply(respon, faktora, length)
    a   <- length(n.a) #jumlah faktor a
    n.b <- tapply(respon, faktorb, length)
    b   <- length(n.b) #jumlah faktor b
    n.t <- length(respon)/(a*b) #jumlah n, n nya harus sama semua di setiap katagori
    
    sum.t <- (sum(respon)^2)/(a*b*n.t) #T***^2 / abn
    SST   <- (sum(respon[]^2))-sum.t #cari jkt sigma^3 x^2
    sum.a <- tapply(respon, faktora, sum) 
    SSA   <- (sum(sum.a[]^2)/(b*n.t))-sum.t #sigma T*j*^2/bn - T***^2/abn
    sum.b <- tapply(respon, faktorb, sum)
    SSB   <- (sum(sum.b[]^2)/(a*n.t))-sum.t #sigma Ti**^2/an - T***^2/abn
    sum.p <- tapply(respon, list(faktora,faktorb),sum)
    SSP   <- (sum(sum.p[]^2)/n.t)-sum.t 
    SSAB  <- SSP-SSA-SSB
    SSE   <- SST-SSP
    
    MSP  <- SSP/((a*b)-1) #SSP/db
    MSA  <- SSA/(a-1) #SSA/db
    MSB  <- SSB/(b-1) #SSB/db
    MSAB <- SSAB/((a-1)*(b-1)) #SSAB/db
    MSE  <- SSE/(a*b*(n.t-1)) #SSE/db
    
    #Statistik Uji
    fhit.a  <- MSA/MSE # F hitung faktor a
    fhit.b  <- MSB/MSE # F hitung faktor b
    fhit.ab <- MSAB/MSE # F hitung inteaksi faktor a dan b
    pval.a <- pf(fhit.a, df1=a-1, df2=a*b*(n.t-1)) # p-value faktor a
    pval.b <- pf(fhit.b, df1=b-1, df2=a*b*(n.t-1)) # p-value faktor b
    pval.ab<- pf(fhit.ab, df1=(a-1)*(b-1), df2=a*b*(n.t-1)) # p-value interaksi faktor a dan b
    
    #Bikin Tabel
    Sumberkeragaman <- c(" Faktor A", " Faktor B", " FaktorA:B", "Error", "Total")
    SS <- c(SSA,SSB,SSAB,SSE,SST) #Sum Square
    df <- c(a-1, b-1, (a-1)*(b-1), a*b*(n.t-1), (a*b*n.t)-1) #drajat bebas
    MS <- c(MSA,MSB,MSAB,MSE,NA) # Mean Square
    Fhit <- c(fhit.a, fhit.b, fhit.ab, NA, NA) #F hitung
    Pval <- c(pval.a, pval.b, pval.ab, NA, NA) #p-value
    
    tabelanova <- data.frame(Sumberkeragaman,df,SS,MS,Fhit,Pval)
    
    cat("ANOVA DUA FAKTOR \n")
    cat("1. Hipotesis : \n")
    cat("    H0 : (AB)ij = 0 \n")
    cat("    H0 : A1= ... = A",a,"\n")
    cat("    H0 : B1= ... = B",b,"\n")
    cat("\n")
    cat("2. Alpha = ",alfa,"\n")
    cat("\n")
    cat("3. ANOVA DUA FAKTOR (dengan replikasi, dengan interaksi, model tetap) \n")
    print(tabelanova)
    cat("\n")
    cat("Kesimpulan: \n")
    if (pval.ab>alfa) {
        cat("1. Gagal tolak H0. Cukup bukti bahwa tidak ada interaksi antara faktor A dan Faktor B.")
    } else {
        cat("1. Tolak H0. Cukup bukti bahwa ada interaksi antara faktor A dengan faktor B.")
    }
    
    cat("\n")
    if (pval.a>alfa) {
        cat("2. Gagal tolak H0. Cukup bukti bahwa tidak ada efek dari perbedaan faktor A.")
    } else {
        cat("2. Tolak H0. Cukup bukti bahwa ada efek dari perbedaan faktor A")
    }
    
    cat("\n")
    if (pval.b>alfa) {
        cat("3. Gagal tolak H0. Cukup bukti bahwa tidak ada efek dari perbedaan faktor B.")
    } else {
        cat("3. Tolak H0. Cukup bukti bahwa ada efek dari perbedaan faktor B")
    }
    cat("\n")
    cat("=========================================================================================")
    
}


box_m <- function(table, resp, group) {
    n_tot <- nrow(table)
    p <- length(resp)
    g <- nlevels(group)
    group_dat <- split(table, group)
    group_dat <- lapply(group_dat, function(x) x[,resp])
    
    n_list <- lapply(group_dat, nrow)
    S_list <- lapply(group_dat, cov)
    temp_list <- mapply(function(n, S) (n-1) * S, n_list, S_list,
                        SIMPLIFY = F)
    S_pool <- (1 / (n_tot-g)) * Reduce('+', temp_list)
    pool_det <- det(S_pool)
    S_det <- mapply(function(x, n) (det(x)/pool_det)^((n-1)/2),
                    S_list, n_list)
    wilks_lmd <- prod(S_det)
    M <- -2 * log(wilks_lmd)
    u <- (sum(sapply(n_list, function(n) 1/(n-1))) 
          - 1/(n_tot-g)) * ((2*p^2 + 3*p -1)/(6*(p+1)*(g-1)))
    C = (1-u) * M
    df = 0.5 * p * (p+1) * (g-1)
    p_val <- pchisq(C, df, lower.tail = FALSE)
    cat('Wilks', wilks_lmd, sep = '\t')
    cat('\n')
    cat('M', M, sep = '\t')
    cat('\n')
    cat('u', u, sep = '\t')
    cat('\n')
    cat('p-val Chisq: ', p_val, ' ')
    cat(if(p_val <= 0.05) 'Reject Ho' else 'Fail to reject Ho', 'with alpha 0.05')
    cat('\n')
    
    list(wilks = wilks_lmd, Mval = M, uval = u, pval = p_val)
}

twoway_manova <- function(table, resp, group1, group2) {
    #check if factor
    if(!is.factor(group1) | !is.factor(group2)) {
        stop('grouping variable is not a factor')
    }
    
    g <- nlevels(group1)
    b <- nlevels(group2)
    p <- length(resp)
    n <- nrow(table) / (g*b)
    xbar <- apply(table[,resp], 2, mean)
    
    tmp <- split(table, group1)
    xbar_l <- sapply(tmp, function(x) colMeans(x[,resp]))
    diff_xl <- sweep(xbar_l, 1, xbar)
    sum_g1 <- b*n * diff_xl %*% t(diff_xl)
    
    tmp <- split(table, group2)
    xbar_k <- sapply(tmp, function(x) colMeans(x[,resp]))
    diff_xk <- sweep(xbar_k, 1, xbar)
    sum_g2 <- g*n * diff_xk %*% t(diff_xk)
    
    mtx_intrc <- aggregate(table[,resp], list(group1,group2), FUN=mean)
    tmp_xl <- xbar_l[,rep(levels(group1), times = b)]
    tmp_xk <- xbar_k[,rep(levels(group2), each = g)]
    tmp_xbar <- matrix(xbar, nrow = p, ncol = b*g)
    mtx_intrc <- t(mtx_intrc[,resp]) - tmp_xl - tmp_xk + tmp_xbar
    sum_intrc <- n * mtx_intrc %*% t(mtx_intrc)
    
    diff_lkr <- sweep(t(table[,resp]), 1, xbar)
    sum_tot <- diff_lkr %*% t(diff_lkr)
    
    sum_res <- sum_tot - sum_g1 - sum_g2 - sum_intrc
    
    det_res <- det(sum_res)
    div_wilks <- c(det(sum_intrc+sum_res), det(sum_g1+sum_res),
                   det(sum_g2+sum_res))
    vec_wilks <- det_res / div_wilks
    
    # calculate chi-square value
    vec_x2 <- rep(0, times = 3)
    vec_x2[1] <- -(g*b*(n-1) - (p+1-(g-1)*(b-1))/2) * log(vec_wilks[1])
    vec_x2[2] <- -(g*b*(n-1) - (p+1-(g-1))/2) * log(vec_wilks[2])
    vec_x2[3] <- -(g*b*(n-1) - (p+1-(b-1))/2) * log(vec_wilks[3])
    
    # find p-value
    vec_df <- c((g-1)*(b-1)*p, (g-1)*p, (b-1)*p)
    vec_p.val <- mapply(
        function(x2,df) pchisq(x2, df, lower.tail = FALSE),
        vec_x2, vec_df
    )
    
    # conclusion
    vec_hypo <- sapply(
        vec_p.val, function(p.val) 
        ifelse(p.val <= 0.05, 'Reject Ho', 'Fail to reject')
    )
    
    # text
    vec_text <- c('Itrct', 'Fac1', 'Fac2')
    
    cat('sov','wilks','chi-sq','df','p-val','conclusion (alpha=0.05)',
        sep='\t')
    cat('\n')
    for (i in 1:3) {
        cat(vec_text[i], round(vec_wilks[i], 2),round(vec_x2[i], 2),
            vec_df[i], round(vec_p.val[i], 3), vec_hypo[i],sep='\t')
        cat('\n')
    }
    
    list(wilks = vec_wilks, x2 = vec_x2, df = vec_df, pval = vec_p.val,
         hyp = vec_hypo)
}

# box.m_result <- box_m(my_data, c('penerimaan', 'pengeluaran'), my_data$pend_krt)
# manova_result <- twoway_manova(
#     my_data,
#     c('penerimaan', 'pengeluaran'),
#     my_data$pend_krt, my_data$jumlah_art
# )