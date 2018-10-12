epl_dat <- read.csv('/storage/Code/data-r/apg/epldata_final.csv')
epl_dat$fpl_sel <- as.numeric(gsub('%', '', epl_dat$fpl_sel))/100
process_col <- c(
    'market_value', 'page_views', 'fpl_value', 'fpl_sel'
)
epl_process <- epl_dat[,process_col]
epl_process <- scale(epl_process)
clusters <- kmeans(epl_process, 5, iter.max = 10000, nstart = 1000)
epl_dat$cluster <- as.factor(clusters$cluster)
epl_res <- epl_dat[order(epl_dat$cluster),]
