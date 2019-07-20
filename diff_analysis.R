args = commandArgs(T)
if(length(args) != 5){
	cat("R script for wilcox and t test and statistical result output\n\n")
	cat("Rscript diff_analysis.R mat.table group.list W prefix outdir\n W for wilcox.test() T for t.test()\n\n")
	cat("mat.table format: column is samples, row is vars\n\n")
	cat("group.list format: first column is samples, second is groups\n\n")
	quit()
}
mat <- args[1]
map <- args[2]
prefix <- args[4]
outdir <- args[5]
m <- args[3]

w.res.get <- function(group, mat, method = c('T', 'W'), adjust = TRUE){

  #mat type: column is vars and row is samples
  #group type: first column is sample, second is groups
  
  method = match.arg(method)

  unqgs <- unique(group[,2])
  
  idx <- seq_len(length(group[,1]))
  
  g1idx <- idx[group[,2] == unqgs[1]]
  
  g2idx <- idx[group[,2] == unqgs[2]]
  
  n.bm <- ncol(mat)
  
  w <- p <- diff.type <- FC <- diff.res <- M1 <- M2 <- vector(length = n.bm)

  for(i in seq_along(1:n.bm)){
    allv = as.numeric(mat[, i])
    if(method == 'W'){
	  w.res <- wilcox.test(allv[g1idx], allv[g2idx], exact = T)
	}else{
	  w.res <- t.test(allv[g1idx], allv[g2idx])
	}
    w[i] <- w.res$statistic
    p[i] <- w.res$p.value
    diff.type[i] <- paste0(unqgs[1],'-',unqgs[2])
    M1[i] <- mean(allv[g1idx])
    M2[i] <- mean(allv[g2idx])
    # cat(M1[i], "\t", M2[i], "\n")
    diff.res[i] <- mean(allv[g1idx]) - mean(allv[g2idx])
    FC[i] <- mean(allv[g1idx])/mean(allv[g2idx])
  }
  
  q <- p.adjust(p, method = 'BH')
  
  diff_R <- data.frame(var = colnames(mat), M1, M2, 
                       diff.type, diff.res, FC = FC, w, p, q)
 
  Stats <- ifelse(method == 'W', 'W-stats', 'T-stats') 
  colnames(diff_R) <- c('ID', paste0(unqgs[1], '_mean'),
                        paste0(unqgs[2], '_mean'),
                        'vs_group', 'Diff', 'FC', 
                        Stats, 'p_value', 'q_value')
  return(diff_R)
  
}


micro <- read.table(mat, sep = '\t', 
                    header = T, row.names = 1, stringsAsFactors = F)

group <- read.table(map, sep = '\t', header = F, stringsAsFactors = F)


adj.group <- group[match(colnames(micro), group$V1), ]


w.micro <- w.res.get(group = adj.group, mat = t(micro), method = m, adjust = TRUE)

if(!file.exists(outdir))
	dir.create(outdir)

out <- ifelse(m == 'W', paste0(outdir, '/', prefix, '.wilcox-test.result.xls'), 
			  paste0(outdir, '/', prefix, '.t-test.result.xls'))

write.table(w.micro, file = out, quote = F, row.names = F, sep = '\t')


