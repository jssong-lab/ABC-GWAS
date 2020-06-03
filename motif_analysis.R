library(data.table)
library(motifbreakR)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqLogo)

global.start <- Sys.time()
package.time <- Sys.time() - package.start
no_cores <- 22
cat("Loading motifs...\n")
data(motifbreakR_motif)
motifdb.hsapiens.set <- query(MotifDb, 'hsapiens')
values(motifbreakR_motif)[values(motifbreakR_motif)$dataSource %in% "FactorBook", "geneSymbol"] <- values(motifbreakR_motif)[values(motifbreakR_motif)$dataSource %in% "FactorBook", ]$providerId
motifdb.hsapiens <- c(motifdb.hsapiens.set, motifbreakR_motif)
motifdb.hsapiens <- motifdb.hsapiens[!is.na(values(motifdb.hsapiens)$geneSymbol) & !values(motifdb.hsapiens)$dataSource %in% c("HOCOMOCO", "JASPAR_2014", "jaspar2016")]

gwas.snp.list <- unique(fread("./data/brca_table_incl_ccv_final.txt")$LD_snp)
snp.list.new <- gwas.snp.list
b <- c(0.3, 0.2, 0.2, 0.3) #A, C, G, T
cat("Loading variants...\n")
variants <- snps.from.rsid(rsid = snp.list.new, dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37, search.genome = BSgenome.Hsapiens.UCSC.hg19)

cat("Analyzing motifs...\n")
motifbreakr.results.pval <- motifbreakR(snpList = variants, pwmList = motifdb.hsapiens, filterp = TRUE, threshold = 1e-3, method = "log", bkg = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3), verbose = TRUE, show.neutral = FALSE)

motifbreakr.results <- motifbreakr.results.pval
#save(list = c("variants", "motifbreakr.results"), file = paste0("./results/motifbreak_results_abc-gwas_1e-3.RData"))
#load(file = "./results/motifbreak_results_abc-gwas_1e-3.RData")
#motifbreakr.results <- motifbreakr.results[400001:length(motifbreakr.results), ]

#Compute p-values
cat("Calculating permutation test p-values...\n")
cl <- makeCluster(no_cores)
clusterExport(cl, c("b", "motifdb.hsapiens", "query", "reverseComplement", "DNAString"))
pval.start <- Sys.time()
get.score <- function(snp.id.list, seq.list, data.source.list, motif.name.list, motif.id.list, strand.vec, motif.pos.list, ref.allele.list, alt.allele.list, score.ref.list, score.alt.list) {
seq.score <- clusterMap(cl, function(snp.id, seq, data.source, motif.name, motif.id, strand, query.pos, ref.allele, alt.allele, score.ref, score.alt) {
  if (motif.name == "JUN (var.2)")
    motif.name <- "MA0489.1"
  if (motif.name == "JUND (var.2)")
    motif.name <- "MA0492.1"
  if (data.source == "HOMER")
    motif.name <- motif.id
  
  pwm <- query(query(query(motifdb.hsapiens, data.source), paste0("^", motif.name, "$")), paste0("^", motif.id, "$"))[[1]]

  pwm2 <- pwm + 0.05
  for (j in 1:ncol(pwm2)) {
    pwm2[1:4, j] <- pwm2[1:4, j]/sum(pwm2[1:4, j])
  }
  pwm <- pwm2
  rownames(pwm) <- c("A", "C", "G", "T")
  seq.match <- seq
  if (score.ref > score.alt) {
    query.nuc <- ref.allele
    mut.nuc <- alt.allele
  } else {
    splitseq <- strsplit(seq.match, split = "")[[1]]
    splitseq[query.pos] <- alt.allele
    seq.match <- paste0(splitseq, collapse = "")
    query.nuc <- alt.allele
    mut.nuc <- ref.allele
  }
  if (strand == "-") {
    seq.match <- as.character(reverseComplement(DNAString(toupper(seq))))
    if (score.ref > score.alt) {
      query.nuc <- as.character(reverseComplement(DNAString(ref.allele)))
      mut.nuc <- as.character(reverseComplement(DNAString(alt.allele)))
    } else {
      query.nuc <- as.character(reverseComplement(DNAString(alt.allele)))
      mut.nuc <- as.character(reverseComplement(DNAString(ref.allele)))
      splitseq <- strsplit(seq.match, split = "")[[1]]
      splitseq[query.pos] <- query.nuc
      seq.match <- paste0(splitseq, collapse = "")
    }

  }
  ref.score <- log(pwm[query.nuc, query.pos]/pwm[mut.nuc, query.pos])

  ic <- apply(pwm, 2, function(col) {
    ic.temp <- sum(col*log(col/b))
    return(ic.temp)
  })
  inverseIC <- 1/ic
  inverseIC <- inverseIC/sum(inverseIC)

  N <- 5000
  Nrep <- 1
  pval <- 1:Nrep
  out <- sapply(1:Nrep, function(m) {
    position.sample <- sample(1:length(rand.seq), N, prob = inverseIC, replace = T)
    vals <- sapply(1:N, function(j) {
      rand.seq <- strsplit(seq.match, "")[[1]]
      old <- rand.seq[position.sample[j]]
      rand.seq[position.sample[j]] <- sample(c("A", "C", "G", "T"), 1, prob = pwm[, position.sample[j]])
      result <- log(pwm[old, position.sample[j]]/pwm[rand.seq[position.sample[j]], position.sample[j]])
      return(result)
    })
    pval[m] <<- sum(vals >= ref.score)/N
  })
  Pvalue <- sum(pval)/Nrep
  
  if (file.exists("./results/pvals_all.txt"))
    write.table(cbind(rsID = snp.id, dataSource = data.source, motifName = motif.name, motifID = motif.id, seqMatch = seq, permPvalue = Pvalue), "./results/pvals_all.txt", row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  else
    write.table(cbind(rsID = snp.id, dataSource = data.source, motifName = motif.name, motifID = motif.id, seqMatch = seq, permPvalue = Pvalue), "./results/pvals_all.txt", row.names = F, col.names = T, quote = F, sep = "\t")
  Pvalue
}, snp.id.list, seq.list, data.source.list, motif.name.list, motif.id.list, strand.vec, motif.pos.list, ref.allele.list, alt.allele.list, score.ref.list, score.alt.list, SIMPLIFY = TRUE)
}
start <- Sys.time()
pvals <- get.score(names(motifbreakr.results), sub("\\s+$", "", motifbreakr.results$seqMatch), motifbreakr.results$dataSource, motifbreakr.results$providerName, motifbreakr.results$providerId, as.vector(strand(motifbreakr.results)), motifbreakr.results$motifPos, as.character(motifbreakr.results$REF), as.character(motifbreakr.results$ALT), motifbreakr.results$scoreRef, motifbreakr.results$scoreAlt)
save.image(file = paste0("./results/motifbreak_results_pvals_abc-gwas_1e-3.RData"))
stopCluster(cl)