##eQTL analysis on TCGA-BRCA data for a given set of SNPs

library(data.table)

load("./data/tcga_rna_seq_tables.RData")
load("./data/tcga_rna_seq_ids.RData")
colnames(tcga.tumor.er.pos.rsem)[1] <- "gene"
ld.snp.only.bed.final <- fread("./data/ccv/brca_table_incl_ccv_final.txt")
gwas.snp.id.all <- unique(ld.snp.only.bed.final$GWAS_snp)
impute.results.dir <- "./TCGA/BRCA2/Genotype/processed/impute"
output.dir <- "./data/eQTL"

file.map <- fread("./TCGA/BRCA2/Genotype/FILE_SAMPLE_MAP_final.txt")
file.map.nonffpe <- file.map[!filename %like% "FFPE"]

tcga.tumor.er.pos.rsem[, "gene" := gsub("^(.*)[|].*$", "\\1", gene)]

gencode.v19.name <- fread("./data/misc/wgEncodeGencodeBasicV19_red.txt")
gene.info <- fread("./data/misc/Homo_sapiens.gene_info")
gencode.v22.attr <- fread("./data/misc/wgEncodeGencodeAttrsV22.txt")
gencode.v19.name[, "GeneID" := gene.info[match(name2, gene.info$Symbol), GeneID]]
gencode.v22.hg19 <- fread("./data/misc/wgEncodeGencodeBasicV22_hg19_liftover.bed")
colnames(gencode.v22.hg19) <- c("chrom", "txStart", "txEnd", "geneId")
tcga.gene.table.gaf <- fread("./data/misc/TCGA.hg19.gene.table.gaf")

brca.clinical <-  fread("./data/shared/TCGA/BRCA/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
brca.clinical.er.pos <- brca.clinical[er_status_by_ihc == "Positive"]
brca.clinical.er.pos.pr.pos <- brca.clinical.er.pos[pr_status_by_ihc == "Positive"]
brca.clinical.er.pos.pr.neg <- brca.clinical.er.pos[pr_status_by_ihc == "Negative"]
brca.clinical.er.pos.pr.pos.her2.pos <- brca.clinical.er.pos.pr.pos[her2_status_by_ihc == "Positive"]
brca.clinical.er.pos.pr.pos.her2.neg <- brca.clinical.er.pos.pr.pos[her2_status_by_ihc == "Negative"]

window.length <- 2e6
impute.server <- "michigan"
start <- 1


eqtl.results <- sapply(gwas.snp.id.all, function(gwas.snp.id) {
  setwd(impute.results.dir)
  gwas.snp.chr <- ld.snp.only.bed.final[LD_snp == gwas.snp.id & GWAS_snp == gwas.snp.id, chr]
  gwas.snp.pos <- ld.snp.only.bed.final[LD_snp == gwas.snp.id & GWAS_snp == gwas.snp.id, end]
  snp.list.for.eqtl <- ld.snp.only.bed.final[GWAS_snp == gwas.snp.id, LD_snp]
  snp.pos.list.for.eqtl <- ld.snp.only.bed.final[GWAS_snp == gwas.snp.id, end]
  window.genes <- tcga.gene.table.gaf[txStart >= gwas.snp.pos - window.length & txStart <= gwas.snp.pos + window.length & chrom == gwas.snp.chr, Gene]
  window.genes <- gsub("^(.*)[|].*$", "\\1", window.genes)
  gwas.imputed.gt.all <- fread(paste0("./TCGA/BRCA2/Genotype/processed/impute/", gwas.snp.chr, "/output_raw/", impute.server, "_server/", gwas.snp.chr, "_", gwas.snp.id, "_locus2.vcf"), autostart = 104)
  
  gwas.imputed.gt <- gwas.imputed.gt.all[POS %in% snp.pos.list.for.eqtl]

  start <- Sys.time()
  pval.list.temp <- sapply(snp.list.for.eqtl, function(snp.for.eqtl) {
    library(data.table)
    snp.chr <- ld.snp.only.bed.final[LD_snp == snp.for.eqtl, chr][1]
    snp.pos <- ld.snp.only.bed.final[LD_snp == snp.for.eqtl, end][1]
    if (snp.pos %in% gwas.imputed.gt$POS) {
      gt.threshold <- 0.9
      ld.snp.full <- gwas.imputed.gt[POS %in% snp.pos]
      ld.snp.gt.full <- ld.snp.full[, 10:ncol(ld.snp.full), with = F]
      ld.snp.gt <- as.data.table(lapply(ld.snp.gt.full, function(str) gsub("^(...)[:].+", "\\1", str)))
      ld.snp.dos <- as.data.table(lapply(ld.snp.gt.full, function(str) gsub("^.+[:](.+)", "\\1", str)))
      bool.conf <- lapply(strsplit(as.character(ld.snp.dos[1]), ","), max) >= gt.threshold
      
      snp.gt.full <- ld.snp.gt[1, bool.conf, with = F]
      file.name.full <- paste0(colnames(snp.gt.full), ".birdseed.data.txt")
      if (is.na(snp.gt.full[1, 1, with = F]) == FALSE) {
        gene.pval <- t(sapply(window.genes, function(gene.name) {
          library(data.table)
          gene.exp.rsem <- round(mean(as.numeric(tcga.tumor.er.pos.rsem[gene == gene.name, 2:ncol(tcga.tumor.er.pos.rsem), with = F])), digits = 2)
          gt.calls <- data.table(sapply(colnames(snp.gt.full), function(name) {
            if (snp.gt.full[, name, with = F] == "0|0")
              return(0)
            else if (snp.gt.full[, name, with = F] == "1|0" | snp.gt.full[, name, with = F] == "0|1")
              return(1)
            else
              return(2)
          }))
          rownames(gt.calls) <- file.name.full
          colnames(gt.calls) <- "Call"
          pval.no.seg <- NA
          pval <- NA
          if (!is.na(gene.exp.rsem)) {
            if (gene.exp.rsem >= 1) {
              eqtl.result <- analyze.eqtl.tcga(data = "BRCA", normalization = "RSEM", sample = "tumor", gene.name = gene.name, snp.id = snp.for.eqtl, patients.rm = NULL, gt.calls, isoform = FALSE, hg38 = FALSE)
              total.exp.all <- as.numeric(c(eqtl.result[[1]], eqtl.result[[2]], eqtl.result[[3]]))
              total.exp <- total.exp.all[total.exp.all != 0]
              snp.gt <- c(rep(0, length(eqtl.result[[1]])), rep(1, length(eqtl.result[[2]])), rep(2, length(eqtl.result[[3]])))[total.exp.all != 0]
              print(table(snp.gt))
              AA.patient.id <- eqtl.result$id.AA
              AB.patient.id <- eqtl.result$id.AB
              BB.patient.id <- eqtl.result$id.BB
              pat.id <- c(AA.patient.id, AB.patient.id, BB.patient.id)[total.exp.all != 0]
              if (length(pat.id) > 20 & length(AA.patient.id) >= 5 & length(AB.patient.id) >= 5 & length(BB.patient.id) >= 5) {
                gene.chr <- tcga.gene.table.gaf[Gene %like% paste0("^", gene.name, "[|]"), chrom]
                gene.tss <- tcga.gene.table.gaf[Gene %like% paste0("^", gene.name, "[|]"), txStart]
                gene.strand <- tcga.gene.table.gaf[Gene %like% paste0("^", gene.name, "[|]"), strand]
                txStart <- gene.tss
                txEnd <- tcga.gene.table.gaf[Gene %like% paste0("^", gene.name, "[|]"), txEnd]
                loh.seg.files <- sapply(pat.id, function(id) {
                  file <- paste0(gsub("^(.+)\\.birdseed.*$", "\\1", file.map.nonffpe[substr(`barcode(s)`, 1, 15) %in% paste0(id, "-01"), filename]), ".hg19.seg.txt")
                  if (length(file) > 1)
                    file <- file[1]
                  return(file)
                })
                setwd(copy.number.dir)
                cnv.seg.vals <- sapply(loh.seg.files, function(seg.name) {
                  if (file.exists(seg.name))
                    seg.file <- fread(seg.name)
                  else
                    return(list(NULL))
                  list(seg.file[Chromosome == substr(gene.chr, 4, 5)])
                })
                locus.seg.vals <- unlist(sapply(cnv.seg.vals, function(seg) {
                  if(!is.null(seg)) {
                    seg.mean.table <- seg[Start <= txStart & End >= txEnd]
                    if (nrow(seg.mean.table) == 0)
                      seg.mean.table <- seg[(Start <= txStart & End <= txEnd & End >= txStart) | (Start >= txStart & End <= txEnd) | (Start >= txStart & End >= txEnd & Start <= txEnd)]
                    #seg.mean <- seg[Start <= gene.tss & End >= gene.tss]$Segment_Mean
                    if (nrow(seg.mean.table) == 0)
                      return(0)
                    if (nrow(seg.mean.table) == 1)
                      return(seg.mean.table$Segment_Mean)
                    if (nrow(seg.mean.table) > 1) {
                      weight <- 0
                      total.length <- 0
                      sapply(1:nrow(seg.mean.table), function(row) {
                        if (row == 1) {
                          len1 <- seg.mean.table[row, End] - txStart
                          weight <<- weight + len1 * seg.mean.table[row, Segment_Mean]
                          total.length <<- total.length + len1
                        }
                        else if (row == nrow(seg.mean.table)) {
                          len2 <- txEnd - seg.mean.table[row, Start]
                          weight <<- weight + len2 * seg.mean.table[row, Segment_Mean]
                          total.length <<- total.length + len2
                        }
                        else {
                          len3 <- seg.mean.table[row, End] - seg.mean.table[row, Start]
                          weight <<- weight + len3 * seg.mean.table[row, Segment_Mean]
                          total.length <<- total.length + len3
                        }
                      })
                      seg.mean <- weight/total.length
                    }
                    return(seg.mean)
                  }
                  else {
                    return(0)
                  }
                }))
                lin.reg.res <- summary(lm(log2(total.exp+1) ~ snp.gt + locus.seg.vals))
                norm.exp <- resid((lm(log2(total.exp+1) ~ locus.seg.vals)))
                if (nrow(lin.reg.res$coefficients) == 2) {
                  pval <- 1
                } else {
                  pval <- signif(lin.reg.res$coefficients[2, 4], digits = 3)
                }
                if (nrow(lin.reg.res$coefficients) != 1)
                  if (is.na(lin.reg.res$coefficients[2, 4]))
                    pval <- 1
                setwd(output.dir)
                if (pval <= 0.05) {
                  filename1 = paste0("./expression/RSEM/gtcn_mich_", snp.for.eqtl, "_eqtl__", gene.name, "__expr_rsem.txt")
                  filename2 = paste0("./pval/RSEM/gtcn_mich_", snp.for.eqtl, "_eqtl_rsem.txt")
                  result.write1 <- cbind(gene = gene.name, ref_ref = paste(round(norm.exp[1:length(eqtl.result[[1]][eqtl.result[[1]] != 0])], 2), collapse = ","), ref_alt = paste(round(norm.exp[(length(eqtl.result[[1]][eqtl.result[[1]] != 0])+1):(length(eqtl.result[[1]][eqtl.result[[1]] != 0])+length(eqtl.result[[2]][eqtl.result[[2]] != 0]))], 2), collapse = ","), alt_alt = paste(round(norm.exp[(length(eqtl.result[[1]][eqtl.result[[1]] != 0])+length(eqtl.result[[2]][eqtl.result[[2]] != 0])+1):length(pat.id)], 2), collapse = ","))
                  result.write2 <- cbind(chrom = gene.chr, gene = gene.name, coeffGT = round(lin.reg.res$coefficients[2, 1], digits = 2), pvalueGT = pval, meanExp = round(mean(total.exp), 2), numPatients = paste(c(length(AA.patient.id), length(AB.patient.id), length(BB.patient.id)), collapse = ", "))
                  write.table(result.write1, filename1, row.names = F, col.names = T, quote = F, sep = "\t")
                  if (file.exists(filename2))
                    write.table(result.write2, filename2, append = T, row.names = F, col.names = F, quote = F, sep = "\t")
                  else
                    write.table(result.write2, filename2, row.names = F, col.names = T, quote = F, sep = "\t")
                }
              }
            }
          }
        }, simplify = FALSE))
      }
      setwd(output.dir)
      filename2 <- paste0("./pval/RSEM/gtcn_mich_", snp.for.eqtl, "_eqtl_rsem.txt")
      if (file.exists(filename2)) {
        data <- fread(filename2)
        data <- data[order(pvalueGT), ]
        write.table(data, filename2, row.names = F, col.names = T, quote = F, sep = "\t")
      }
    }}, simplify = FALSE)
  end <- Sys.time()
  end - start
}, simplify = FALSE)
end <- Sys.time()
end - start
