#compute correlation between target genes and transcription factors
library(data.table)
library(Biobase)

setwd("./data/")
load("./data/tcga_rna_seq_tables.RData")
colnames(tcga.tumor.er.pos.rsem)[1] <- "gene"
tcga.tumor.er.pos.rsem[, "gene" := gsub("^(.*)[|].*$", "\\1", gene)]
snps.list <- fread("./data/ccv/brca_table_incl_ccv_final.txt")
colnames(snps.list) <- c("chr", "start", "end", "LD_snp", "GWAS_snp")

gwas.snps.list <- unique(snps.list$GWAS_snp)
motif.tables <- list.files("./motifs/new_snp_tables/")
motif.analysis.snps <- gsub("^(.*?)[_].*$", "\\1", motif.tables)

no_cores <- 14

count <- 1
temp.out <- sapply(gwas.snps.list, function(query.snp) {
  print(count)
  ld.snp.all <- snps.list[GWAS_snp == query.snp, LD_snp]
  ld.snp.all <- ld.snp.all[ld.snp.all %in% motif.analysis.snps]
  out <- mclapply(ld.snp.all, function(ld.snp) {
    library(data.table)
    motif.result <- fread(paste0("./motifs/new_snp_tables/", ld.snp, "_motif_table.txt"))
    if (nrow(fread(paste0("./eQTL/pval/RSEM/gtcn_mich_", query.snp, "_eqtl_rsem.txt"))) > 0) {
      eqtl.genes <- fread(paste0("./eQTL/pval/RSEM/gtcn_mich_", query.snp, "_eqtl_rsem.txt"))$gene
      gene.set.1 <- eqtl.genes
      gene.set.2 <- unique(motif.result$geneSymbol)

      gene.exp.2 <- sapply(gene.set.2, function(gene.name) {
        log2(mean(as.numeric(tcga.tumor.er.pos.rsem[gene == gene.name, 2:ncol(tcga.tumor.er.pos.rsem), with = F][1]))+1)
      })
      gene.exp.raw.2 <- gene.exp.2
      gene.set.2 <- unique(gene.set.2[!is.na(gene.exp.raw.2) & gene.exp.raw.2 >= 1])
  
      compute.corr <- function(set.1, set.2, corr.method, rna.seq.data, gt.id) {
        rna.seq.cols.select <- c("gene", colnames(rna.seq.data)[2:ncol(rna.seq.data)][substr(colnames(rna.seq.data)[2:ncol(rna.seq.data)], 1, 12) %in% gt.id])
        tcga.rna.seq.temp <- rna.seq.data[, rna.seq.cols.select, with = F]
        corr.coeff <- sapply(set.1, function(gene1) {
          cor.temp <- t(sapply(set.2, function(gene2) {
            gene1.exp <- log2(as.numeric(tcga.rna.seq.temp[gene == gene1, 2:ncol(tcga.rna.seq.temp), with = F][1])+1)
            gene2.exp <- log2(as.numeric(tcga.rna.seq.temp[gene == gene2, 2:ncol(tcga.rna.seq.temp), with = F][1])+1)
            if (sum(is.na(gene1.exp)) == 0 & sum(is.na(gene2.exp)) == 0) {
              gene.cor <- cor.test(gene1.exp, gene2.exp, method = corr.method)
              return(cbind(cor=gene.cor$estimate, exp=median(gene2.exp), pval=gene.cor$p.value))
            } else {
              return(list(NA, NA, NA))
            }
          }, simplify = T))
          colnames(cor.temp) <- c(gene1, "exp", "pval")
          return(cor.temp)
        }, simplify = F)
        return(corr.coeff[[1]])
      }
      
      sapply(gene.set.1, function(target.gene) {
        brca.temp <- tcga.tumor.er.pos.rsem[gene %in% c(target.gene, gene.set.2), ]
        AA.patient.id <- fread(paste0("./data/eQTL/patient_id/patient_id_for_website_corr_", ld.snp, "_AA.txt"), header = F)$V1
        AA.patient.id <- AA.patient.id[which(as.numeric(tcga.tumor.er.pos.rsem[gene == target.gene, match(AA.patient.id, substr(colnames(tcga.tumor.er.pos.rsem), 1, 12)), with = F][1]) > 0)]
        AB.patient.id <- fread(paste0("./data/eQTL/patient_id/patient_id_for_website_corr_", ld.snp, "_AB.txt"), header = F)$V1
        AB.patient.id <- AB.patient.id[which(as.numeric(tcga.tumor.er.pos.rsem[gene == target.gene, match(AB.patient.id, substr(colnames(tcga.tumor.er.pos.rsem), 1, 12)), with = F][1]) > 0)]
        BB.patient.id <- fread(paste0("./data/eQTL/patient_id/patient_id_for_website_corr_", ld.snp, "_BB.txt"), header = F)$V1
        BB.patient.id <- BB.patient.id[which(as.numeric(tcga.tumor.er.pos.rsem[gene == target.gene, match(BB.patient.id, substr(colnames(tcga.tumor.er.pos.rsem), 1, 12)), with = F][1]) > 0)]

        if (length(AA.patient.id) >= 5 & length(AB.patient.id) >= 5 & length(BB.patient.id) >= 5) {
          corr.coeff.AA <- compute.corr(target.gene, gene.set.2, "pear", brca.temp, AA.patient.id)
          corr.coeff.AB <- compute.corr(target.gene, gene.set.2, "pear", brca.temp, AB.patient.id)
          corr.coeff.BB <- compute.corr(target.gene, gene.set.2, "pear", brca.temp, BB.patient.id)

          out <- sapply(1:length(gene.set.2), function(i) {
            AA.corr <- round(corr.coeff.AA[gene.set.2[i], 1], 2)
            AB.corr <- round(corr.coeff.AB[gene.set.2[i], 1], 2)
            BB.corr <- round(corr.coeff.BB[gene.set.2[i], 1], 2)
            all.corr <- c(AA.corr, AB.corr, BB.corr)
            out <- sapply(1:3, function(j) {
              if (!file.exists(paste0("./correlation/RSEM/", query.snp, "_", ld.snp, "_", target.gene, "_correlation.txt")))
                write.table(cbind(tf_id = i, group = j, value = all.corr[j], tf = gene.set.2[i]), paste0("./correlation/RSEM/", query.snp, "_", ld.snp, "_", target.gene, "_correlation.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
              else
                write.table(cbind(tf_id = i, group = j, value = all.corr[j], tf = gene.set.2[i]), paste0("./correlation/RSEM/", query.snp, "_", ld.snp, "_", target.gene, "_correlation.txt"), append = T, col.names = F, row.names = F, quote = F, sep = "\t")
            })
          })

        }
      })
    }
  }, mc.cores = no_cores)
  count <<- count + 1
})