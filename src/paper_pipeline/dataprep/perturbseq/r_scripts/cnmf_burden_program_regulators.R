args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop(
    paste(
      "usage: cnmf_burden_program_regulators.R <gep.txt> <gene_map.txt> <regulation_dir>",
      "<lof.tsv> <shet.txt> <k> <out_regulators.txt> <out_programs.txt>",
      "<random_iterations> <top_n>"
    )
  )
}

options("scipen" = 10)

gep_path <- args[1]
gene_map_path <- args[2]
regulation_dir <- args[3]
lof_path <- args[4]
shet_path <- args[5]
K <- as.numeric(args[6])
out_reg_path <- args[7]
out_program_path <- args[8]
random_iterations <- as.numeric(args[9])
top_n <- as.numeric(args[10])

GEP <- read.table(gep_path, header = TRUE, stringsAsFactor = FALSE)
GEP <- t(GEP)
colnames(GEP) <- paste0("P", c(1:K))

corresp <- read.table(gene_map_path, header = FALSE, stringsAsFactor = FALSE)
corresp <- corresp[!duplicated(corresp[, 1]), ]

# gene_map_path is expected to be "<ENSG> <gene_symbol>".
# Keep explicit maps in both directions to avoid mixing ENSG-indexed
# lookups with symbol-indexed lookups later in the script.
ensg_to_gene <- corresp[, 2]
names(ensg_to_gene) <- corresp[, 1]

symbol_map <- corresp[!duplicated(corresp[, 2]), ]
symbol_to_ensg <- symbol_map[, 1]
names(symbol_to_ensg) <- symbol_map[, 2]

GEP <- data.frame(GENE = unname(ensg_to_gene[row.names(GEP)]), GEP)

GEP_reg_beta <- data.frame()
for (i in 1:K) {
  tmp <- read.table(
    file.path(regulation_dir, paste0("K", K, "_program", i, "_perturb_effects.txt")),
    header = TRUE
  )
  tmp1 <- tmp[, c(1, 2)]
  colnames(tmp1) <- c("GENE", paste0("P", i))

  if (i == 1) {
    GEP_reg_beta <- tmp1
  } else {
    GEP_reg_beta <- merge(GEP_reg_beta, tmp1, by = "GENE")
  }
}

LOF <- read.table(lof_path, sep = "\t", quote = "", header = TRUE, stringsAsFactor = FALSE)
LOF$post_mean[is.element(LOF$post_mean, "Inf")] <- max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")] <- min(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF <- data.frame(LOF, gene = unname(ensg_to_gene[as.character(LOF$ensg)]))

GEP2 <- data.frame(ensg = row.names(GEP), GEP)
df <- merge(GEP2, LOF, by = "ensg")

shet <- read.table(shet_path, header = TRUE, stringsAsFactor = FALSE)
shet <- shet[is.element(shet$ensg, LOF$ensg), ]

summary_reg <- data.frame()
summary_pro <- data.frame()

for (Program in c(1:K)) {
  tmp <- GEP_reg_beta[, c("GENE", paste0("P", Program))]
  colnames(tmp) <- c("GENE", "perturb_beta")
  tmp <- data.frame(tmp, ensg = unname(symbol_to_ensg[as.character(tmp$GENE)]))
  df <- merge(tmp, LOF, by = "ensg")

  P_pearson_all <- cor.test(df$perturb_beta, df$post_mean, method = "pearson")$p.value
  R_pearson_all <- cor.test(df$perturb_beta, df$post_mean, method = "pearson")$estimate
  R_pearson_all_CIlower <- cor.test(df$perturb_beta, df$post_mean, method = "pearson")$conf.int[1]
  R_pearson_all_CIupper <- cor.test(df$perturb_beta, df$post_mean, method = "pearson")$conf.int[2]

  df <- merge(df, shet, by = "ensg")
  df$post_mean <- scale(df$post_mean)
  df$perturb_beta <- scale(df$perturb_beta)

  fit <- lm(post_mean ~ perturb_beta + shet, data = df)
  P_withShet <- summary(fit)$coefficients[2, 4]
  beta_withShet <- summary(fit)$coefficients[2, 1]
  betaSE_withShet <- summary(fit)$coefficients[2, 2]

  hoge <- data.frame(
    FILE = basename(lof_path),
    Program = Program,
    P_pearson_all = P_pearson_all,
    R_pearson_all = R_pearson_all,
    R_pearson_all_CIlower = R_pearson_all_CIlower,
    R_pearson_all_CIupper = R_pearson_all_CIupper,
    P_withShet = P_withShet,
    beta_withShet = beta_withShet,
    betaSE_withShet = betaSE_withShet
  )
  summary_reg <- rbind(summary_reg, hoge)

  tmp <- GEP[, c("GENE", paste0("P", Program))]
  tmp <- tmp[order(tmp[, 2], decreasing = TRUE), ]
  colnames(tmp)[2] <- "GEPscore"
  tmp <- data.frame(tmp, ensg = row.names(tmp))
  df <- merge(tmp, LOF, by = "ensg")
  df2 <- df[order(df$GEPscore, decreasing = TRUE), ][1:top_n, ]

  mean_gamma_top <- mean(df2[, "post_mean"])

  shet_tmp <- shet[is.element(shet$ensg, df2$ensg), ]
  shet_tmp_all <- shet[is.element(shet$ensg, df$ensg), ]
  A <- table(shet_tmp$shet_BIN)

  random <- c()
  for (I in 1:random_iterations) {
    set.seed(I)
    genes <- c()
    for (p in 1:length(A)) {
      genes <- c(
        genes,
        sample(shet_tmp_all$ensg[is.element(shet_tmp_all$shet_BIN, names(A)[p])], A[p])
      )
    }
    random <- c(random, mean(df[is.element(df$ensg, genes), "post_mean"], na.rm = TRUE))
  }

  random_mean <- mean(random)
  P1 <- (rank(c(mean_gamma_top, random))[1] / length(random)) * 2
  P2 <- ((length(random) + 2 - (rank(c(mean_gamma_top, random))[1])) / length(random)) * 2
  program_P <- min(P1, P2)

  hoge <- data.frame(
    FILE = basename(lof_path),
    Program = Program,
    MEANgamma_top100 = mean_gamma_top,
    shet_adjusted_random_mean = random_mean,
    MEANgamma_top100_shet_adjusted_P = program_P
  )
  summary_pro <- rbind(summary_pro, hoge)
}

dir.create(dirname(out_reg_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_program_path), showWarnings = FALSE, recursive = TRUE)
write.table(summary_reg, out_reg_path, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(summary_pro, out_program_path, row.names = FALSE, sep = "\t", quote = FALSE)
