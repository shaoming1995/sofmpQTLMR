#' @title PQTL和在线结局MR分析
#' @param savefile 结果输出的文件夹名称
#' @param PATH 切分好的文件夹
#' @param GWASID 结局的GWASID号
#' @param outname 结局的名称
#' @param kb 聚类的距离
#' @param r2 聚类的r2
#' @export
PQTL_IEU<-function (savefile, PATH, GWASID, outname, kb, r2) {
  library(TwoSampleMR)
  dir.create(savefile)
  filename <- data.frame(dir(PATH))
  A_temp <- c()
  B_temp <- c()
  C_temp <- c()
  D_temp <- c()
  E_temp <- c()
  for (i in filename[, 1]) {
    ipath <- paste0(PATH, "/", i)
    exp_temp <- read.csv(ipath, header = T)
    test2 <- (try(exp_temp1 <- clump_data(exp_temp, clump_kb = kb,
                                          clump_r2 = r2)))
    if (class(test2) != "try-error") {
      test3<-(try(OUT <- extract_outcome_data(snps = exp_temp1$SNP,
                                  outcomes = GWASID, proxies = T, maf_threshold = 0.01,
                                  access_token = NULL)))
      if (class(test3)!="NULL") {
        OUT$id.outcome <- outname
        OUT$outcome <- outname
        OUT <- OUT[!duplicated(OUT$SNP), ]
        exp_temp_out <- merge(exp_temp1, OUT, by = "SNP",
                              all = F)
        exp <- exp_temp_out[, c("SNP", "effect_allele.exposure",
                                "other_allele.exposure", "beta.exposure", "se.exposure",
                                "pval.exposure", "id.exposure", "exposure",
                                "eaf.exposure")]
        out <- exp_temp_out[, c("SNP", "effect_allele.outcome",
                                "other_allele.outcome", "eaf.outcome", "beta.outcome",
                                "se.outcome", "pval.outcome", "id.outcome",
                                "outcome", "eaf.outcome")]
        dat <- harmonise_data(exposure_dat = exp, outcome_dat = out,
                              action = 2)
        test4<-(try(data_h <- dat %>% subset(dat$mr_keep == TRUE)))
        if (dim(data_h)[[1]]!=0){
          res <- mr(data_h)

          data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure) *
            (data_h$beta.exposure/data_h$se.exposure)
          data_h_TableS1 <- data_h[, c("exposure", "SNP",
                                       "effect_allele.exposure", "other_allele.exposure",
                                       "beta.exposure", "se.exposure", "Fvalue", "pval.exposure",
                                       "beta.outcome", "se.outcome", "pval.outcome")]
          res$cluster <- 1
          mr_OR <- generate_odds_ratios(res)
          mr_OR$or <- round(mr_OR$or, 3)
          mr_OR$or_lci95 <- round(mr_OR$or_lci95, 3)
          mr_OR$or_uci95 <- round(mr_OR$or_uci95, 3)
          mr_OR$OR_CI <- paste0(mr_OR$or, "(", mr_OR$or_lci95,
                                "-", mr_OR$or_uci95, ")")
          het <- mr_heterogeneity(dat)
          ple <- mr_pleiotropy_test(dat)
          A_temp <- rbind(mr_OR, A_temp)
          B_temp <- rbind(het, B_temp)
          C_temp <- rbind(ple, C_temp)
          print(paste0("当前运行到", i, "文件"))
          Aname <- paste0(savefile, "/", "PQTL与",
                          outname, "的MR结果.csv")
          Bname <- paste0(savefile, "/", "PQTL与",
                          outname, "的异质性结果.csv")
          Cname <- paste0(savefile, "/", "PQTL与",
                          outname, "的多效性结果.csv")
          Dname <- paste0(savefile, "/", "PQTL与",
                          outname, "的SNPs情况.csv")
          write.csv(A_temp, Aname, row.names = F)
          write.csv(B_temp, Bname, row.names = F)
          write.csv(C_temp, Cname, row.names = F)
          write.csv(D_temp, Dname, row.names = F)
        }else{E_temp<-rbind(i,E_temp)
        write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}}
      else{E_temp<-rbind(i,E_temp)
      write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}}
    else {
      cat(i, "由于网络502问题未完成clump")
      test5<-(try(OUT <- extract_outcome_data(snps = exp_temp$SNP,
                                  outcomes = GWASID, proxies = T, maf_threshold = 0.01,
                                  access_token = NULL)))
      if (class(test5)!="NULL") {
        OUT$id.outcome <- outname
        OUT$outcome <- outname
        OUT <- OUT[!duplicated(OUT$SNP), ]
        exp_temp_out <- merge(exp_temp, OUT, by = "SNP",
                              all = F)
        exp_temp_out$eaf.exposure <- exp_temp_out$eaf.outcome
        exp <- exp_temp_out[, c("SNP", "effect_allele.exposure",
                                "other_allele.exposure", "beta.exposure", "se.exposure",
                                "pval.exposure", "id.exposure", "exposure",
                                "eaf.exposure")]
        out <- exp_temp_out[, c("SNP", "effect_allele.outcome",
                                "other_allele.outcome", "eaf.outcome", "beta.outcome",
                                "se.outcome", "pval.outcome", "id.outcome",
                                "outcome", "eaf.outcome")]
        dat <- harmonise_data(exposure_dat = exp, outcome_dat = out,
                              action = 2)
        test6<-(try(data_h <- dat %>% subset(dat$mr_keep == TRUE)))
        if (dim(data_h)[[1]]!=0){
          res <- mr(dat)

          data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure) *
            (data_h$beta.exposure/data_h$se.exposure)
          data_h_TableS1 <- data_h[, c("exposure", "SNP",
                                       "effect_allele.exposure", "other_allele.exposure",
                                       "beta.exposure", "se.exposure", "Fvalue", "pval.exposure",
                                       "beta.outcome", "se.outcome", "pval.outcome")]
          res$cluster <- 0
          mr_OR <- generate_odds_ratios(res)
          mr_OR$or <- round(mr_OR$or, 3)
          mr_OR$or_lci95 <- round(mr_OR$or_lci95, 3)
          mr_OR$or_uci95 <- round(mr_OR$or_uci95, 3)
          mr_OR$OR_CI <- paste0(mr_OR$or, "(", mr_OR$or_lci95,
                                "-", mr_OR$or_uci95, ")")
          het <- mr_heterogeneity(dat)
          ple <- mr_pleiotropy_test(dat)
          A_temp <- rbind(mr_OR, A_temp)
          B_temp <- rbind(het, B_temp)
          C_temp <- rbind(ple, C_temp)
          print(paste0("当前运行到", i, "文件"))
          Aname <- paste0(savefile, "/", "PQTL与",
                          outname, "的MR结果.csv")
          Bname <- paste0(savefile, "/", "PQTL与",
                          outname, "的异质性结果.csv")
          Cname <- paste0(savefile, "/", "PQTL与",
                          outname, "的多效性结果.csv")
          Dname <- paste0(savefile, "/", "PQTL与",
                          outname, "的SNPs情况.csv")
          write.csv(A_temp, Aname, row.names = F)
          write.csv(B_temp, Bname, row.names = F)
          write.csv(C_temp, Cname, row.names = F)
          write.csv(D_temp, Dname, row.names = F)
        } else{E_temp<-rbind(i,E_temp)
        write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}}
      else{E_temp<-rbind(i,E_temp)
      write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}
    }

  }
  cat("当前分析已完成")
}

