#' @title PQTL和本地结局MR分析
#' @param name 学号
#' @param key 密码
#' @param savefile 结果输出的文件夹名称
#' @param PATH 切分好的文件夹
#' @param GWASsummay 结局的GWASsummary数据
#' @param outname 结局的名称
#' @param kb 聚类的距离
#' @param r2 聚类的r2
#' @export
PQTL_local<-function (name, key, savefile, PATH, GWASsummay, outname, kb, r2){
  filename <- data.frame(dir(PATH))
  A_temp <- c()
  B_temp <- c()
  C_temp <- c()
  D_temp <- c()
  E_temp <- c()
  for (i in filename[, 1]) {
    ipath <- paste0(PATH, "/", i)
    exp_temp <- read.csv(ipath, header = T)
    test2 <- (try(exp_temp <- clump_data(exp_temp, clump_kb = kb,
                                         clump_r2 = r2)))

    if (class(test2) == "try-error") {
      cat(i, "由于网络502问题未完成clump")
      test3<-(try(total <- merge(GWASsummay, exp_temp, by = "SNP")))
      if(class(test3)!="NULL"){
        exp3 <- total[, c("SNP", "effect_allele.exposure",
                          "other_allele.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "id.exposure", "exposure",
                          "eaf.exposure")]
        out3 <- total[, c("SNP", "effect_allele.outcome",
                          "other_allele.outcome", "eaf.outcome", "beta.outcome",
                          "se.outcome", "pval.outcome", "id.outcome",
                          "outcome", "eaf.outcome")]
        dat <- harmonise_data(exposure_dat = exp3, outcome_dat = out3,
                              action = 2)
        test4<-(try(data_h <- dat %>% subset(dat$mr_keep == TRUE)))
        if (class(test4)!="NULL"){
          data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure) *
            (data_h$beta.exposure/data_h$se.exposure)
          data_h_TableS1 <- data_h[, c("exposure", "SNP",
                                       "effect_allele.exposure", "other_allele.exposure",
                                       "beta.exposure", "se.exposure", "Fvalue", "pval.exposure",
                                       "beta.outcome", "se.outcome", "pval.outcome")]
          data_h_TableS1$cluster <- 0
          res <- mr(data_h)
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
          D_temp <- rbind(data_h_TableS1, D_temp)
          print(paste0("当前运行到", i, "文件"))
          Aname <- paste0(savefile, "/", "pQTL与",
                          outname, "的MR结果.csv")
          Bname <- paste0(savefile, "/", "pQTL群与",
                          outname, "的异质性结果.csv")
          Cname <- paste0(savefile, "/", "pQTL与",
                          outname, "的多效性结果.csv")
          Dname <- paste0(savefile, "/", "pQTL群与",
                          outname, "的SNPs情况.csv")
          write.csv(A_temp, Aname, row.names = F)
          write.csv(B_temp, Bname, row.names = F)
          write.csv(C_temp, Cname, row.names = F)
          write.csv(D_temp, Dname, row.names = F)}else{E_temp<-rbind(i,E_temp)
          write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}}
      else{E_temp<-rbind(i,E_temp)
      write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}}
    else {
      test5<-(try(total <- merge(GWASsummay, exp_temp, by = "SNP")))
      if(class(test5)!="NULL"){
        exp3 <- total[, c("SNP", "effect_allele.exposure",
                          "other_allele.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "id.exposure", "exposure",
                          "eaf.exposure")]
        out3 <- total[, c("SNP", "effect_allele.outcome",
                          "other_allele.outcome", "eaf.outcome", "beta.outcome",
                          "se.outcome", "pval.outcome", "id.outcome",
                          "outcome", "eaf.outcome")]
        dat <- harmonise_data(exposure_dat = exp3, outcome_dat = out3,
                              action = 2)
        test6<-(try(data_h <- dat %>% subset(dat$mr_keep == TRUE)))
        if (class(test6)!="NULL"){
          data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure) *
            (data_h$beta.exposure/data_h$se.exposure)
          data_h_TableS1 <- data_h[, c("exposure", "SNP",
                                       "effect_allele.exposure", "other_allele.exposure",
                                       "beta.exposure", "se.exposure", "Fvalue", "pval.exposure",
                                       "beta.outcome", "se.outcome", "pval.outcome")]
          data_h_TableS1$cluster <- 1
          res <- mr(data_h)
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
          D_temp <- rbind(data_h_TableS1, D_temp)
          print(paste0("当前运行到", i, "文件"))
          Aname <- paste0(savefile, "/", "pQTL与",
                          outname, "的MR结果.csv")
          Bname <- paste0(savefile, "/", "pQTL群与",
                          outname, "的异质性结果.csv")
          Cname <- paste0(savefile, "/", "pQTL与",
                          outname, "的多效性结果.csv")
          Dname <- paste0(savefile, "/", "pQTL群与",
                          outname, "的SNPs情况.csv")
          write.csv(A_temp, Aname, row.names = F)
          write.csv(B_temp, Bname, row.names = F)
          write.csv(C_temp, Cname, row.names = F)
          write.csv(D_temp, Dname, row.names = F)}else{E_temp<-rbind(i,E_temp)
          write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}}

      else{E_temp<-rbind(i,E_temp)
      write.csv(E_temp, "未匹配到IVs.csv", row.names = F)}
    }
  }
  cat("当前分析已完成")
}
