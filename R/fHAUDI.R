
# Creates multiplicative penalty vector based on annotation data
anno_phi <- function(anno_info_df, tg_param){

  anno_mat <- as.matrix(anno_info_df[!(anno_info_df$anc_ref),-c(1:7)]) #make it into matrix to make it easier to work with
  #row.names(anno_mat) <- anno_info_df %>% filter(!anc_ref) %>% pull("rsid")

  if (n_anno == 1){
    anno_mult <- ((1 - anno_mat) %*% tg_param[-1]) + anno_mat # lambda_s * (1 - R_js) + R_js
  } else{
    anno_mult <- (1 - anno_mat) %*% diag(tg_param[-1]) + anno_mat # lambda_s * (1 - R_js) + R_js
  }

  gamma_vec <- anno_info_df %>% filter(!anc_ref) %>% mutate(gamma = ifelse(anc == "all", 1, tg_param[1])) %>% pull(gamma)
  anno_mult <-  cbind(anno_mult, gamma_vec) #add gamma line
  anno_phi_vec <- apply(anno_mult, 1, prod) # product for each snp

  return(anno_phi_vec)
}


##' Fit fHAUDI model
##'
##' @title fHAUDI
##' @param fbm_obj object of FBM class
##' @param fbm_anno_info info data frame containing information FBM (chrom/pos/samples/etc.) merged with annotation data
##' @param y Vector of responses
##' @param gamma_anno_s Vector specifying multiplicative penalty for gamma and then annotation specific penalties
##' @param ind_train Vector of indices specifying the rows to use for training the model
##' @param family Either "gaussian" or "binomial" (default is gaussian)
##' @param snps Vector of SNPs to include in model
##' @param ... additional arguments to pass to big_spLinReg / big_spLogReg
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Brian Chen
##' @import bigstatsr
##' @export

fHAUDI <- function(fbm_obj, fbm_anno_info, y, gamma_anno_s, family="gaussian", ind_train = NULL, snps = NULL, ...) {


  ## specify training data
  if (is.null(ind_train)) {
    ind_train <- seq_len(nrow(fbm_obj))
  }

  ## specify SNPs to retain
  if (is.null(snps)) {
    col_keep <- rep(TRUE, ncol(fbm_obj))
  } else {
    col_keep <- fbm_anno_info$rsid %in% snps
  }
  ## remove reference-ancestry columns
  col_keep <- col_keep & (!fbm_anno_info$anc_ref)
  ind_col <- which(col_keep) # subset SNPs


  #multiplicative penalty
  anno_phi(anno_info, gamma_anno_s)

  if (family == "gaussian") {
    model <- bigstatsr::big_spLinReg(
      X = fbm_obj,
      y.train = y[ind_train],
      ind.train = ind_train,
      pf.X = pf_x,
      ind.col = ind_col,
      ...
    )
  } else if (family == "binomial") {
    model <- bigstatsr::big_spLogReg(
      X = fbm_obj,
      y01.train = y[ind_train],
      ind.train = ind_train,
      pf.X = pf_x,
      ind.col = ind_col,
      ...
    )
  }


  return(model)

}



##' Extract population-specific SNP coefficients from a HAUDI model
##'
##' @title get_beta_haudi
##' @param fbm_anno_info data frame containing information
##' for FBM (chrom/pos/samples/etc.)
##' @param haudi_model an object of class big_sp_list,
##' returned by `haudi` or `lasso`
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @import bigstatsr
##' @import data.table
##' @export
get_beta_fhaudi <- function(fbm_anno_info, haudi_model) {
  anc <- snp <- beta_all <- `:=` <- NULL # due to R CMD check

  dt_snp <- data.table::data.table(
    snp = fbm_anno_info$rsid[attr(haudi_model, "ind.col")],
    beta = summary(haudi_model)$beta[[1]],
    anc = fbm_anno_info$anc[attr(haudi_model, "ind.col")]
  )

  ancestries <- unique(fbm_anno_info$anc)
  ancestries <- ancestries[ancestries != "all"]

  dt_ref <- data.table::data.table(dt_snp[dt_snp$anc == "all", ])
  dt_ref$anc <- NULL
  dt_ref$beta_all <- dt_ref$beta

  for (ancestry in ancestries) {
    x <- dt_snp[dt_snp$anc == ancestry, ]
    idx <- match(dt_ref$snp, x$snp)
    beta_diff <- rep(0, length(idx))
    beta_diff[!is.na(idx)] <- x[idx[!is.na(idx)], ]$beta
    new_col <- paste0("beta_", ancestry)
    dt_ref <- data.table::set(dt_ref,
                              j = new_col,
                              value = dt_ref$beta_all + beta_diff
    )
  }

  dt_ref$beta_all <- dt_ref$beta <- NULL
  return(dt_ref)
}



