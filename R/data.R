#' Expression data from breast cancer patients
#'
#' This matrix contains a subset of RNA-seq gene expression data from patients
#' with breast cancer. The full data is available in the Gene Expression
#' Omnibus (GEO) repository under accession number GSE96058. The format of the
#' expression value is log2(FPKM + 0.1) and was rounded to one digit to reduce
#' the data size. The corresponding sample information is stored in the
#' data object "samplesBC".
#'
#' @format A matrix with 12,787 rows (genes) and 400 columns (samples).
#' @usage data(exprBC)
#'
"exprBC"

#' Phenotype data from breast cancer patients.
#'
#' This data frame contains a subset of phenotype data from patients with
#' breast cancer. Full data are available in the Gene Expression Omnibus
#' (GEO) repository under accession number GSE96058. The corresponding
#' expression data is stored in the data object "exprBC".
#'
#' @format A data frame containing 11 variables:
#' \describe{
#'  \item{title}{Sample title.}
#'  \item{age_at_diagnosis}{Age at diagnosis.}
#'  \item{chemo_treated}{Whether the patient received adjuvant chemotherapy
#'  (1 for "yes", 0 for "no").}
#'  \item{endocrine_treated}{Whether the patient received adjuvant endocrine
#'  therapy (1 for "yes", 0 for "no").}
#'  \item{er_status}{Estrogen receptor status (1 for "positive", 0 for
#'  "negative").}
#'  \item{her2_status}{Human epidermal growth factor receptor 2 (HER2) status
#'  (1 for "positive", 0 for "negative").}
#'  \item{pgr_status}{Progesterone receptor status (1 for "positive", 0 for
#'  "negative").}
#'  \item{nhg}{Nottingham histological grade.}
#'  \item{pam50_subtype}{PAM50 subtype.}
#'  \item{overall_survival_days}{Days of overall survival.}
#'  \item{overall_survival_event}{Overall survival event (1 for "deceased", 0
#'  for "alive").}
#' }
#' @usage data(samplesBC)
#'
"samplesBC"

#' A subset of genes involved in early estrogen response
#'
#' This vector contains a subset of genes involved in early estrogen response.
#'
#' @format A character vector
#' @usage data(genesER)
"genesER"
