#' Subject-linked longitudinal infant faecal microbiome, mother milk microbiome and mother milk metabolome.
#'
#' A subset of the MAINHEALTH project data.
#'
#' @format ## `Jakobsen2025`
#' A list object with the following information:
#' \describe{
#'   \item{Zbmi}{Preprocessed and prepared data object for pre-pregnancy BMI.}
#'   \item{Zwhz}{Preprocessed and prepared data object for infant 6-month WHZ.}
#'   \item{subjectMeta_BMI}{Homogenized subject metadata that is ordered the same as the rows in each data block in Zbmi.}
#'   \item{subjectMeta_WHZ}{Homogenized subject metadata that is ordered the same as the rows in each data block in Zwhz.}
#'   \item{infant_faecal_microbiome_taxonomy}{Taxonomy of the processed ASVs, ordered the same way as the columns in block 1.}
#'   \item{milk_microbiome_taxonomy}{Taxonomy of the processed ASVs, ordered the same way as the columns in block 2.}
#'   \item{milk_metabolomics_features}{Metabolite names, ordered the same way as the columns in block 3.}
#' }
#' @source <https://bmjopen.bmj.com/content/12/11/e059552>
"Jakobsen2025"
