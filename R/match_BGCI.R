#' match_original_to_BGCI()
#'
#' @param original_report  original_report
#' @param BGCI  BGCI
#' @param taxon_name_col  taxon_name_col
#'
#' @return no gardens
#' @export
#'
match_original_to_BGCI <- function(original_report, BGCI,
                                   taxon_name_col = 'TaxonName'){

  #Implies original_report and wcvp are already in the workspace.
  cli::cli_h1("Matching names to Garden Count (BGCI)")
  no_records = nrow(original_report)
  cli::cli_alert_info("{.var {no_records}} records found.")

  ################################################
  # 1) Setup original report. (only look at unique taxon name / taxon name full and add is_autonym)
  ################################################
  # Get taxon name and taxon author out of the original report.
  cli::cli_h2("Extracting taxon names and authors from the original report")

  taxon_name = original_report[,match(taxon_name_col, names(original_report))]

  # Restrict to only unique taxon name and author combinations.
  unique_taxon_name =unique(taxon_name)

  report_match = match(taxon_name,unique_taxon_name)

  taxon_name = unique_taxon_name

  no_unique = length(taxon_name)
  cli::cli_alert_info("{.var {no_unique}} unique taxon names found.")


  ################################################
  # 2) Setup outputs.
  ################################################
  taxon_match_full = rep(NA,nrow(original_report))
  taxon_name_story_full = rep(NA,nrow(original_report))
  taxon_match = rep(NA, length(taxon_name))
  taxon_name_story = taxon_name
  index_to_find_matches = 1:length(taxon_name)
  index_complete = NULL

  ################################################
  # 3) Sanitise taxon names.
  ################################################
  if(length(index_to_find_matches) > 0){
    cli::cli_h2("Sanitise taxon names")

    santise_taxon_name = unlist(pbapply::pblapply(taxon_name[index_to_find_matches],LivingCollectionDynamics::sanitise_name))
    original_taxon_name = taxon_name
    taxon_name[index_to_find_matches] = santise_taxon_name

    indices_require_sanitise=which(taxon_name != original_taxon_name)
    taxon_name_story[indices_require_sanitise] = paste0(taxon_name_story[indices_require_sanitise],
                                                        ' -> (Sanitise name) -> ',
                                                        taxon_name[indices_require_sanitise])
    cli::cli_alert_success("Sanitising required for {length(indices_require_sanitise)} taxon names")
  }

  ################################################
  # 3) Match to BGCI.
  ################################################
  taxon_match = match(taxon_name, BGCI$full_name)


  cli::cli_h2("Matching Complete")
  taxon_match_full = taxon_match[report_match]
  no_gardens = BGCI$count.gardenid.[taxon_match_full]
  return(list(no_gardens = no_gardens))
}
