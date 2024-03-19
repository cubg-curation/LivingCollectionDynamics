#'  Matching criterions
#'
#'  Custom matching criterions used within `match_all_issue()` and `get_match_from_multiple()`.
#'
#' Within BGSmartR we have two in-built custom matching criterions.
#'
#'  - `no_additional_matching()` is used when no custom matching wants to be performed.
#'
#'  - `additional_wcvp_matching()` is the default for additional matching to the World Checklist of Vascular Plants (WCVP) database. This method includes looking ahead to see if potential matches lead to the same accepted plant name, and matching on the taxon_status of the records (Accepted or Synonym).
#'
#'  These can be used within [match_all_issue()] and [get_match_from_multiple()] via the input `matching_criterion`. See vignette `Method of Matching taxonomic records` for further details on custom matching and how to create new matching criterions.
#'
#' @param enrich_database_extract extract of enrichment database (often records with identical taxonomic names).
#' @param message Matching message.
#'
#' @return Always returns a list of 2 with:
#' - `$row` The rows of `enrich_database_extract` corresponding to the best match/s according to the criterion.
#' - `$message` A combination of the input message the method used to determine a match.
#' @export
#'
no_additional_matching <- function(enrich_database_extract, message = ''){
  # No matching criteria
  matched = 1:nrow(enrich_database_extract)
  message = paste0('(', message, 'unclear, do not match',')',collapse = '')

  return(list(row = matched, message = message))
}


#' @rdname no_additional_matching
#' @export
additional_wcvp_matching <- function(enrich_database_extract, message = ''){

  #Check that the desired columns used exist in enrich_database_extract otherwise throw an error.
  if(!all(c('plant_name_id', 'accepted_plant_name_id','taxon_status') %in% colnames(enrich_database_extract))){
    stop('Error in choose_best_from_enriched_database required columns not in enrich_database!')
  }
  match_flag = FALSE

  accepted_plant_id = enrich_database_extract$accepted_plant_name_id

  # Check if all exact matches point to the same accepted name
  if(identical(accepted_plant_id, rep(accepted_plant_id[1], length(accepted_plant_id)))){
    # match to accepted if one exists if not the first plant that matches.
    taxon_accepted = enrich_database_extract$taxon_status == 'Accepted'
    if(any(taxon_accepted)){
      matched = which(taxon_accepted)[1]
      message = paste0('(', message, 'all point to same accepted plant',')',collapse = '')
    }
    else{
      matched = 1
      message = paste0('(', message, 'all point to same accepted plant',')',collapse = '')

    }
    match_flag = TRUE
  }

  # Check for differences in taxon_status.
  if(!match_flag){
    taxon_status = enrich_database_extract$taxon_status
    taxon_status_match = match(taxon_status , c('Accepted', 'Synonym'))
    if(all(is.na(taxon_status_match))){
      matched = 1:nrow(enrich_database_extract)
      message = paste0('(', message, 'no accepted or synonym',')',collapse = '')
      match_flag = TRUE
    }
    else{
      chosen_record = which(taxon_status_match == min(taxon_status_match, na.rm = T))
      if(length(chosen_record) == 1){
        matched = chosen_record
        message = paste0('(', message, 'choose via taxon_status',')',collapse = '')
        match_flag = TRUE
      }
      else if(length(chosen_record) >1){
        matched = chosen_record
        message = paste0('(', message, 'multiple best taxon status',')',collapse = '')
        match_flag = TRUE
      }
    }

  }

  return(list(row = matched, message = message))
}
