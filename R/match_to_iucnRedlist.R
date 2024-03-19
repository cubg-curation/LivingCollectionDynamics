#' Match collection to IUCN Red List of Threatened Species via taxonomic name
#'
#' @param collection A data frame containing a collection.
#' @param iucnRedlist IUCN Red List of Threatened Species database, obtained using the function XXXX.
#' @param taxon_name_column The name of the column in the `collection` corresponding to taxonomic names.
#' @param taxon_name_full_column The name of the column in the `collection` corresponding to joined taxonomic names and authors.
#' @param taxon_author_column The name of the column in the `collection` corresponding to the authors of the taxonomic names.
#' @param typo_method Either `'All'`, `'Data frame only'`,`'Data frame + common'`, detailing the level of typo finding required.
#' @param do_add_split Flag (TRUE/FALSE) for whether we search for missing f./var./subsp.
#' @param do_fix_hybrid Flag (TRUE/FALSE) for whether we search for hybrid issues.
#' @param do_rm_autonym Flag (TRUE/FALSE) for whether we try removing autonyms.
#' @param ... Arguments (i.e., attributes) used in the matching algorithm (passed along to nested functions). Examples include `enrich_taxon_authors_column`, `enrich_display_in_message_column` and `enrich_plant_identifier_column`.
#' @param enrich_taxon_name_column The name of the column in the `iucnRedlist` corresponding to taxonomic names.Default value is `scientific_name`.
#' @param enrich_taxon_authors_column The name of the column in `enrich_database` that corresponds to the authors of taxonomic names. Default value is `sanitise_author`.
#' @param enrich_display_in_message_column The name of the column in `iucnRedlist` that contains values to show in the matching messages. Default value is `taxonid`.
#' @param enrich_plant_identifier_column The name of the column in `iucnRedlist` that corresponds to record identifier. Default value is `taxonid`.
#' @param try_hybrid Flag (TRUE/FALSE) whether we want to look at hybrid fixes across all fixing methods.
#' @param matching_criterion A function used to chose the best method from extracts of the `iucnRedlist`.
#'
#' @details
#' This function allows matching of a collection's database to IUCN Red List of Threatened Species database. This function relies of matching functions found that are documented in  [match_single()], and is broadly similar to the function [match_collection_to_wcvp()].
#'
#' Note that by default the matching functions use column names from (WCVP) therefore these will often require changing prior to matching unless you change the column names in iucnRedlist to concure with wcvp column names.
#'
#' Within the algorithm there exists methods to improve the matching such as trying to change infraspecific levels (e.g var. to subsp.) or adding hybridisation. These methods can be turned on and off using `do_add_split`, `do_fix_hybrid`, `do_rm_autonym` and `typo_method`.
#'
#'
#' @return A list of length seven containing:
#' - `$match` the index of the record in `iucnRedlist` which matches the record in the collection database.
#' - `$details_short` a simplified message detailing the match.
#' - `$match_taxon_name` a longer format message detailing the match.
#' - `$original_authors` The author/s (extracted) from the `collection` database.
#' - `$match_authors` The author/s of the matched record in `iucnRedlist`.
#' - `$author_check` Either `Identical`, `Partial` or `Different`  (`No Match` if a match to iucnRedlist cannot be found). A message informing the similarity of the collection's taxon authors and the authors found in `iucnRedlist`. Author similarity is found using the function  [author_check()].
#'
#'
#' match = taxon_match_full,
#' @export
match_collection_to_iucnRedlist <- function(collection, iucnRedlist,
                                     taxon_name_column = 'TaxonName',
                                     taxon_name_full_column = NA,
                                     taxon_author_column = NA,
                                     typo_method = 'All',
                                     do_add_split = TRUE, do_fix_hybrid = TRUE,
                                     do_rm_autonym = TRUE,
                                     ...,
                                     enrich_taxon_name_column = 'scientific_name',
                                     enrich_taxon_authors_column = 'sanitise_author',
                                     enrich_display_in_message_column = 'taxonid',
                                     enrich_plant_identifier_column = 'taxonid',
                                     matching_criterion = LivingCollectionDynamics::no_additional_matching,
                                     try_hybrid = FALSE
                                    ){
  if(!typo_method %in% c('All', 'Data frame only','Data frame + common', 'no')){
    stop('Invalid typo_method input!')
  }
  #Implies collection and iucnRedlist are already in the workspace.
  no_records = nrow(collection)
  cli::cli_alert_info("{.var {no_records}} records found.")

  ################################################
  # 1) Setup original report. (only look at unique taxon name / taxon name full and add is_autonym)
  ################################################
  # Get taxon name and taxon author out of the original report.
  cli::cli_h2("Extracting taxon names and authors from the original report")

  taxon_name = collection[,match(taxon_name_column, names(collection))]
  do_taxon_author = TRUE
  if(!is.na(taxon_author_column)){
    taxon_author = collection[,match(taxon_author_column, names(collection))]
    if(all(taxon_author == '')){
      do_taxon_author = FALSE
    }
  }
  else if(!is.na(taxon_name_full_column)){
    taxon_name_full = collection[,match(taxon_name_full_column, names(collection))]
    taxon_author = author_from_taxon_name_full(taxon_name, taxon_name_full)
    if(all(taxon_author == '')){
      do_taxon_author = FALSE
    }
  }
  else{
    taxon_author = rep('', length(taxon_name))
    do_taxon_author = FALSE
  }
  if(!do_taxon_author){
    #Simplify taxon_author names.
    cli::cli_h2("Taxon authors not provided.")
    cli::cli_h2("Reducing to unique taxon names")
  }
  else{
    taxon_author = stringi::stri_trans_general(taxon_author, id = "Latin-ASCII")
    cli::cli_h2("Reducing to unique taxon name and author combinations")
  }

  # Restrict to only unique taxon name and author combinations.
  taxon_name_and_author = data.frame(taxon_name = taxon_name, taxon_author = taxon_author)
  unique_taxon_name_and_author =unique(taxon_name_and_author)
  taxon_name_and_author_combined = do.call(paste, c(taxon_name_and_author, sep='-'))
  unique_taxon_name_and_author_combined = do.call(paste, c(unique_taxon_name_and_author, sep='-'))

  report_match = match(taxon_name_and_author_combined,unique_taxon_name_and_author_combined)

  taxon_name = unique_taxon_name_and_author$taxon_name
  taxon_author =  unique_taxon_name_and_author$taxon_author

  no_unique = length(taxon_name)
  cli::cli_alert_info("{.var {no_unique}} unique taxon names/ taxon name author combinations found.")


  ################################################
  # 2) Setup outputs.
  ################################################
  taxon_match_full = rep(NA,nrow(collection))
  taxon_name_story_full = rep(NA,nrow(collection))
  taxon_match = rep(NA, length(taxon_name))
  taxon_name_story = taxon_name
  index_to_find_matches = 1:length(taxon_name)
  index_complete = NULL


  ################################################
  # 3) Remove cultivars and indeterminates. (set taxon match to -1)
  ################################################
  cli::cli_h2("Removing cultivars and indeterminates from {length(index_to_find_matches)} name{?s}")

  match_info = no_match_cultivar_indet(taxon_name[index_to_find_matches])

  taxon_match[index_to_find_matches] = match_info$match
  taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

  index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
  no_found = length(index_to_find_matches[!is.na(match_info$match)])
  cli::cli_alert_success("Found {no_found} cultivars and indeterminates")

  index_to_find_matches = index_to_find_matches[is.na(match_info$match)]

  ################################################
  # 5) Match original report to all unique taxon names in POWO.
  ################################################
  if(length(index_to_find_matches) > 0){
    cli::cli_h2("Matching {length(index_to_find_matches)} name{?s} to unique taxon names")

    single_indices = which(iucnRedlist$single_entry == TRUE)
    match_info = match_single(taxon_names = taxon_name[index_to_find_matches],
                              enrich_database = iucnRedlist,
                              enrich_database_search_index = single_indices,
                              enrich_taxon_name_column = enrich_taxon_name_column,
                              enrich_display_in_message_column = enrich_display_in_message_column,
                              ...)

    taxon_match[index_to_find_matches] = match_info$match
    taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

    index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
    no_found = length(index_to_find_matches[!is.na(match_info$match)])
    cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")

    index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
  }

  ################################################
  # 6) Match original report to all taxon names with a multiple entry in iucnRedlist.
  ################################################
  if(length(index_to_find_matches) > 0){
    cli::cli_h2("Matching {length(index_to_find_matches)} name{?s} to non-unique taxon names")

    mult_indices = which(iucnRedlist$single_entry == FALSE)
    match_info = match_multiple(taxon_names = taxon_name[index_to_find_matches],
                                taxon_authors = taxon_author[index_to_find_matches],
                                enrich_database = iucnRedlist,
                                enrich_database_search_index = mult_indices,
                                enrich_taxon_name_column = enrich_taxon_name_column,
                                enrich_taxon_authors_column = enrich_taxon_authors_column,
                                enrich_display_in_message_column = enrich_display_in_message_column,
                                enrich_plant_identifier_column = enrich_plant_identifier_column,
                                matching_criterion = matching_criterion,
                                ...)

    taxon_match[index_to_find_matches] = match_info$match
    taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

    index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
    no_found = length(index_to_find_matches[!is.na(match_info$match)])
    cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")

    index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
  }

  ################################################
  # 7) If the matched author is different after exact matching try to find a similar taxon name.
  ################################################
  if(do_taxon_author){
    # Find the records which have disagreeing authors.
    with_match = index_complete[which(!grepl('Do not attempt matching',taxon_name_story[index_complete]))]

    # Do we have at least one match.
    if(length(with_match) > 0){
      original_authors = taxon_author[with_match]
      matched_authors = rep(NA,length(original_authors))
      match_bigger_0 = taxon_match[with_match] > 0
      matched_authors[match_bigger_0] = iucnRedlist[[enrich_taxon_authors_column]][taxon_match[with_match][match_bigger_0]]

      author_checked = rep(NA, length(original_authors))
      for(i in 1:length(author_checked)){
        author_checked[i] = author_check(original_authors[i],  matched_authors[i])
      }

      # Get the indices of those matched to a single record with author checked = 'Different'.
      diff_index = with_match[which(author_checked == 'Different')]

      if(length(diff_index) > 0){
        # Get matches trying to fix taxon name.
        match_info = match_all_issue(taxon_names = taxon_name[diff_index],
                                     taxon_authors = taxon_author[diff_index],
                                     enrich_database = iucnRedlist,
                                     do_add_split = do_add_split,
                                     do_fix_hybrid = do_fix_hybrid,
                                     do_rm_autonym = do_rm_autonym,
                                     matching_criterion = matching_criterion,
                                     enrich_taxon_name_column = enrich_taxon_name_column,
                                     enrich_taxon_authors_column = enrich_taxon_authors_column,
                                     enrich_plant_identifier_column = enrich_plant_identifier_column,
                                     enrich_display_in_message_column = enrich_display_in_message_column,
                                     try_hybrid = try_hybrid,
                                     ...)

        proposed_authors = rep(NA, length(taxon_author[diff_index]))
        current_match = match_info$match
        current_match[is.na(current_match)] = -3
        proposed_authors[current_match > 0] = iucnRedlist[[enrich_taxon_authors_column]][current_match[current_match >0]]

        # Check if proposed author is same/similar to original.
        compare_author_new = rep(NA, length(proposed_authors))
        for(i in 1:length(diff_index)){
          compare_author_new[i] = author_check(taxon_author[diff_index[i]],  proposed_authors[i])
        }
        compare_author_new[current_match < 0] = 'No Match'

        improved_author_index = which(compare_author_new %in% c('Identical', 'Partial'))
        if(length(improved_author_index)>0){
          taxon_match[diff_index[improved_author_index]] = match_info$match[improved_author_index]
          taxon_name_story[diff_index[improved_author_index]] =
            paste0(taxon_name_story[diff_index[improved_author_index]],
                   ' -> (Author differ)',
                   match_info$message[improved_author_index])
        }
      }
    }

  }


  ################################################
  # 8) Try fixing taxon name issues.
  ################################################
  if(length(index_to_find_matches) > 0){
    cli::cli_h2("Testing and matching taxon name issues for {length(index_to_find_matches)} name{?s}")

    match_info = match_all_issue(taxon_names = taxon_name[index_to_find_matches],
                                 taxon_authors = taxon_author[index_to_find_matches],
                                 enrich_database = iucnRedlist,
                                 do_add_split = do_add_split,
                                 do_fix_hybrid = do_fix_hybrid,
                                 do_rm_autonym = do_rm_autonym,
                                 matching_criterion = matching_criterion,
                                 enrich_taxon_name_column = enrich_taxon_name_column,
                                 enrich_taxon_authors_column = enrich_taxon_authors_column,
                                 enrich_display_in_message_column = enrich_display_in_message_column,
                                 enrich_plant_identifier_column = enrich_plant_identifier_column,
                                 try_hybrid = try_hybrid,
                                 ...)

    taxon_match[index_to_find_matches] = match_info$match
    taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

    index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
    no_found = length(index_to_find_matches[!is.na(match_info$match)])
    cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")

    index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
  }

  ################################################
  # 9) Try to find typo and then match.
  ################################################
  if(typo_method %in%  c('All', 'Data frame only','Data frame + common') & length(index_to_find_matches) > 0){
    cli::cli_h2("Testing and matching typos for {length(index_to_find_matches)} name{?s}")

    match_info = match_typos(taxon_names = taxon_name[index_to_find_matches],
                             taxon_authors = taxon_author[index_to_find_matches],
                             enrich_database = iucnRedlist,
                             single_indices = single_indices,
                             mult_indices = mult_indices,
                             typo_method = typo_method,
                             enrich_taxon_name_column = enrich_taxon_name_column,
                             enrich_taxon_authors_column = enrich_taxon_authors_column,
                             enrich_display_in_message_column = enrich_display_in_message_column,
                             enrich_plant_identifier_column = enrich_plant_identifier_column,
                             matching_criterion = matching_criterion,
                             ...)

    taxon_match[index_to_find_matches] = match_info$match
    taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

    index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
    no_found = length(index_to_find_matches[!is.na(match_info$match)])
    cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")


    index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
  }

  ################################################
  # 10) Check if the authors match from the name match (before going to accepted name)
  ################################################
  # Get the proposed authors.
  proposed_authors = rep(NA, length(taxon_author))
  current_match = taxon_match
  current_match[is.na(current_match)] = -3
  proposed_authors[current_match > 0] = iucnRedlist[[enrich_taxon_authors_column]][current_match[current_match >0]]
  matched_name = rep(NA, length(taxon_author))
  matched_name[current_match > 0] = iucnRedlist[[enrich_taxon_name_column]][current_match[current_match >0]]

  # Check if proposed author is same/similar to original.
  author_checked = rep(NA, length(proposed_authors))
  for(i in 1:length(author_checked)){
    author_checked[i] = author_check(taxon_author[i],  proposed_authors[i])
  }
  author_checked[current_match < 0] = 'No Match'

  ################################################
  # 12) Set remaining taxon_match to -3 and add story.
  ################################################
  taxon_match[index_to_find_matches] = -3
  author_checked[index_to_find_matches] = 'No Match'
  taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], ' -> (No match found)')

  ################################################
  # 13) Create a shortened version on the match details
  ################################################
  match_short = shorten_message(taxon_name_story)

  ################################################
  # 15) return match (to original report) and details of the matches (for unique plants).
  ################################################
  cli::cli_h2("Matching Complete")
  taxon_match_full = taxon_match[report_match]
  taxon_name_story_full = taxon_name_story[report_match]
  match_short_full = match_short[report_match]
  matched_name_full = matched_name[report_match]
  proposed_authors_full = proposed_authors[report_match]
  author_checked_full = author_checked[report_match]
  original_authors_full = taxon_author[report_match]
  return(list(match = taxon_match_full,
              details = taxon_name_story_full,
              details_short = match_short_full,
              match_taxon_name = matched_name_full,
              original_authors = original_authors_full,
              match_authors = proposed_authors_full,
              author_check = author_checked_full))
}

