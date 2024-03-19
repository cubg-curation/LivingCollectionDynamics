#' Prepare enrichment database for matching
#'
#' @details
#' This function adds columns to the enrichment database used when matching to taxonomic names. By default this includes:
#'
#' - `sanitise_name` The sanitised taxonomic name.
#' - `sanitise_author` The sanitised author of the taxonomic name.
#' - `require_sanitise` A logical column (TRUE/FALSE) for whether the taxonomic name required sanitising.
#' = `author_parts` Words found in the taxonomic author (removing initials and punctuation) used when performing partial author matching.
#' - `taxon_length` The length of the taxonomic name (string length) used if typo searching if performed in the matching algorithm.
#' - `single_entry` A logical column (TRUE/FALSE) for whether the taxonomic name is unique in `enrich_database`, used to restrict records when performing either single matching or multiple matching.
#' - `ID` A column of unique identifiers for each record in `enrich_database`. Used for referencing the matches.
#'
#' Moreover the function sorts the enrich database in alphabetical order of the taxonomic names.
#'
#' These additional column can be switched on/off using the inputs  `do_sanitise`, `do_taxon_length`, `do_single_entry`, `do_author_parts` and `do_sort`.
#'
#' Note that if sanitising if performed then sorting and single_entry will be performed with the sanitised taxonoic names and not the inputted names (`enrich_database[[enrich_taxon_name_column]]`). Similarly, sanitised authors will be used when creating author parts.
#'
#' @param enrich_database A data frame of enriching information.
#' @param enrich_taxon_name_column The name of the column in `enrich_database` that corresponds to taxonomic names. Default value is `taxon_names`.
#' @param enrich_taxon_authors_column The name of the column in `enrich_database` that corresponds to the authors of taxonomic names.
#' @param enrich_taxon_name_full_column The name of the column in the `enrich_database` corresponding to joined taxonomic names and authors.
#' @param do_sanitise  Flag (TRUE/FALSE) detailing whether to add the columns `sanitise_name`, `sanitise_author` and  `require_sanitise` corresponding to the sanitised taxonomic name, sanitised author and a flag (TRUE/FALSE) of whether the taxonomic name needed sanitising.
#' @param do_taxon_length Flag (TRUE/FALSE) detailing whether to add the column `taxon_length` containing the string length of the taxonomic names.
#' @param do_single_entry Flag (TRUE/FALSE) detailing whether to add the column `single_entry` containing whether the taxonomic name appears multiple times in `enrich_database`.
#' @param do_author_parts Flag (TRUE/FALSE) detailing whether to add the column `author_parts` containing the words extracted from the taxonomic author. This is used for partial author matching.
#' @param do_add_id Flag (TRUE/FALSE) detailing whether to add the column `ID` containing a unique identifer for each record in the `enrich_database`.
#' @param do_sort Flag (TRUE/FALSE) detailing whether to  alphabetically sort the taxonomic names in  `enrich database`.
#' @param console_message Flag (TRUE/FALSE) detailing whether to show messages in the console.
#'
#'
#' @return `enrich_database` with additional columns used by matching to taxonomic names.
#' @export
#'
#' @examples
#' taxon_names = c('Abies taxifolia', 'ABIES taxifolia',
#'  'Acalypha gracilens', 'Eupatorium magdalenae', 'Adina racemosa')
#' taxon_authors = c('Duhamel', 'Poir.', 'A.Gray', 'Stehlé', '(Siebold & Zucc.) Miq.')
#'
#' enrich_database = data.frame(taxon_names,taxon_authors, value = runif(5))
#' prepare_enrich_database(enrich_database,
#' enrich_taxon_name_column = 'taxon_names',
#' enrich_taxon_authors_column = 'taxon_authors')
prepare_enrich_database <- function(enrich_database,
                                        enrich_taxon_name_column = 'taxon_names',
                                        enrich_taxon_authors_column = NA,
                                        enrich_taxon_name_full_column = NA,
                                        do_sanitise = TRUE,
                                        do_taxon_length = TRUE,
                                        do_single_entry = TRUE,
                                        do_author_parts = TRUE,
                                        do_add_id = TRUE,
                                        do_sort = TRUE,
                                        console_message = FALSE){

  ### A) Sanitising taxonomic names and authors.
  if(do_sanitise){
    if(console_message){
      cli::cli_h1("Sanitising taxonomic names and author")
    }
    sanitised = LivingCollectionDynamics::sanitise_names_authors_report(
      enrich_database,
      taxon_name_column = enrich_taxon_name_column,
      taxon_name_full_column = enrich_taxon_name_full_column,
      taxon_author_column = enrich_taxon_authors_column,
      console_message = console_message)

    enrich_database$sanitise_name = sanitised$taxon_name
    enrich_database$sanitise_author = sanitised$author
    enrich_database$require_sanitise = sanitised$sanitised
  }

  ### B) Add taxon name length column.
  if(do_taxon_length){
    if(console_message){
      cli::cli_h1("Add taxon name length column")
      taxon_length = unlist(pbapply::pblapply(enrich_database[[enrich_taxon_name_column]], stringr::str_length))
    }else{
      taxon_length = unlist(lapply(enrich_database[[enrich_taxon_name_column]], stringr::str_length))
    }

    enrich_database$taxon_length = taxon_length

  }

  ### C)  Add column detailing whether the taxonomic name is a single entry.
  if(do_single_entry){
    if(console_message){
      cli::cli_h1("Add single entry column")
    }
    ### Decide whether to use the taxonomic name column or the sanitised name column.
    if('sanitise_name' %in% names(enrich_database)){
      wanted_column = 'sanitise_name'
    }else{
      wanted_column = enrich_taxon_name_column
    }

    name_freq = table(enrich_database[[wanted_column]])
    single_entry = rep(NA, nrow(enrich_database))
    single_entry[enrich_database[[wanted_column]] %in% names(name_freq)[as.numeric(name_freq) == 1]] = TRUE
    single_entry[enrich_database[[wanted_column]] %in% names(name_freq)[as.numeric(name_freq) > 1]] = FALSE
    enrich_database$single_entry =  single_entry
  }

  ### D) Add column containing words from the taxonomic author.
  if(do_author_parts){

    ### D.1) Decide whether we use the inputted author or extracted authors.
    if('sanitise_author' %in% names(enrich_database)){
      wanted_author_column = 'sanitise_author'
    }else{
      wanted_author_column = enrich_taxon_authors_column
    }

    ### D.2) Use author words to get the words from each author.
    if(console_message){
      cli::cli_h1("Add author words column")
      author_word = unlist(pbapply::pblapply(enrich_database[[wanted_author_column]],function(x){paste0(LivingCollectionDynamics::author_words(x),collapse =', ')}))

    }else{
      author_word = unlist(lapply(enrich_database[[wanted_author_column]],function(x){paste0(LivingCollectionDynamics::author_words(x),collapse =', ')}))
    }

    enrich_database$author_parts = author_word

  }

  ### E) sort the records into alphabetical order.
  if(do_sort){
    if(console_message){
      cli::cli_h1("Sort the records into alphabetical order")
    }
    enrich_database = enrich_database[order(enrich_database[[wanted_column]]),]
  }

  ### F) Add identifer column 'ID'.
  if(do_add_id){
    enrich_database = data.frame(ID = 1:nrow(enrich_database), enrich_database)
  }
  return(enrich_database)
}

#' Match collection to an enrichment database via taxonomic names
#'
#' @param collection A data frame containing a collection.
#' @param enrich_database A data frame of enriching information.
#' @param taxon_name_column The name of the column in the `collection` corresponding to taxonomic names.
#' @param taxon_name_full_column The name of the column in the `collection` corresponding to joined taxonomic names and authors.
#' @param taxon_author_column The name of the column in the `collection` corresponding to the authors of the taxonomic names.
#' @param typo_method Either `'All'`, `'Data frame only'`,`'Data frame + common'`, `no`; detailing the level of typo finding required.
#' @param do_rm_cultivar_indeterminates Flag (TRUE/FALSE) for whether we remove cultivars and indeterminates prior to taxonomic name matching.
#' @param do_match_single Flag (TRUE/FALSE) for whether we do matching to unique taxonomic names in `enrich_database`.
#' @param do_match_multiple Flag (TRUE/FALSE) for whether we do matching to non-unique taxonomic names in `enrich_database`.
#' @param do_fix_taxon_name Flag (TRUE/FALSE) for whether attempt to fix common issues in taxonomic names to aid matching. Sections of common issue fixes can also be turned on/off using the inputs `do_add_split`, `do_fix_hybrid`, `do_rm_autonym`.
#' @param do_add_split Flag (TRUE/FALSE) for whether we search for missing f./var./subsp.
#' @param do_fix_hybrid Flag (TRUE/FALSE) for whether we search for hybrid issues.
#' @param do_rm_autonym Flag (TRUE/FALSE) for whether we try removing autonyms.
#' @param ... Arguments (i.e., attributes) used in the matching algorithm (passed along to nested fuctions). Examples include, `enrich_display_in_message_column` and `enrich_plant_identifier_column`.
#' @param matching_criterion A function used to chose the best method from extracts of the `enrich_database`.
#' @param enrich_taxon_authors_column The name of the column in `enrich_database` that corresponds to the authors of taxonomic names.
#' @param enrich_taxon_name_column The name of the column in `enrich_database` that corresponds to the taxonomic names.

#'
#' @details
#' This function allows matching of a collection's database to an enrichment database.
#'
#' By default the function uses all the steps of our matching algorithm, for details of this see the vignette `Matching.Rmd` (Method of Matching taxonomic records). If parts of the algorithm are not required these can be switched off using  `typo_method`, `do_add_split`, `do_fix_hybrid`, `do_rm_autonym`, `do_rm_cultivar_indeterminates`, `do_match_single`, `do_match_multiple` and `do_fix_taxon_name.` Moreover, by default no custom matching if performed. A user inputted custom matching criterion (function) can be added via the input `matching_criterion`.
#'
#' To perform the matching you must specify the columns name of the taxon name in the enrichment database (`enrich_taxon_name_column`). If author matching is required then this column must also be specified for the enrichment database (`enrich_taxon_authors_column`).
#'
#'  The enrichment database must have some columns required for matching (`single_entry`, `taxon_length`, etc), we advice using [prepare_enrich_database()] to add these columns.
#'
#' Similarly, you must specify the columns name of the taxon name in the collection database (`taxon_name_column`). If author matching is desired then you have two choices:
#'
#' - specify the taxon author column `taxon_author_column`.
#' - Specify the combined taxon name and author column, `taxon_name_full_column` which removes words found in the taxon names from taxon names full to extract the authors.
#'
#' Note if both are specified then the authors from `taxon_author_column` are used.
#'
#' @return A list of length seven containing:
#' - `$match` the index of the record in `enrich_database` which matches the record in the collection database.
#' - `$details_short` a simplified message detailing the match.
#' - `$match_taxon_name` a longer format message detailing the match.
#' - `$original_authors` The author/s (extracted) from the `collection` database.
#' - `$match_authors` The author/s of the matched record in `enrich_database`.
#' - `$author_check` Either `Identical`, `Partial` or `Different`  (`No Match` if a match to `enrich_database` cannot be found). A message informing the similarity of the collection's taxon authors and the authors found in `enrich_database`. Author similarity is found using the function  [author_check()].
#'
#' @export
match_collection_to_enrich_database <- function(collection, enrich_database,
                                     taxon_name_column = NA,
                                     taxon_name_full_column = NA,
                                     taxon_author_column = NA,
                                     enrich_taxon_name_column = NA,
                                     enrich_taxon_authors_column = NA,
                                     typo_method = 'All',
                                     do_add_split = TRUE,
                                     do_fix_hybrid = TRUE,
                                     do_rm_autonym = TRUE,
                                     do_rm_cultivar_indeterminates = TRUE,
                                     do_match_single = TRUE,
                                     do_match_multiple = TRUE,
                                     do_fix_taxon_name = TRUE,
                                     matching_criterion = LivingCollectionDynamics::no_additional_matching,
                                     ...){
  ################################################
  # 1) Check we have required inputs.
  ################################################
  if(!typo_method %in% c('All', 'Data frame only','Data frame + common')){
    stop('Invalid typo_method input!')
  }
  if(!taxon_name_column %in% names(collection)){
    stop('"taxon_name_column" is not in "collection"!')
  }
  if(!taxon_name_column %in% names(collection)){
    stop('"enrich_taxon_name_column" is not in "enrich_database"!')
  }



  ################################################
  # 2) Setup original report. (only look at unique taxon name / taxon name full and add is_autonym)
  ################################################
  # 2.1) Get the number of records and print to console.
  no_records = nrow(collection)
  cli::cli_alert_info("{.var {no_records}} records found.")

  cli::cli_h2("Extracting taxon names and authors from the collection")

  # 2.2) Get the taxon names and authors from the collection.
  taxon_name = collection[[taxon_name_column]]
  do_taxon_author = TRUE
  if(!is.na(taxon_author_column)){
    if(!taxon_author_column %in% names(collection)){
      stop('"taxon_author_column" is not in "collection"!')
    }
    taxon_author = collection[[taxon_author_column]]
    if(all(taxon_author == '')){
      do_taxon_author = FALSE
    }
  }
  else if(!is.na(taxon_name_full_column)){
    if(!taxon_name_full_column %in% names(collection)){
      stop('"taxon_name_full_column" is not in "collection"!')
    }
    taxon_name_full = collection[[taxon_name_full_column]]
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
    cli::cli_h2("Taxon authors not provided.")
    cli::cli_h2("Reducing to unique taxon names")
  }
  else{
    taxon_author = stringi::stri_trans_general(taxon_author, id = "Latin-ASCII")
    cli::cli_h2("Reducing to unique taxon name and author combinations")
  }

  # 2.3) Restrict to only unique taxon name and author combinations.
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
  # 3) Setup outputs.
  ################################################
  taxon_match_full = rep(NA,nrow(collection))
  taxon_name_story_full = rep(NA,nrow(collection))
  taxon_match = rep(NA, length(taxon_name))
  taxon_name_story = taxon_name
  index_to_find_matches = 1:length(taxon_name)
  index_complete = NULL

  ################################################
  # 4) Remove cultivars and indeterminates. (set taxon match to -1)
  ################################################
  if(do_rm_cultivar_indeterminates){
    cli::cli_h2("Removing cultivars and indeterminates from {length(index_to_find_matches)} name{?s}")

    match_info = no_match_cultivar_indet(taxon_name[index_to_find_matches])

    taxon_match[index_to_find_matches] = match_info$match
    taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

    index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
    no_found = length(index_to_find_matches[!is.na(match_info$match)])
    cli::cli_alert_success("Found {no_found} cultivars and indeterminates")

    index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
  }

  ################################################
  # 5) Match collection to all single taxon names.
  ################################################
  if(do_match_single){
    if(length(index_to_find_matches) > 0){
      cli::cli_h2("Matching {length(index_to_find_matches)} name{?s} to unique taxon names")

      single_indices = which(enrich_database$single_entry == TRUE)
      match_info = match_single(taxon_names = taxon_name[index_to_find_matches],
                                enrich_database = enrich_database,
                                enrich_database_search_index = single_indices,
                                enrich_taxon_name_column = enrich_taxon_name_column,
                                ...)

      taxon_match[index_to_find_matches] = match_info$match
      taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

      index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
      no_found = length(index_to_find_matches[!is.na(match_info$match)])
      cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")

      index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
    }
    }

  ################################################
  # 6) Match collection to all taxon names with a multiple entries.
  ################################################
  if(do_match_multiple){
    if(length(index_to_find_matches) > 0){
      cli::cli_h2("Matching {length(index_to_find_matches)} name{?s} to non-unique taxon names")

      mult_indices = which(enrich_database$single_entry == FALSE)
      match_info = match_multiple(taxon_names = taxon_name[index_to_find_matches],
                                  taxon_authors = taxon_author[index_to_find_matches],
                                  enrich_database = enrich_database,
                                  enrich_database_search_index = mult_indices,
                                  enrich_taxon_name_column = enrich_taxon_name_column,
                                  enrich_taxon_authors_column = enrich_taxon_authors_column,
                                  ...)

      taxon_match[index_to_find_matches] = match_info$match
      taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

      index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
      no_found = length(index_to_find_matches[!is.na(match_info$match)])
      cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")

      index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
    }

  }

  ################################################
  # 7) If the matched author is different after exact matching try to find a similar taxon name.
  # THIS NEEDS CHECKING.
  ################################################
  if(do_taxon_author){
    # Find the records which have disagreeing authors.
    with_match = index_complete[which(!grepl('Not in POWO',taxon_name_story[index_complete]))]

    # Do we have at least one match.
    if(length(with_match) > 0){
      original_authors = taxon_author[with_match]
      matched_authors = rep(NA,length(original_authors))
      match_bigger_0 = taxon_match[with_match] > 0
      matched_authors[match_bigger_0] = enrich_database[[enrich_taxon_authors_column]][taxon_match[with_match][match_bigger_0]]

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
                                     enrich_database = enrich_database,
                                     do_taxon_author = do_taxon_author,
                                     do_add_split = do_add_split,
                                     do_fix_hybrid = do_fix_hybrid,
                                     do_rm_autonym = do_rm_autonym,
                                     matching_criterion = matching_criterion,
                                     enrich_taxon_name_column = enrich_taxon_name_column,
                                     enrich_taxon_authors_column = enrich_taxon_authors_column)

        proposed_authors = rep(NA, length(taxon_author[diff_index]))
        current_match = match_info$match
        current_match[is.na(current_match)] = -3
        proposed_authors[current_match > 0] = enrich_database[[enrich_taxon_authors_column]][current_match[current_match >0]]

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
                   ' -> (Author differ) -> (Try fixing taxon name)',
                   match_info$message[improved_author_index])
        }
      }
    }

  }

  ################################################
  # 8) Try fixing taxon name issues.
  ################################################
  if(do_fix_taxon_name){
    if(length(index_to_find_matches) > 0){
      cli::cli_h2("Testing and matching taxon name issues for {length(index_to_find_matches)} name{?s}")

      match_info = match_all_issue(taxon_names = taxon_name[index_to_find_matches],
                                   taxon_authors = taxon_author[index_to_find_matches],
                                   enrich_database = enrich_database,
                                   do_add_split = do_add_split,
                                   do_fix_hybrid = do_fix_hybrid,
                                   do_rm_autonym = do_rm_autonym,
                                   matching_criterion = matching_criterion,
                                   enrich_taxon_name_column = enrich_taxon_name_column,
                                   enrich_taxon_authors_column = enrich_taxon_authors_column)

      taxon_match[index_to_find_matches] = match_info$match
      taxon_name_story[index_to_find_matches] = paste0(taxon_name_story[index_to_find_matches], match_info$message)

      index_complete = c(index_complete, index_to_find_matches[!is.na(match_info$match)])
      no_found = length(index_to_find_matches[!is.na(match_info$match)])
      cli::cli_alert_success("Found {no_found} of {length(index_to_find_matches)} names")

      index_to_find_matches = index_to_find_matches[is.na(match_info$match)]
    }
  }

  ################################################
  # 9) Try to find typo and then match.
  ################################################
  if(typo_method %in%  c('All', 'Data frame only','Data frame + common') & length(index_to_find_matches) > 0){
    cli::cli_h2("Testing and matching typos for {length(index_to_find_matches)} name{?s}")

    # Check that single_indices and mult_indices exist, if not create them.
    if(!exists('single_indices')){
      single_indices = which(enrich_database$single_entry == TRUE)
    }
    if(!exists('mult_indices')){
      mult_indices = which(enrich_database$single_entry == FALSE)
    }

    match_info = match_typos(taxon_names = taxon_name[index_to_find_matches],
                             taxon_authors = taxon_author[index_to_find_matches],
                             enrich_database = enrich_database,
                             single_indices = single_indices,
                             mult_indices = mult_indices,
                             typo_method = typo_method,
                             enrich_taxon_name_column = enrich_taxon_name_column,
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
  proposed_authors[current_match > 0] = enrich_database[[enrich_taxon_authors_column]][current_match[current_match >0]]
  matched_name = rep(NA, length(taxon_author))
  matched_name[current_match > 0] = enrich_database[[enrich_taxon_name_column]][current_match[current_match >0]]

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


#' Enrich a collection using enrichment database
#'
#' @param collection A data frame containing a collection.
#' @param enrich_database A data frame of enriching information.
#' @param ... Arguments (i.e., attributes) used in the matching algorithm (passed along to nested functions). See [match_collection_to_enrich_database()].
#' @param columns_to_enrich Column names of `enrich_database` that are enriched to `collection`.
#' @param add_to_column_name A string to prepend to new columns in `collection` created by enrichment.
#' @param add_match_details  A flag (TRUE/FALSE) for whether the match details also want to be added to the enriched `collection`.
#'
#' @details
#'
#' This function takes a `collection` and "enriches" it by adding new information taken from `enrich_database` by matching via taxonomic names. The matching is performed by [match_collection_to_enrich_database()], and inputs can be passed to this function via `...`.
#'
#' The enrichment columns to add can be specified by `columns_to_enrich`, and you can include th matching details via `add_match_details`. To make sure the enrichment column names don't clash with names already in the collection each is prepended with `add_to_column_name`.
#'
#' @return `collection` enriched with new columns after matching to `enrich_database`
#' @export
enrich_collection_from_enrich_database <- function(collection, enrich_database,
                                                   ...,
                                                   columns_to_enrich = NA,
                                                   add_to_column_name = 'Enrich_',
                                                   add_match_details = TRUE){

  ### 1) are columns_to_enrich in the enrich_database.
  if(!all(columns_to_enrich %in% names(enrich_database))){
    stop('Not all columns in `columns_to_enrich` are in `enrich_database`!')
  }

  ### 2) Match collection to enrich_database.
  match_info = match_collection_to_enrich_database(collection, enrich_database,...)

  ### 3) Create a new database called enrich_info containing the information to add to
  enrich_info = data.frame(matrix(NA, nrow = nrow(collection), ncol = length(columns_to_enrich)))
  names(enrich_info) = paste0(add_to_column_name, columns_to_enrich)
  indices = which(!(is.na(match_info$match) | match_info$match < 0))
  enrich_info[indices,] = enrich_database[match_info$match[indices],match(columns_to_enrich,names(enrich_database))]

  ### 4) Do we want to add the match details?
  if(add_match_details){
    match_details = data.frame(match_info)
    match_details = match_details[,-1] # remove the match column (i.e row number in enrich database)
    names(match_details) = paste0(add_to_column_name, names(match_details))

    collection = data.frame(collection, match_details)
  }

  ### 5) Add the enrich information and return
  collection = data.frame(collection, enrich_info)
  return(collection)
}
