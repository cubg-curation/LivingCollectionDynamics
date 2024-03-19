#' Matching functions
#'
#' Functions used to match taxonomic names from a collection to exterior databases (POWO's WCVP, IUCN Redlist)
#'
#'  Below we outline the uses of each function. For further details and examples on matching functions please see the `Matching.Rmd` vignette.
#'
#'  Each of the matching functions generally return the index of the matching record in `enrich_database` and a message detailing how the match was obtained. These function can be used as building blocks to build a custom taxonomic name matching algorithm.
#'
#'  - `match_single()` matches `taxon_names` to `enrich_database` taking only the first match. `enrich_database_search_index` should be used to restrict the enrich database to only 'unique' taxonomic names (i.e taxonomic names that correspond to a single record in the enrich database). For 'non-unique' taxonomic names `match_multiple()` should be used.
#'
#' - `match_multiple()` matches `taxon_names` to `enrich_database` for entries in enrich database that have `non-unique` taxonomic names. For 'unique' taxonomic names `match_single()` should be used. For `non-unique` taxonomic names we first use taxonomic author matching to decide which record to use. This matching is performed to each taxonomic name and author using the function `get_match_from_multiple()`. `get_match_from_multiple()` further depends on a  matching criteria function which can be added using the input `matching_criterion` (passed via `...`). By default this is set to [additional_wcvp_matching()], which uses accepted_plant_name_id and taxon_status to chose the best match (in WCVP).
#'
#' - `match_all_issue()` attempts to fix hybridisation, change infraspecific levels or remove autonyms to find matches to an enriched database. This function depends on the functions:
#'
#'    - `try_rm_autonym()` attempts to find taxonomic names in `enrich_database` by removing autonyms.
#'
#'    - `try_fix_infraspecific_level()` attempts to find taxonomic names in `enrich_database` by adding/changing/removing infraspecific levels (var., f., etc).
#'
#'    - `try_fix_hybrid()`  attempts to find taxonomic names in `enrich_database` by adding/changing/removing hybrid markers (+ or x).
#'
#' - `match_typos()` attempts to find matches by searching for typos in the taxonomic name. This depends on the function:
#'
#'    - `check_taxon_typo()` to check a single taxonomic name for typos found either in a typo list or the enriched database.
#' - `no_match_cultivar_indet()` searches for cultivars and indeterminates and sets their match to `-1` indicating no match.
#'
#' - `shorten_message()` compresses matching message (details of how a match is found) into an easy to read format.
#'
#'
#' @param taxon_names Vector of taxonomic names.
#' @param taxon_authors A vector of full taxon names (corresponding to `taxon_names`)
#' @param enrich_database A data frame of enriching information we want to match `taxon_names` to.
#' @param enrich_database_search_index A vector of indices of `enrich_database` that are desired to be matched to.
#' @param single_indices A vector of indices of `enrich_database` that correspond to the records that have 'unique' taxonomic names.
#' @param mult_indices A vector of indices of `enrich_database` that correspond to the records that have 'non-unique' taxonomic names.
#' @param enrich_taxon_name_column The name of the column in `enrich_database` that corresponds to taxonomic names. Default value is `taxon_names`.
#' @param enrich_display_in_message_column The name of the column in `enrich_database` that contains values to show in the matching messages. Default value is `powo_id` (wcvp identifier).
#' @param match_column either `NA` or the name of the column in `enrich_database`. The default value if `NA` which means the values of the match are the indices of the matched records in the enrich database. If instead a single column of `enrich_database` is desired to be the result of the match the name of the column needs to be provided.
#' @param enrich_taxon_authors_column The name of the column in `enrich_database` that corresponds to the authors of taxonomic names. Default value is `taxon_authors_simp`.
#' @param enrich_taxon_author_words_column The name of the column in `enrich_database` that corresponds to the words contained in the authors of taxonomic names. Default value is `author_parts`.
#' @param enrich_database_taxon_names The taxon names taken from `enrich_database`.
#' @param enrich_plant_identifier_column The name of the column in `enrich_database` that corresponds to record identifier. Default value is `plant_name_id`.
#' @param enrich_database_mult `enrich_database` restricted to the rows that correspond to 'non-unique' taxonomic names.
#' @param do_add_split Flag (TRUE/FALSE) for whether we search for missing f./var./subsp.
#' @param do_fix_hybrid Flag (TRUE/FALSE) for whether we search for hybrid issues.
#' @param do_rm_autonym Flag (TRUE/FALSE) for whether we try removing autonyms.
#' @param typo_method Either `'All'`, `'Data frame only'`,`'Data frame + common'`, detailing the level of typo finding required.
#' @param typo_df A data frame where the first column is a taxonomic name with a typo and the second column is the corrected taxonomic name. By default `LivingCollectionDynamics::typo_list` is used.
#' @param taxon_name_and_author the pair of taxonomic name and combined taxonomic name and author
#' @param taxon_name A single taxonomic name.
#' @param show_progress Flag (TRUE/FALSE) for whether we show progress bar.
#' @param ... Arguments (i.e., attributes) used in the matching algorithm (passed along to nested fuctions). Examples include `enrich_taxon_authors_column`, `enrich_display_in_message_column` and `enrich_plant_identifier_column`.
#' @param matching_criterion The function used to find the best match when we have 'non-unique' taxonomic names. By default the function `LivingCollectionDynamics::get_match_from_multiple()` is used.
#' @param matching_authors The function used to find the best match using the author of taxonomic names. By default the function `LivingCollectionDynamics::match_authors()` is used.

#' @param messages messages detailing how a match is obtained.
#' @param console_message Flag (TRUE/FALSE) detailing whether to show messages in the console.
#' @param try_hybrid Flag (TRUE/FALSE) for whether hybrid fixes are attempted.
#'
#' @export
match_single <- function(taxon_names, enrich_database, enrich_database_search_index,
                         enrich_taxon_name_column = 'taxon_name',
                         enrich_display_in_message_column = 'ID',
                         match_column = NA,...){

  # If no indices given return no match.
  if(length(enrich_database_search_index) == 0){
    return(list(match = rep(NA, length(taxon_names)), message = rep('', length(taxon_names))))
  }

  # A) setup
  enriched_taxon_names = enrich_database[[enrich_taxon_name_column]]
  enriched_display_in_message = enrich_database[,match(enrich_display_in_message_column, names(enrich_database))]
  message = rep('', length(taxon_names))

  # B) Perform the matching.
  match_to_single = match(taxon_names, enriched_taxon_names[enrich_database_search_index])

  # C) Find the indices of the match for both taxon names and enrich_database.
  orep_index_match = (1:length(taxon_names))[!is.na(match_to_single)]
  wcvp_index_match = enrich_database_search_index[match_to_single[!is.na(match_to_single)]]

  # D) Update message.
  message[orep_index_match] = paste0(message[orep_index_match], ' -> (matches record with single entry) -> (', enriched_display_in_message[wcvp_index_match],
                                     ', ', enriched_taxon_names[wcvp_index_match],
                                     ')')

  # E) Set the match. If NA return the index in enriched report otherwise select one column to return.
  if(is.na(match_column)){
    matched = enrich_database_search_index[match_to_single]
  }
  else{
    matched = enrich_database[,match(match_column, names(enrich_database))][enrich_database_search_index[match_to_single]]
  }


  return(list(match = matched, message = message))
}

#' @rdname match_single
#' @export
match_multiple <- function(taxon_names,taxon_authors, enrich_database, enrich_database_search_index,
                           enrich_taxon_name_column = 'taxon_name',
                           enrich_display_in_message_column = 'ID',
                           enrich_plant_identifier_column = 'ID',
                           match_column = NA,
                           ...,
                           show_progress = TRUE){

  #################################
  # If no indices given return no match.
  #################################
  if(length(enrich_database_search_index) == 0){
    return(list(match = rep(NA, length(taxon_names)), message = rep('', length(taxon_names))))
  }

  # A) Setup.
  match_to_multiple = rep(NA,length(taxon_names))
  message = rep('',length(taxon_names))
  enriched_plant_identifier = enrich_database[[enrich_plant_identifier_column]]
  enriched_display_in_message = enrich_database[[enrich_display_in_message_column]]
  enriched_taxon_names = enrich_database[[enrich_taxon_name_column]]

  wcvp_multiple = enrich_database[enrich_database_search_index,]
  wcvp_multiple_taxon_name = wcvp_multiple[[enrich_taxon_name_column]]


  # 1) Find which taxon names are in the restricted wcvp.
  in_wcvp = which(taxon_names %in% wcvp_multiple_taxon_name)

  # 2) Names to find matches for. (list of pairs of taxon name and taxon full)
  to_find_match = Map(c, taxon_names[in_wcvp], taxon_authors[in_wcvp])

  # 3) Find the match.
  if(show_progress){
    match_info = pbapply::pblapply(to_find_match, function(x){
      get_match_from_multiple(taxon_name_and_author = x,
                         enrich_database_mult = wcvp_multiple,
                         ...)
      })
  }
  else{
    match_info = lapply(to_find_match, function(x){
      get_match_from_multiple(taxon_name_and_author = x,
                         enrich_database_mult = wcvp_multiple,
                         ...)})
  }

  match_info_plant_identifier = unlist(lapply(match_info,function(x){x[[1]]}))
  match_info_mess = unlist(lapply(match_info,function(x){x[[2]]}))

  # 4) update match_to_multiple and message.
  match_info_match =  match(match_info_plant_identifier, enriched_plant_identifier)
  match_info_match[is.na(match_info_match)] = -2
  match_to_multiple[in_wcvp] = match_info_match

  #message if we agree to a match
  has_accept_match = match_info_match > 0
  message[in_wcvp][has_accept_match] = paste0(message[in_wcvp][has_accept_match], ' -> ', match_info_mess[has_accept_match], ' -> (',
                                              enriched_display_in_message[match_info_match[has_accept_match]], ', ',
                                              enriched_taxon_names[match_info_match[has_accept_match]],
                                              ')')
  #message if we don't agree to a match
  no_accept_match = match_info_match < 0
  message[in_wcvp][no_accept_match] = paste0(message[in_wcvp][no_accept_match], ' -> ', match_info_mess[no_accept_match])

  #Return match and message
  return(list(match = match_to_multiple, message = message))
}

#' @rdname match_single
#' @export
match_all_issue <- function(taxon_names,
                            taxon_authors = rep(NA,length(taxon_names)),
                            enrich_database,
                            matching_authors = LivingCollectionDynamics::match_authors,
                            matching_criterion = LivingCollectionDynamics::additional_wcvp_matching,
                            do_add_split = TRUE,
                            do_fix_hybrid = TRUE,
                            do_rm_autonym = TRUE,
                            enrich_taxon_name_column = 'taxon_name',
                            enrich_taxon_authors_column = 'taxon_authors_simp',
                            enrich_plant_identifier_column = 'ID',
                            enrich_display_in_message_column = 'ID',
                            ...){

  ##############################
  # 1) Setup
  ##############################
  ### 1.1) If there are no taxon_names return NULLs.
  if(length(taxon_names) == 0){
    return(list(match = NULL, message = NULL))
  }

  ### 1.2) If none of the methods are selected return no matches found.
  if(all(c(!do_add_split, !do_fix_hybrid, !do_rm_autonym))){
    return(list(match = rep(NA, length(taxon_names)), message = rep('', length(taxon_names)) ))
  }

  ### 1.3) Get the quantities needed from enrich_database
  enriched_taxon_names = enrich_database[[enrich_taxon_name_column]]

  ### 1.4) Setup whether we need author matching.
  if(all(is.na(taxon_authors))){
    do_taxon_author = FALSE
  }else{
    do_taxon_author = TRUE
  }

  ##############################
  # 2) Find potential fixed taxonomic names.
  ##############################
  ### 2.1) Try removing autonyms.
  if(do_rm_autonym){
    fix_auto = try_rm_autonym(taxon_names = taxon_names,
                              enrich_database_taxon_names = enrich_database[[enrich_taxon_name_column]],
                              ...)
  }else{
    fix_auto = rep('',length(taxon_names))
  }

  ### 2.2) Try changing/removing/adding infraspecific level (var., f., subsp.)
  if(do_add_split){
    fix_splitter = try_fix_infraspecific_level(taxon_names = taxon_names,
                                               enrich_database_taxon_names = enrich_database[[enrich_taxon_name_column]],
                                               ...)
  }else{
    fix_splitter =  rep('',length(taxon_names))
  }

  ### 2.3) Try adding/changing/removing hybrid markers.
  if(do_fix_hybrid){
    fix_hybrid = try_fix_hybrid(taxon_names = taxon_names,
                                enrich_database_taxon_names = enrich_database[[enrich_taxon_name_column]],
                                ...)
  }else{
    fix_hybrid =  rep('',length(taxon_names))

  }

  ### 2.4) Combine fixed names.
  names_to_try = paste0(fix_auto, ' OR ',fix_splitter, ' OR ',fix_hybrid)
  names_to_try = stringr::str_replace_all(names_to_try, pattern = ' OR  OR ', ' OR ') # clean missing values
  names_to_try = stringr::str_replace_all(names_to_try, pattern = '^ OR | OR $', '') # clean missing values

  ##############################
  # 3) Find best match out of fixed names.
  ##############################
  to_try = data.frame(taxon_names = names_to_try, authors = taxon_authors)
  ### 3.1) Loop over all names_to_try.
  counter = 0
  matches = pbapply::pbapply(to_try, 1, function(name_author){
    counter <<- counter +1
    tax_names = name_author[1]
    tax_author = name_author[2]
    current_message = ''

    ### 3.2) if names == '' return no match.
    if(tax_names == ''){
      return(list(match = NA, message = ''))
    }
    current_message = paste0(current_message, ' -> (Try Fixing taxomonic name) -> ',collapse =' ')
    current_message = paste0(current_message, tax_names, ' -> ',collapse =' ')


    ### 3.3) Extract possible matches from enriched database.
    split_names = stringr::str_split(tax_names, pattern = ' OR ')[[1]]
    enriched_cur = enrich_database[enriched_taxon_names %in% split_names,]

    ### 3.4) If only a single record return it.
    if(nrow(enriched_cur) == 1){
      current_message = paste0(current_message, ' (single fixed record) -> ',collapse =' ')
      current_message = paste0(current_message, ' (',enriched_cur[[enrich_display_in_message_column]],
                               ', ', enriched_cur[[enrich_taxon_name_column]],
                               ')',collapse =' ')
      return(list(plant_identifer = enriched_cur[[enrich_plant_identifier_column]], message = current_message))
    }

    enrich_taxon_authors_cur = enriched_cur[[enrich_taxon_authors_column]]
    ### 3.4) Author matching.
    if(do_taxon_author){
      # Get the author matches.
      matched_by_authors = matching_authors(collection_author = tax_author,
                                         enriched_database_authors = enrich_taxon_authors_cur,
                                         ...)

      # Reduce enriched_cur by author matching.
      enriched_cur = enriched_cur[matched_by_authors$wanted,]
      current_message = paste0(current_message, matched_by_authors$message, collapse =' ')

      # If enriched_cur only has one row then we have found the best match and no further matching required.
      if(nrow(enriched_cur) == 1){
        current_message = paste0(current_message, ' -> (',enriched_cur[[enrich_display_in_message_column]],
                                 ', ', enriched_cur[[enrich_taxon_name_column]],
                                 ')',collapse =' ')
        return(list(plant_identifer = enriched_cur[[enrich_plant_identifier_column]], message = current_message))
      }
    }
    else{
      current_message = paste0(current_message, '(No authors) -> ', collapse =' ')
    }

    ### 3.5) Matching criterion. (dependent on enrich_database)
    match_by_criterion = matching_criterion(enrich_database_extract = enriched_cur, message = '')
    enriched_cur = enriched_cur[match_by_criterion$row,]
    current_message = paste0(current_message, match_by_criterion$message, collapse =' ')
    if(nrow(enriched_cur) == 1){
      current_message = paste0(current_message, ' -> (',enriched_cur[[enrich_display_in_message_column]],
                               ', ', enriched_cur[[enrich_taxon_name_column]],
                               ')',collapse =' ')
      return(list(plant_identifer = enriched_cur[[enrich_plant_identifier_column]], message = current_message))
    }

    ### 3.6) Match by method used.(Fix splitter > fix hybrid < remove autonym)
    remaining_taxon_names = enriched_cur[[enrich_taxon_name_column]]
    best_method = unlist(lapply(remaining_taxon_names, function(x){
      if(x %in% stringr::str_split(fix_splitter[counter], pattern = ' OR ')[[1]] ){
        return(1)
      }
      if(x %in% stringr::str_split(fix_hybrid[counter], pattern = ' OR ')[[1]] ){
        return(2)
      }
      if(x %in% stringr::str_split(fix_auto[counter], pattern = ' OR ')[[1]] ){
        return(3)
      }
    }))
    min_best_method = min(best_method,na.rm = T)
    enriched_cur = enriched_cur[which(best_method == min_best_method),]

    if(nrow(enriched_cur) == 1){
      message = '(Decide on Method: Remove Autonym) -> '
      if(min_best_method == 2){message = '(Decide on Method: Fix hybrid) -> '}
      if(min_best_method == 1){message = '(Decide on Method: Fix infraspecific level) -> '}
      current_message = paste0(current_message, message, collapse =' ')
      current_message = paste0(current_message, ' (',enriched_cur[[enrich_display_in_message_column]],
                               ', ', enriched_cur[[enrich_taxon_name_column]],
                               ')',collapse =' ')
      return(list(plant_identifer = enriched_cur[[enrich_plant_identifier_column]], message = current_message))
    }

    current_message = paste0(current_message, '(Cannot decide via fixing method)', collapse =' ')

    # No method can find the single best record match.
    return(list(plant_identifer = -2, message = current_message))

  })

  counter = 0
  m = unlist(lapply(matches, function(x){x[[1]]}))
  mess = unlist(lapply(matches, function(x){x[[2]]}))

  ### Convert enrich_plant_identifier_column to index in enrich database.
  matched = rep(NA, length(taxon_names))
  matched =match(m, enrich_database[[enrich_plant_identifier_column]])
  matched[which(m == -2)] = -2

  return(list(match = matched, message = mess))
}

#' @rdname match_single
#' @export
match_typos <- function(taxon_names, taxon_authors, enrich_database,
                        enrich_taxon_name_column = 'taxon_name',
                        single_indices = NA,
                        mult_indices = NA,
                        typo_method = 'Data frame only', ...){
  enriched_taxon_names = enrich_database[,match(enrich_taxon_name_column, names(enrich_database))]

  #Check for NA in taxon_names and remove if they exist.
  NAs = which(is.na(taxon_names))
  if(length(NAs) > 1){
    warning('In match_typos(), taxon names contain NA.')
  }

  ########################
  # Setup + find typos + create wcvp indices
  ########################
  # setup output
  out_match = rep(NA, length(taxon_names))
  out_message = rep('',length(taxon_names))

  # Get typos
  if(typo_method == 'Data frame only'){
    fixed_typo = unlist(pbapply::pblapply(taxon_names, function(x){check_taxon_typo(x,NA, typo_method = typo_method)}))
  }
  else{
    wcvp_without_repeated = enrich_database[match(unique(enriched_taxon_names),
                                                  enriched_taxon_names),]
    fixed_typo = unlist(pbapply::pblapply(taxon_names, function(x){check_taxon_typo(x,wcvp_without_repeated, typo_method = typo_method)}))

  }

  # Get index of taxons found to have typos.
  typo_indices =  which(!is.na(fixed_typo))

  # If no typos return.
  if(length(typo_indices) == 0){
    return(list(match = out_match, message = out_message))
  }

  if(is.na(single_indices[1])){
    single_indices = which(enrich_database$single_entry == TRUE)
  }
  if(is.na(mult_indices[1])){
    mult_indices = which(enrich_database$single_entry == FALSE)
  }

  name_to_try = fixed_typo[typo_indices]

  ########################
  # Match to single.
  ########################
  match_info = match_single(name_to_try, enrich_database, single_indices, ...)
  out_match[typo_indices] = match_info$match
  out_message[typo_indices] = paste0(out_message[typo_indices], match_info$message)
  index_complete = typo_indices[!is.na(match_info$match)]
  index_to_find_matches = which(is.na(match_info$match))

  ########################
  # Match to multiple.
  ########################

  if(length(index_to_find_matches) > 0 & length(mult_indices) > 0){
    typo_index = typo_indices[index_to_find_matches]
    match_info = match_multiple(name_to_try[index_to_find_matches], taxon_authors[typo_index],  enrich_database, mult_indices, ...)
    out_match[typo_index] = match_info$match
    out_message[typo_index] = paste0(out_message[typo_index], match_info$message)
  }

  ########################
  # Update match message to include the removed autonym.
  ########################
  with_match = out_message[typo_indices] != ''
  out_message[typo_indices[with_match]] = paste0(" -> (Typo) -> ", name_to_try[with_match],  out_message[typo_indices[with_match]])


  return(list(match = out_match, message = out_message))
}

#' @rdname match_single
#' @export
no_match_cultivar_indet <- function(taxon_names){
  #Setup output
  out_match = rep(NA,length(taxon_names))
  out_message = rep('',length(taxon_names))

  #Find indices of those known to be cultivars or indeterminates
  not_in_wcvp = which(grepl(" sp\\.| gx |'.*?'|\\[|\\(|^Indet| gx|indet$|CV|cv$|cv\\.|Group|unkn|hybrid$|Hybrid |Unknown|[0-9]|\\#",taxon_names))

  # Those that end with "grex" with > 2 words.
  words = stringr::str_count(taxon_names, ' ') +1
  with_grex = which(words > 2 & grepl(' grex$| series$', taxon_names))

  # Those that end with a colour. (checked not in wcvp)
  end_in_colour = which(grepl(' blue$| green$| red$| orange$| pink$| yellow$| white$| purple$| black$| brown$| tall$| short$| dwarf$| form$| mix$| mixture$| mixed$|cultivar| group$|selection|indet\\.$| aggr$| aggr\\.$|varient$|varient | sp$|unidentified|yield| hybrids$|strain| male$| female$|variegated|".*?"',taxon_names))

  not_in_wcvp = unique(c(not_in_wcvp, with_grex, end_in_colour))

  # For cultivars and indeterminates set the match to -1 and create message.
  out_match[not_in_wcvp] = -1
  out_message[not_in_wcvp] = ' -> (Cultivar or Indeterminate <Do not attempt matching>)'

  return(list(match = out_match, message = out_message))
}

#' @rdname match_single
#' @export
get_match_from_multiple <- function(taxon_name_and_author, enrich_database_mult,
                                    matching_authors = LivingCollectionDynamics::match_authors,
                                    matching_criterion = LivingCollectionDynamics::no_additional_matching,
                                    enrich_plant_identifier_column = 'plant_name_id',
                                    enrich_taxon_name_column = 'taxon_name',
                                    enrich_taxon_authors_column = 'taxon_authors_simp',
                                    enrich_taxon_author_words_column = 'author_parts',...){


  ##############################
  # 1) Setup
  ##############################
  ### 1.1) Separate taxon names and taxon author, get enriched taxon names.
  enriched_taxon_names = enrich_database_mult[[enrich_taxon_name_column]]
  taxon_name_current = taxon_name_and_author[1]
  taxon_author_current = taxon_name_and_author[2]
  current_message = '(Multiple records in enriched database) -> '

  ### 1.2) Create `try_author_match` flag to check whether we have author information needed for author matching.
  try_author_match = TRUE
  if(is.na(taxon_author_current) || taxon_author_current == ''){
    try_author_match = FALSE
  }

  ### 1.3) Get the corresponding records in enrich_database_mult.
  enriched_cur = enrich_database_mult[enriched_taxon_names == taxon_name_current,]
  enrich_taxon_authors_cur = enriched_cur[[enrich_taxon_authors_column]]
  enrich_taxon_authors_words_cur = enriched_cur[[enrich_taxon_author_words_column]]

  ### 1.4) If there are no corresponding records in enrich_database_mult (nrow(enriched_cur) = 0) then return no match.
  if(nrow(enriched_cur) == 0){
    list(plant_identifer = -2, message = '')
  }


  ##############################
  # 2) Author Matching (reducing enriched_cur to the best matches)
  ##############################
  if(try_author_match){
    # Get the author matches.
    matched_by_authors = matching_authors(collection_author = taxon_author_current,
                  enriched_database_authors = enrich_taxon_authors_cur,
                  ...)

    # Reduce enriched_cur by author matching.
    enriched_cur = enriched_cur[matched_by_authors$wanted,]
    current_message = paste0(current_message, matched_by_authors$message, collapse =' ')

    # If enriched_cur only has one row then we have found the best match and no further matching required.
    if(nrow(enriched_cur) == 1){
      return(list(plant_identifer = enriched_cur[[enrich_plant_identifier_column]], message = current_message))
    }
  }
  else{
    current_message = paste0(current_message, '(No authors)', collapse =' ')
  }

  ##############################
  # 3) Criterion Matching (dependent on enrich_database)
  ##############################
  match_by_criterion = matching_criterion(enriched_cur)
  current_message = paste0(current_message, match_by_criterion$message, collapse =' ')
  # If matching_criterion finds a best match the corresponding row in enriched_cur is returned otherwise the row is set to -2, to denote no match.
  if(length(match_by_criterion$row) == 1){
    enriched_cur = enriched_cur[match_by_criterion$row,]
    return(list(plant_identifer = enriched_cur[[enrich_plant_identifier_column]], message = current_message))
  }

  return(list(plant_identifer = -2, message = current_message))

}

#' @rdname match_single
#' @export
check_taxon_typo <- function(taxon_name, enrich_database = NA,
                             enrich_taxon_name_column = 'taxon_name',
                             typo_df = LivingCollectionDynamics::typo_list,
                             typo_method = 'Data frame only',...){
  ########################
  # 1) Return NA for non-words and special characters
  ########################
  if(is.null(taxon_name))(return(NA))
  if(is.na(taxon_name)){return(NA)}
  if(taxon_name ==''){return(NA)}
  if(grepl('\\(|\\)|\\*|\\?|\\$|\\^|\\[|//]',taxon_name)){return(NA)} # Since no words in wcvp have '(',')', '?' or '*', '$', '^' we return NA.

  ########################
  # 2) Check if the typo is in typo list.
  ########################
  match_to_typo = match(taxon_name,typo_df[,1])
  if(!is.na(match_to_typo)){
    return(typo_df[match_to_typo,2])
  }
  #If we only want to check the typo list return NA for non-matches at this stage.
  if(typo_method == 'Data frame only'){
    return(NA)
  }

  ########################
  # 3) reduce the wcvp names to check. And split into three vectors for same length, one less and one more.
  ########################
  #     A) Make sure the wcvp names are either the same length or one extra character.
  length_taxon = stringr::str_length(taxon_name)
  wcvp_needed = enrich_database[enrich_database$taxon_length %in% c(length_taxon-1, length_taxon, length_taxon+1),]
  #     B) Make sure we only look at taxons which contain only one word with a change (after removing words with '.' or '+')
  words = stringr::str_split(taxon_name,' ')[[1]]
  words = words[!grepl('\\.|\\+|\u00D7', words)]
  pat = paste0(words,collapse = '|')

  wcvp_needed = wcvp_needed[grepl(pat, wcvp_needed[[enrich_taxon_name_column]]),]
  wcvp_needed_minus_1 = wcvp_needed[[enrich_taxon_name_column]][wcvp_needed$taxon_length == (length_taxon-1)]
  wcvp_needed_same = wcvp_needed[[enrich_taxon_name_column]][wcvp_needed$taxon_length == (length_taxon)]
  wcvp_needed_plus_1 = wcvp_needed[[enrich_taxon_name_column]][wcvp_needed$taxon_length == (length_taxon+1)]

  ########################
  # 4) Find common typos in taxon names
  ########################
  # a) final letter change.
  final_letter_change = matrix(c('i','ii',
                      'i', 'ae',
                      'a', 'um',
                      'a', 'us',
                      'a', 'is',
                      'a', 'os',
                      'a', 'er',
                      'er', 'erus',
                      'ae', 'eae',
                      'e', 'is',
                      'ii', 'us',
                      'is','e',
                      'us', 'is',
                      'ense','iense',
                      'oides', 'ioides',
                      'orum', 'iorum'),byrow = T, ncol = 2)
  final_letter_change = rbind(final_letter_change, final_letter_change[,2:1])
  # Function that loops over all final letter changes and returns if typo is found.
  for(i in 1:nrow(final_letter_change)){
    if(stringr::str_ends(taxon_name,final_letter_change[i,1])){
      fixed = wcvp_needed[[enrich_taxon_name_column]][wcvp_needed[[enrich_taxon_name_column]] == stringr::str_replace(taxon_name,paste0(final_letter_change[i,1],'$',collapse = ''),final_letter_change[i,2])]
      if(length(fixed) >0){
        return(fixed[1])
      }
    }
  }

  # b) Common none simple letter swap.
  letter_change = matrix(c('i','ae'),byrow = T, ncol = 2)
  for(i in 1:nrow(letter_change)){
    locations = stringr::str_locate_all(taxon_name, letter_change)
    for(j in 1:length(locations)){
      current = locations[[j]]
      if(nrow(current) > 0){
        for(k in 1:nrow(current)){
          # original letter in middle.
          if(current[k,1] == 1){
            new_name = paste0(letter_change[i,3-j],
                              stringr::str_sub(taxon_name,current[k,2]+1,-1))
          }
          #original letter at end.
          else if(current[k,2] == length(taxon_name)){
            new_name = paste0(stringr::str_sub(taxon_name,1,current[k,1]-1),
                              letter_change[i,3-j])
          }
          #original letter in middle
          else{
            new_name = paste0(stringr::str_sub(taxon_name,1,current[k,1]-1),
                              letter_change[i,3-j],
                              stringr::str_sub(taxon_name,current[k,2]+1,-1))
          }

          fixed = wcvp_needed[[enrich_taxon_name_column]][wcvp_needed[[enrich_taxon_name_column]] == new_name]
          if(length(fixed) >0){
            return(fixed[1])
          }
        }
      }
    }
  }

  #If we only want to check the typo list return NA for non-matches at this stage.
  if(typo_method == 'Data frame + Common'){
    return(NA)
  }
  ########################
  # 4) More general search of one letter differences.
  ########################
  for(i in (length_taxon-1):1){
    #Check changing a single letter.
    patternA = paste0(stringr::str_sub(taxon_name,1,i),'[a-zA-Z]',stringr::str_sub(taxon_name,i+2,length_taxon))
    patternA = stringr::str_replace_all(patternA,'\\.','\\\\.')
    fixed_typoA = wcvp_needed_same[grepl(patternA, wcvp_needed_same)]
    if(length(fixed_typoA) >0){
      return(fixed_typoA[1])
    }

    #Check adding a single new letter
    patternB = paste0(stringr::str_sub(taxon_name,1,i),'[a-zA-Z-]',stringr::str_sub(taxon_name,i+1,length_taxon))
    patternB = stringr::str_replace_all(patternB,'\\.','\\\\.')
    fixed_typoB = wcvp_needed_plus_1[grepl(patternB, wcvp_needed_plus_1)]
    if(length(fixed_typoB) >0){
      return(fixed_typoB[1])
    }
    #Check removing a single new letter
    patternC = paste0(stringr::str_sub(taxon_name,1,i),stringr::str_sub(taxon_name,i+2,length_taxon))
    patternC = stringr::str_replace_all(patternC,'\\.','\\\\.')
    fixed_typoC = wcvp_needed_minus_1[grepl(patternC, wcvp_needed_minus_1)]
    if(length(fixed_typoC) >0){
      return(fixed_typoC[1])
    }
  }

  # If no typo is found by this stage we haven't managed to find a fix and return NA.
  return(NA)
}

#' @rdname match_single
#' @export
shorten_message <- function(messages){
  ### 1) Create Match options with phrases found in match_detail and how we want to shorten the message.
  match_options = c("(Sanitise)", 'SANITISE',
                    'Try Fixing taxomonic name', 'FIX',
                    "(Typo)", 'TYPO',
                    "Multiple records in enriched database", 'MULTIPLE',
                    "(matches record with single entry)", 'SINGLE',
                    "Author differ", "AUTHOR_DIFF",
                    "(Exact author match)", 'AUTHOR',
                    "all point to same accepted plant", 'SAME_ACC',
                    "Partial author", 'PARTIAL',
                    "choose via taxon_status", "TAXON_STATUS",
                    # "multiple best taxon status, do not match", "UNCLEAR",
                    # "no accepted or synonym", "UNCLEAR",
                    "Decide on Method", 'METHOD',
                    "(Cultivar or Indeterminate <Do not attempt matching>)", 'CULT/INDET',
                    "Remove autonym", 'AUTONYM',
                    "(No match found)", 'NO_MATCH',
                    "single fixed record", 'SINGLE',
                    "(Go to accepted name)", 'ACCEPTED',
                    "Infrageneric level update", 'INFRA',
                    "Hybrid fix", 'HYBRID'
  )
  match_options = stringr::str_replace_all(match_options, '\\(', '\\\\(')
  match_options = stringr::str_replace_all(match_options, '\\)', '\\\\)')
  match_options = matrix(match_options, byrow = T, ncol=2)

  ### 2) For each phrase location the position in the message.
  locations = matrix(NA, nrow = length(messages), ncol = nrow(match_options))
  for(i in 1:nrow(match_options)){
    loc = stringr::str_locate(messages, match_options[i,1])[,1]
    locations[,i] = loc
  }

  ### 3) Use locations to decide the order of the shortened message.
  match_short = apply(locations, 1, function(x){
    has_value_index = which(!is.na(x))
    order_index = has_value_index[order(x[has_value_index])]
    paste0(match_options[order_index,2], collapse = ', ')
  })

  #If we have multiple of EXACT, PARTIAL, TAXON_STATUS only keep the worst level of matching.
  match_short = stringr::str_replace_all(match_short,'EXACT, PARTIAL','PARTIAL')
  match_short = stringr::str_replace_all(match_short,'PARTIAL, TAXON_STATUS','TAXON_STATUS')
  match_short = stringr::str_replace_all(match_short,'AUTHOR, TAXON_STATUS','TAXON_STATUS')
  match_short = stringr::str_replace_all(match_short,'PARTIAL, SAME_ACC','SAME_ACC')
  match_short = stringr::str_replace_all(match_short,'AUTHOR, SAME_ACC','SAME_ACC')
  match_short = stringr::str_replace_all(match_short,'EXACT, UNCLEAR','UNCLEAR')
  match_short = stringr::str_replace_all(match_short,'PARTIAL, UNCLEAR','UNCLEAR')
  match_short = stringr::str_replace_all(match_short,'AUTHOR, METHOD','METHOD')
  match_short = stringr::str_replace_all(match_short,'PARTIAL, METHOD','METHOD')

  # Find the unclear matches, these are when we don't have NO_MATCH or CULT/INDET but we do not point to a matched record.
  possible_index = which(!match_short %in% c('NO_MATCH', 'CULT/INDET'))
  no_code = !grepl('[0-9][0-9\\-]+',messages[possible_index])
  match_short[possible_index[no_code]] = paste0(match_short[possible_index[no_code]], ', ', 'UNCLEAR')
  return(match_short)
}



#' @rdname match_single
#' @export
try_rm_autonym <- function(taxon_names, enrich_database_taxon_names,
                           console_message = TRUE, ...
){
  ##############################
  # 1) Setup
  ##############################
  ### 1.1) Check for NA in taxon names and give warning if there are.
  if(any(is.na(taxon_names))){
    warning('In try_rm_autonym(), taxon names contain NA.')
  }

  ### 1.2) If no taxon_names provided return NULLs.
  if(length(taxon_names) == 0){
    return(list(match = NULL, message = NULL))
  }

  if(console_message){
    cli::cli_alert_info(text = "Trying removing autonyms from taxon names")
  }

  ##############################
  # 2) Find which taxon_names are autonyms and create vector of taxon names with autonym removed that are in enrich database.
  ##############################
  ### 2.1) As is_autonym column to the taxon_names.
  taxon_names_and_autonym = add_is_autonym(data.frame(TaxonName = taxon_names))
  autonym_indices =  which(taxon_names_and_autonym$is_autonym)

  ### 2.2) If there exist no_autonyms exist function.
  if(length(autonym_indices) == 0){
    return(rep('',length(taxon_names)))
  }

  ### 2.3) extract the base of the autonyms and save in autonym_base.
  autonym_base = rep(NA,length(taxon_names))
  autonym_base[autonym_indices] = stringr::str_extract(string = taxon_names[autonym_indices],
                                                       pattern = '^.*?(?= var\\. | subsp\\. | f\\. | ssp\\. | nothosubsp\\. )')


  ### 2.4) Reduce autonym base to NA, if the base is not found in enrich_database.
  autonym_base[!autonym_base %in% enrich_database_taxon_names] = ''

  return(autonym_base)
}

#' @rdname match_single
#' @export
try_fix_infraspecific_level <- function(taxon_names, enrich_database_taxon_names,
                                        try_hybrid = TRUE,
                                        console_message = TRUE,
                                        ...){
  ##############################
  # 1) Setup
  ##############################
  ### 1.1) Check for NA in taxon names and give warning if there are.
  if(any(is.na(taxon_names))){
    warning('In try_rm_autonym(), taxon names contain NA.')
  }

  ### 1.2) If no taxon_names provided return NULLs.
  if(length(taxon_names) == 0){
    return(list(match = NULL, message = NULL))
  }

  ### 1.3) Define the infraspecific levels.
  splitters = c('subsp.', 'var.', 'f.', 'nothosubsp.')
  splitters_grepl = ' subsp\\. | var\\. | f\\. | nothosubsp\\. '

  ### 1.4) Find how many word each taxon_name has.
  no_words = stringr::str_count(taxon_names, ' ')+1

  ### 1.5) Define output variable
  out_fixed_names = rep('',length(taxon_names))

  ### 1.6) Set lapply to pbapply::pblapply to show console progression.
  if(console_message){
    lapply = pbapply::pblapply
  }

  ### 1.7) Get reduced enriched taxon names that contain infraspecific levels.
  enrich_taxon_names_w_split = enrich_database_taxon_names[grepl(splitters_grepl,enrich_database_taxon_names)]


  ### 1.8) If there are no enriched taxon names with infraspecific levels return emply out_fixed_names.
  if(length(enrich_taxon_names_w_split) == 0){
    return(out_fixed_names)
  }

  ##############################
  # 2) Try adding a infraspecific level. (taxon name needs 3 words)
  ##############################
  ### 2.1) Get the indices of taxon_names with 3 words
  index_words_3 = which(no_words == 3 & !grepl('\u00D7|\\+',taxon_names))

  if(length(index_words_3)>0){
    if(console_message){
      cli::cli_alert_info("Trying adding infraspecific level to {length(index_words_3)} name{?s}")
    }
    ### 2.2) Get the taxon names to try with adding infraspecific (/hybrid) markers. whilst checking if in enrich_database.
    words_3 = stringr::str_split(taxon_names[index_words_3], ' ')
    if(try_hybrid){
      new_taxon_names = unlist(lapply(words_3, function(x){
        just_splitter = paste(x[1], x[2], splitters, x[3])
        taxon_namesA = just_splitter[just_splitter %in% enrich_taxon_names_w_split]
        hybrid_and_splitter =  paste(x[1], '\u00D7', x[2], splitters, x[3])
        taxon_namesB = hybrid_and_splitter[hybrid_and_splitter %in% enrich_taxon_names_w_split]
        new_taxon_names = c(taxon_namesA,taxon_namesB)
        return(paste0(new_taxon_names, collapse = ' OR '))
      }))

    }
    else{
      new_taxon_names = unlist(lapply(words_3, function(x){
        just_splitter = paste(x[1], x[2], splitters, x[3])
        new_taxon_names = just_splitter[just_splitter %in% enrich_taxon_names_w_split]

        return(paste0(new_taxon_names, collapse = ' OR '))
      }))

    }
    out_fixed_names[index_words_3] = paste0(out_fixed_names[index_words_3], new_taxon_names, sep = '')
  }


  ##############################
  # 3) Try adding a infraspecific level when hybrid marker exists. (taxon name needs 4 words)
  ##############################
  if(try_hybrid){
    ### 3.1) Get the indices of taxon_names with 4 words and a single hybrid marker in the best position.
    index_words_4 = which(no_words == 4 & stringr::str_count(taxon_names, '\u00D7|\\+') == 1)
    hybrid_position = unlist(lapply(stringr::str_split(taxon_names[index_words_4], ' '),function(x){which(grepl('\u00D7|\\+',x))}))
    index_words_4 = index_words_4[which(hybrid_position %in% c(1,2))]

    if(length(index_words_4)>0){
      if(console_message){
        cli::cli_alert_info("Trying adding infraspecific level (taxon with hybrid markers) to {length(index_words_4)} name{?s}")
      }
      ### 3.2) Get the taxon names to try with adding infraspecific (/hybrid) markers. whilst checking if in enrich_database.
      words_4 = stringr::str_split(taxon_names[index_words_4], ' ')
      new_taxon_names = unlist(lapply(words_4, function(x){
        just_splitter = paste(x[1], x[2], x[3], splitters, x[4])
        taxon_names_with_split = just_splitter[just_splitter %in% enrich_taxon_names_w_split]

        return(paste0(taxon_names_with_split, collapse = ' OR '))
      }))
      out_fixed_names[index_words_4] = paste0(out_fixed_names[index_words_4], new_taxon_names, sep = '')
    }


  }

  ##############################
  # 4) Try changing an infraspecific level.
  ##############################
  ### 4.1) Get the indices of taxon_names with 3 words
  index_words_with_splitter = which(grepl(splitters_grepl, taxon_names))

  if(length(index_words_with_splitter) >0){
    if(console_message){
      cli::cli_alert_info("Trying changing infraspecific level to {length(index_words_with_splitter)} name{?s}")
    }
    new_taxon_names = unlist(lapply(taxon_names[index_words_with_splitter], function(x){
      A=stringr::str_replace(x,pattern = splitters_grepl, ' var\\. ')
      B=stringr::str_replace(x,pattern = splitters_grepl, ' f\\. ')
      C=stringr::str_replace(x,pattern = splitters_grepl, ' subsp\\. ')
      D=stringr::str_replace(x,pattern = splitters_grepl, ' nothosubsp\\. ')
      options = c(A,B,C,D)
      new_taxon_names = options[-match(x, options)]

      new_taxon_names = new_taxon_names[new_taxon_names %in% enrich_taxon_names_w_split]

      return(paste0(new_taxon_names, collapse = ' OR '))
    }))
    out_fixed_names[index_words_with_splitter] = paste0(out_fixed_names[index_words_with_splitter], new_taxon_names, sep = '')
  }

  return(out_fixed_names)
}

#' @rdname match_single
#' @export
try_fix_hybrid <- function(taxon_names, enrich_database_taxon_names,
                           try_hybrid = TRUE,
                           console_message = TRUE,
                           ...){
  ##############################
  # 1) Setup
  ##############################
  ### 1.1) Check for NA in taxon names and give warning if there are.
  if(any(is.na(taxon_names))){
    warning('In try_rm_autonym(), taxon names contain NA.')
  }

  ### 1.2) If no taxon_names provided return NULLs.
  if(length(taxon_names) == 0){
    return(NULL)
  }

  ### 1.2) If we don't want to fix hybrid things.
  if(!try_hybrid){
    return(rep('',length(taxon_names)))
  }

  ### 1.3) Set lapply to pbapply::pblapply to show console progression.
  if(console_message){
    lapply = pbapply::pblapply
  }

  ### 1.4) Define the hybrid markers
  hybrids = c('\u00D7', '+')
  hybrids_grepl = ' \u00D7 | \\+  '

  ### 1.5) Find how many word each taxon_name has.
  no_words = stringr::str_count(taxon_names, ' ')+1

  ### 1.6) Define output variable
  out_fixed_names = rep('',length(taxon_names))

  ### 1.7) Get reduced enriched taxon names that contain hybrids.
  enrich_taxon_names_w_hybrid = enrich_database_taxon_names[grepl(hybrids_grepl,enrich_database_taxon_names)]

  ##############################
  # 2) Try hybrid at start. (taxon name needs 1 word)
  ##############################
  ### 2.1) Get the indices of taxon_names with 1 word
  index_words_1 = which(no_words == 1)

  if(length(index_words_1)>0){
    if(console_message){
      cli::cli_alert_info("Trying fixing hybrid for taxon names with 1 words {length(index_words_1)} name{?s}")
    }
    ### 2.2) Get the taxon names By adding hybrid before first word.
    words_1 = taxon_names[index_words_1]
    new_taxon_names = unlist(lapply(words_1, function(x){
      new_taxon_names = c(paste('+',x,collapse =' '), paste('\u00D7', x, collapse = ' '))
      new_taxon_names = new_taxon_names[new_taxon_names %in% enrich_taxon_names_w_hybrid]

      return(paste0(new_taxon_names, collapse = ' OR '))
    }))

    out_fixed_names[index_words_1] = paste0(out_fixed_names[index_words_1], new_taxon_names, sep = '')

  }

  ########################
  # 3) 2/3/4 words try hybrid at start or after first word.
  ########################
  ### 3.1) Get the indices of taxon_names with 2/3/4 words
  index_words_2_3_4 = which(no_words %in% c(2,3,4), !grepl('\u00D7|\\+',taxon_names))

  if(length(index_words_2_3_4) > 0){
    if(console_message){
      cli::cli_alert_info("Trying fixing hybrid for taxon names with 2/3/4 words {length(index_words_2_3_4)} name{?s}")
    }

    words_2_3_4 = stringr::str_split(taxon_names[index_words_2_3_4], ' ', n=2)
    new_taxon_names = unlist(lapply(words_2_3_4, function(x){
      ### 3.2) Get the taxon names to try with adding infraspecific (/hybrid) markers. whilst checking if in enrich_database.
      before_first_word = paste(hybrids, x[1],x[2])
      after_first_word = paste(x[1], hybrids, x[2])
      new_taxon_names = c(before_first_word,after_first_word)
      new_taxon_names = new_taxon_names[new_taxon_names %in% enrich_taxon_names_w_hybrid]

      return(paste0(new_taxon_names, collapse = ' OR '))
    }))

    out_fixed_names[index_words_2_3_4] = paste0(out_fixed_names[index_words_2_3_4], new_taxon_names, sep = '')


  }

  ########################
  # 4) Change/remove hybrid (one hybrid marker in the taxon name)
  ########################
  ### 3.1) Get the indices of taxon_names with hybrid characters.
  index_words_with_hybrid = which(stringr::str_count(taxon_names, '\u00D7|\\+') == 1)

  if(length(index_words_with_hybrid) > 0){
    if(console_message){
      cli::cli_alert_info("Trying changing/removing hybrid for taxon names {length(index_words_with_hybrid)} name{?s}")
    }

    # Get the words with the hybrid changed or removed.
    new_taxon_names = unlist(lapply(taxon_names[index_words_with_hybrid], function(x){
      A=stringr::str_replace(x,pattern = '\u00D7|\\+', '')
      B=stringr::str_replace(x,pattern = '\u00D7|\\+', '\\+')
      C=stringr::str_replace(x,pattern = '\u00D7|\\+', '\u00D7')
      options = c(A,B,C)
      options = stringr::str_squish(options)
      new_taxon_names = options[-match(x, options)]
      new_taxon_names = new_taxon_names[new_taxon_names %in% enrich_database_taxon_names]

      return(paste0(new_taxon_names, collapse = ' OR '))
    }))

    out_fixed_names[index_words_with_hybrid] = paste0(out_fixed_names[index_words_with_hybrid], new_taxon_names, sep = '')

  }


  return(out_fixed_names)

}

#' Find best author matches
#'
#' @param collection_author The author from the collection wanted to be matched
#' @param enriched_database_authors Author options from the enrich_database.
#' @param partial_method Either `'most words` or `'any words'`, defining the method used to find partial matches.
#' @param ... Arguments (i.e., attributes) used in the matching algorithm (passed along to nested fuctions).
#'
#' @return a list of length 2 with:
#'  - `$wanted` is a logical (TRUE/FALSE) vector with length `length(enriched_database_authors)` corresponding to the enriched database authors that most match the collection author.
#'  - `$message` detailing whether the author match was exact, partial or no match was found.
#' @export
#'
#' @examples
#' collection_author = 'Schult'
#' enriched_database_authors = c("(Lour.) Schult", "Borhidi & E.Martinez")
#' match_authors(collection_author, enriched_database_authors)
match_authors <- function(collection_author, enriched_database_authors, partial_method = 'most words', ...){

  ### 1) Exact matching
  exact_match = enriched_database_authors  == collection_author
  exact_match[is.na(exact_match)] = FALSE # Set NA values to FALSE
  if(any(exact_match)){
    return(list(wanted = exact_match, message = '(Exact author match)'))
  }

  ### 2) Partial matching
  if(partial_method == 'most words'){
    # Get the words for the collections and databases authors.
    enriched_database_authors_words = lapply(enriched_database_authors, author_words)
    collection_author_words = author_words(collection_author)

    #Find number of words in database's authors found in the collection's author
    no_database_author_in_collection = unlist(lapply(enriched_database_authors_words,function(words){
      words = words[words != '']
      contain_words = unlist(lapply(words, function(x){grepl(x,collection_author)}))
      return(sum(contain_words))
    }))

    #Find number of words in collection's author found in the database's authors
    no_collection_author_in_database = rowSums(data.frame(lapply(collection_author_words, function(x){grepl(x,enriched_database_authors)})))

    # Combine above and find the maximum shared words.
    total_match_word_count = rowSums(cbind(no_database_author_in_collection,no_collection_author_in_database))
    max_words_found =  max(total_match_word_count, na.rm=T)

    # If maximum shared words > 0 return those with the maximum number of shared words.
    if(max_words_found > 0){
      match_author_words = total_match_word_count == max_words_found
      match_author_words[is.na(match_author_words)] = FALSE # Set NA values to FALSE
      return(list(wanted = match_author_words, message = '(Partial author <most words>)'))
    }
  }
  else if(partial_method == 'any words'){
    author_checks = unlist(lapply(enriched_database_authors,function(x){author_check(collection_author,x)}))
    if(any(author_checks == 'Partial')){
      partial_authors = author_checks == 'Partial'
      partial_authors[is.na(partial_authors)] = FALSE # Set NA values to FALSE
      return(list(wanted = partial_authors, message = '(Partial author <any words>)'))
    }
  }else{
    stop('Invalid partial_method input in match_authors()!')
  }


  ### 3) No matching
  return(list(wanted = rep(TRUE, length(enriched_database_authors)), message = '(No authors match)'))
}

