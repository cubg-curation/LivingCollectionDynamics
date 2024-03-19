## These scripts are for obtaining the information from POWO (WCVP) used for enrichment.

#' Get accepted plant details from POWO online
#'
#'This function queries POWO online for the accepted name of a plant and the corresponding Plants of the World Online identifier.
#'
#'If the connection to POWO doesn't work (or exist) then the function returns (NA, NA).
#'
#' @param powo_id Plants of the World Online identifier.
#'
#' @return A vector of length two containing the accepted name and accepted powo_id.
#' @export
#'
#' @examples
#' get_accepted_plant('582307-1')
#'
get_accepted_plant <- function(powo_id) {
  out <- tryCatch(
    {
      powo_info = taxize::pow_lookup(paste0('urn:lsid:ipni.org:names:',powo_id))
      new_accepted_name = powo_info$meta$accepted$name
      new_accepted_id = stringr::str_remove(powo_info$meta$accepted$fqId, 'urn:lsid:ipni.org:names:')
      return(c(new_accepted_name, new_accepted_id))

    },
    error=function(cond) {

      return(c(NA,NA))
    },
    warning=function(cond) {
      message(paste("URL caused a warning for powo:", powo_id))
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    }



  )
}


#' Information from wcvp_names
#'
#' We allow four three for uploading the information from wcvp_names:
#' - By choosing a file path to a csv file containing the information using `filepath`.
#' - By loading an R save (file type `.RData`, `.rdata`, `.rda` or `.rds`) via `filepath`.
#' - By using the package  `rWCVPData`, by setting `use_rWCVPData` = TRUE.
#'
#' The columns extracted always include:"plant_name_id", "taxon_name", "taxon_authors", "taxon_rank", "accepted_plant_name_id" and "powo_id". Further columns can be extracted using `wanted_columns` parameter, where either columns names or indices can be used.
#'
#' We add two new columns:
#' - 'taxon_length' which contains the number of characters in a given taxon_name used when searching for typos in garden databases.
#' - 'is_autonym' contains a logical flag for whether the taxon name is an autonym. This is used to match autonyms in garden databases.
#'
#' Moreover, the code performs fixes to the data provided namely:
#' - Checks for '[*]' and '[**]' in taxon_name and taxon_rank, which should be 'var.' and 'f.' instead.
#' - There are potentially plants that are synonyms without an accepted form in wcvp_names. However, in some cases the accepted form can be found on POWO online. Therefore, in this case we scrap POWO online for an accepted form and add to the loaded information in wcvp_names.
#' - No author is provided for autonyms, so we extract the corresponding author from matching GenusSpecies record.
#'
#' @param filepath path to wcvp_names.csv file.
#' @param use_rWCVPdata Flag for whether we use rWCVPdata to get wcvp_names
#' @param wanted_columns specify extra columns to extract from `wcvp_names.csv`.
#'
#' @return Data frame of desired information from wcvp_names
#' @export
#'
import_wcvp_names <- function(filepath=NULL, use_rWCVPdata = FALSE, wanted_columns = c("plant_name_id", "taxon_rank", "taxon_status", "family", "genus", "species", "lifeform_description", "climate_description", "taxon_name", "taxon_authors", "accepted_plant_name_id", "powo_id")){
  cli::cli_h1("Importing wcvp information (wcvp_names)")

  #################################
  # 1) Load wcvp_names (using method as described in the download)
  #################################
  cli::cli_h2("Loading wvcp_names")
  if(!is.null(filepath)){
    # We have a filepath, try and load depending on file format.
    if(grepl('\\.csv$',filepath)){
      wcvp_names = utils::read.table(filepath, sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

    }
    else if(grepl('\\.RData$|\\.rdata$|\\.rda$',filepath)){
      load(filepath)
      wcvp_names = data.frame(wcvp_names) # To reformat to standard data frame not a tibble.
    }
    else if(grepl('\\.rds$',filepath)){
      readRDS(filepath)
      wcvp_names = data.frame(wcvp_names) # To reformat to standard data frame not a tibble.
    }
    else{
      stop('filepath is not a csv or R data file!')
    }

    } else if(use_rWCVPdata){
    if(system.file(package='rWCVPdata') == ''){
      stop('Package rWCVPdata is not installed please use another method to load wcvp_names or install the package')
    }
      if(requireNamespace('rWCVPdata', quietly = T)){
        wcvp_names <- rWCVPdata::wcvp_names

      }

  }

  if(!exists('wcvp_names')){
    stop('Have not successfully loaded wcvp_names, please make sure you provide either: filepath, filepath or set use_rWCVPdata = TRUE')
  }

  columns  = match(c("plant_name_id", "taxon_name", "taxon_authors", "taxon_rank", "taxon_status", "accepted_plant_name_id", "powo_id"), names(wcvp_names))
  if(is.character(wanted_columns)){
    wanted_columns = match(wanted_columns, names(wcvp_names))
  }
  columns = unique(c(columns,wanted_columns))
  wcvp_names = wcvp_names[,columns]

  # Change the taxon name for one special case.
  wcvp_names$taxon_name[wcvp_names$taxon_name == 'xx viridissimus var. viridissimus'] = 'Trigonostemon viridissimus var. viridissimus'


  #################################
  # 2) Prepare wcvp_names.
  #################################
  wcvp_names <- prepare_enrich_database(wcvp_names,
                          enrich_taxon_name_column = 'taxon_name',
                          enrich_taxon_authors_column = 'taxon_authors',
                          do_sanitise = TRUE,
                          do_taxon_length = TRUE,
                          do_single_entry = TRUE,
                          do_author_parts = TRUE,
                          do_add_id = FALSE,
                          do_sort = FALSE,
                          console_message = TRUE)
  indices_require_sanitise = which(wcvp_names$require_sanitise)

  #################################
  # 2) Sanitise taxon names.
  #################################
  # cli::cli_h2("Sanitising taxon names")
  # santise_taxon_name = unlist(pbapply::pblapply(wcvp_names$taxon_name,LivingCollectionDynamics::sanitise_name))
  # original_taxon_name = wcvp_names$taxon_name
  # wcvp_names$taxon_name = santise_taxon_name
  # indices_require_sanitise=which(santise_taxon_name != original_taxon_name)
  # cli::cli_alert_success("Sanitising required for {length(indices_require_sanitise)} taxon names")

  #################################
  # 3) Update author of autonyms to match the none-autonymed species.
  #################################
  # If no author is provided for an autonym, extract the Genus Species then:
  # - Use the author of genus species if only a single record exists.
  # - IF > 1 record exists use the accepted record, otherwise return NA.
  cli::cli_h2("Updating author for autonyms...")

  # Create a data frame of taxon name, is_autonym, authors, and the GenusSpecies.
  GenusSpecies = paste0(wcvp_names$genus, ' ', wcvp_names$species)
  autonym_names = LivingCollectionDynamics::add_is_autonym(data.frame(taxon = wcvp_names$taxon_name),taxon_name_column = 'taxon',progress_bar = T)
  autonym_names = data.frame(autonym_names,
                             author = wcvp_names$taxon_authors,
                             GenusSpecies = GenusSpecies)

  # Those we want to try and update.
  issue_index = which(autonym_names$is_autonym == TRUE & is.na(autonym_names$author))
  cli::cli_alert_success("Updating required for {length(issue_index)} records")

  new_author = unlist(pbapply::pblapply(issue_index, function(x){
    genusSpecies = autonym_names$GenusSpecies[x]
    POWO_records = wcvp_names[wcvp_names$taxon_name == genusSpecies,c(3,5)]
    if(nrow(POWO_records) == 1){
      return(POWO_records$taxon_authors)
    }
    else{
      accepted_record = which(POWO_records$taxon_status == 'Accepted')
      if(!is.na(accepted_record[1])){
        return(POWO_records$taxon_authors[accepted_record])
      }
      else{return(NA)}
    }
  }))
  wcvp_names$taxon_authors[issue_index] = new_author

  #################################
  # 3) Sanitise authors.
  #################################
  # cli::cli_h2("(2/6) Creating author parts...")
  # authors = wcvp_names$taxon_authors
  # authors = stringi::stri_trans_general(authors, id = "Latin-ASCII")
  # remove_initials = unlist(pbapply::pblapply(authors,function(x){paste0(author_words(x),collapse =', ')}))
  # wcvp_names$taxon_authors_simp = authors
  # wcvp_names$author_parts = remove_initials

  #################################
  # 4) Create new column for the length of the taxonName (used if we search for typos
  #################################
  # cli::cli_h2("(3/6) Adding taxon_length column")
  # taxon_length = unlist(pbapply::pblapply(wcvp_names$taxon_name, stringr::str_length))
  # wcvp_names$taxon_length = taxon_length

  #################################
  # 5) Add flag for whether there exists multiple taxon names.
  #################################
  # cli::cli_h2("(4/6) Splitting wcvp_names into those with multiple Taxon names and those with only one...")
  # name_freq = table(wcvp_names$taxon_name)
  # single_entry = rep(NA, nrow(wcvp_names))
  # single_entry[wcvp_names$taxon_name %in% names(name_freq)[as.numeric(name_freq) == 1]] = TRUE
  # single_entry[wcvp_names$taxon_name %in% names(name_freq)[as.numeric(name_freq) > 1]] = FALSE
  # wcvp_names = data.frame(wcvp_names, single_entry)

  #################################
  # 6) Where possible update records that are synonyms without an accepted form to include an accepted form via powo.
  #################################
  cli::cli_h2("Checking for accepted form issues with synonyms...")
  #find the indices of the entries that have missing accepted plant name id when their taxon status is a synonym.
  issue_index = which(wcvp_names$taxon_status == 'Synonym' & is.na(wcvp_names$accepted_plant_name_id))
  non_accepted_synonyms = wcvp_names[issue_index,]

  powo_ids = non_accepted_synonyms$powo_id
  cli::cli_alert_danger("Found {length(powo_ids)}  synonyms without accepted form...")
  cli::cli_text("searching powo online for accepted forms...")
  pbapply::pboptions(type = "timer")
  accepted_details = pbapply::pblapply(powo_ids,get_accepted_plant)
  new_accepted_name = unlist(lapply(accepted_details, function(x){x[1]}))
  new_accepted_powo_id = unlist(lapply(accepted_details, function(x){x[2]}))
  new_accepted_plant_name_id = wcvp_names$accepted_plant_name_id[match(new_accepted_powo_id, wcvp_names$powo_id)]

  wcvp_names$accepted_plant_name_id[issue_index] = new_accepted_plant_name_id

  remaining_NA = sum(is.na(new_accepted_plant_name_id))
  fixed_synonyms = 100 - remaining_NA/length(new_accepted_plant_name_id)*100
  cli::cli_alert_success("Fixed {round(fixed_synonyms,digits=1)}% of missing accepted forms...")
  cli::cli_alert_danger("{remaining_NA} remaining unexplained missing accepted forms......")

  #################################
  # 7)  Check for common hybrids, cultivars in wcvp to check we can definitely exclude these initally. Find exceptions that are included in wcvp_names.
  #################################
  cli::cli_h2("Checking for common hybrids, cultivar, etc symbols...")

  # A)  Note current wcvp has a bug where [*] and [**] are sometimes used instead if var. and f.
  # Fix this here
  with_sqbracket = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "\\[") >= 1,]
  with_sqbracket_original = with_sqbracket
  with_sqbracket$taxon_name = stringr::str_replace(with_sqbracket$taxon_name,'\\[\\*\\*\\]', 'f\\.')
  with_sqbracket$taxon_name = stringr::str_replace(with_sqbracket$taxon_name,'\\[\\*\\]', 'var\\.')
  with_sqbracket$taxon_rank = stringr::str_replace(with_sqbracket$taxon_rank,'\\[\\*\\*\\]', 'Form')
  with_sqbracket$taxon_rank = stringr::str_replace(with_sqbracket$taxon_rank,'\\[\\*\\]', 'Variety')
  wcvp_names[stringr::str_count(wcvp_names$taxon_name, "\\[") >= 1,] = with_sqbracket
  with_sqbracket2 = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "\\[") >= 1,]

  # B) Contains multiple apostrophe. Cultivar.
  with_mult_apostrophe = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "'") >= 2,]
  with_CV =  wcvp_names[stringr::str_count(wcvp_names$taxon_name, "CV") >= 1,]
  with_cv_end =  wcvp_names[stringr::str_count(wcvp_names$taxon_name, "$cv") >= 1,]
  with_cv_stop =  wcvp_names[stringr::str_count(wcvp_names$taxon_name, "cv\\.") >= 1,]

  # C) Is indeterminant by ' sp. '.
  with_spdot = wcvp_names[stringr::str_count(wcvp_names$taxon_name, " sp\\. ") >= 1,]

  # D) Is a hybrid by ' gx '.
  with_gx = wcvp_names[stringr::str_count(wcvp_names$taxon_name, " gx ") >= 1,]
  with_gx_no_space = wcvp_names[stringr::str_count(wcvp_names$taxon_name, " gx") >= 1,]

  # E) is indeterminant by starting with 'Indet ' or ending with 'indet'.
  with_indet = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "^Indet ") >= 1,]
  with_indet_end = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "$indet") >= 1,]

  # F) other.
  with_group = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "Group") >= 1,]
  with_unkn = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "unkn") >= 1,]
  with_hybrid_end = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "$hybrid") >= 1,]
  with_Hybrid_space = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "Hybrid ") >= 1,]
  with_Unknown = wcvp_names[stringr::str_count(wcvp_names$taxon_name, "Unknown") >= 1,]


  # Join all exceptions into one data.frame.
  exceptions = rbind(with_sqbracket2,
                     with_CV, with_cv_end, with_cv_stop, with_mult_apostrophe,
                     with_spdot,
                     with_gx, with_gx_no_space,
                     with_indet, with_indet_end,
                     with_group, with_unkn, with_hybrid_end, with_Hybrid_space, with_Unknown)

  #################################
  # 8) Log the changes.
  #################################
  changes = NULL

  # #add taxon name sanitising.
  # if(length(indices_require_sanitise)>0){
  #   changes = data.frame(issue = rep('Taxon name required sanitising',length(indices_require_sanitise)),
  #                        powo_id = wcvp_names$powo_id[indices_require_sanitise],
  #                        taxon_name = wcvp_names$taxon_name[indices_require_sanitise],
  #                        issue_entry = paste0('Original taxon_name = ', original_taxon_name[indices_require_sanitise]),
  #                        fix = paste0('New taxon_name = ', wcvp_names$taxon_name[indices_require_sanitise]))
  # }

  #add missing accepted form.
  indices = !is.na(new_accepted_name)
  if(sum(indices)>0){
    changesB = data.frame(issue = rep('synonym with missing accepted form',sum(indices)),
                         powo_id = non_accepted_synonyms$powo_id[indices],
                         taxon_name = non_accepted_synonyms$taxon_name[indices],
                         issue_entry = paste0('accepted_plant_name_id = ',non_accepted_synonyms$accepted_plant_name_id[indices]),
                         fix = paste0('accepted_plant_name_id = ',new_accepted_plant_name_id[indices], ' (',new_accepted_name[indices],')'))
    changes = rbind(changes, changesB)
  }

  #add square bracket issue.
  indices = with_sqbracket$taxon_name != with_sqbracket_original$taxon_name
  if(sum(indices)>0){
    changesB = data.frame(issue = rep('error in taxon name and taxon rank with [*] or [**]',sum(indices)),
                          powo_id = with_sqbracket_original$powo_id[indices],
                          taxon_name = with_sqbracket_original$taxon_name[indices],
                          issue_entry = paste0('taxon_name = ',with_sqbracket_original$taxon_name[indices], ' and taxon_rank = ',with_sqbracket_original$taxon_rank[indices]),
                          fix = paste0('taxon_name = ',with_sqbracket$taxon_name[indices], ' and taxon_rank = ',with_sqbracket$taxon_rank[indices]))
    changes = rbind(changes, changesB)
  }

  #################################
  # 9) return
  #################################
  return(list(wcvp_names = wcvp_names,
              exceptions = exceptions,
              changes = changes))
}

#' generate_labels()
#'
#'Generate shortened distribution details that can be used on labels in collections.
#'
#' @param data distribution data for plants taken from wcvp_distribution.
#'
#' @return labels
#' @export
#'
generate_labels <- function(data){

  #Vector to store the label.
  labels = rep(NA, nrow(data))

  ###
  #1) Reformat areas/regions/continents and count how many we have of each for all plants
  ###
  no_areas = unlist(lapply(data$area, function(x){
    stringr::str_count(x,', ') +1
  }))
  unique_regions= lapply(data$region, function(x){
    sort(unique(stringr::str_split(x,', ')[[1]]))
  })
  no_regions = unlist(lapply(unique_regions,length))
  unique_continents= lapply(data$continent, function(x){
    sort(unique(stringr::str_split(x,', ')[[1]]))
  })
  no_continents = unlist(lapply(unique_continents,length))

  ###
  # 2) Get the labels going through the following steps
  ###

  #A) If one regions.
  # Return: Endemic to region.
  indices = which(no_areas == 1)
  labels[indices] = paste0('Endemic to ',data$area[indices])

  #C) If two regions.
  # Return: RegionA, RegionB.
  indices_remain = which(is.na(labels))
  indices = indices_remain[no_areas[indices_remain] == 2]
  labels[indices] = data$area[indices]

  #D) If we have more than two areas but they all have the same region.
  # Return: the Region.
  indices_remain = which(is.na(labels))
  indices = indices_remain[no_regions[indices_remain] == 1]
  labels[indices] = unlist(unique_regions[indices])

  #E) If we have more than two areas and there is a total of 2 regions.
  # Return the two regions with ', ' in between.
  indices_remain = which(is.na(labels))
  indices = indices_remain[no_regions[indices_remain] == 2]
  labels[indices] = unlist(lapply(unique_regions[indices],function(x){paste0(x,collapse = ', ') }))

  #F) If we have more than two regions and they are contained in a one or two continent.
  # Return the continent or ContinentA, ContinentB
  indices_remain = which(is.na(labels))
  indices_to_POWO_remark = indices_remain
  indices = indices_remain[no_continents[indices_remain] < 2.5]
  labels[indices] = unlist(lapply(unique_continents[indices],function(x){paste0(x,collapse = ', ')}))
  #if label == "NORTHERN AMERICA, SOUTHERN AMERICA" set to OLD WORLD.
  labels[labels == "NORTHERN AMERICA, SOUTHERN AMERICA"] = 'NEW WORLD'

  #G) If we have more than two regions and they are contained in two continent.
  # Return the two continent with ', ' in between.
  indices_remain = which(is.na(labels))
  indices = indices_remain[no_continents[indices_remain] < 5]
  labels[indices] = unlist(lapply(unique_continents[indices],function(x){paste0(unique(x),collapse = ', ') }))
  # Any combination of 3 or more of below we set to OLD WORLD.
  Old_world = c("AFRICA", "ASIA-TEMPERATE", "ASIA-TROPICAL", "EUROPE")
  Old_world_combinations = unlist(lapply(combn(Old_world,3, simplify = F), function(x){
    paste0(sort(x),collapse = ', ')
  }))
  labels[labels %in% Old_world_combinations] = 'OLD WORLD'

  #H) Finally just return worldwide for the remaining (>5 continents)
  # Return Worldwide.
  indices_remain = which(is.na(labels))
  labels[indices_remain] = rep('WORLDWIDE', length(indices_remain))

  return(labels)
}



#' Add wcvp_distributions information to wcvp_names
#'
#' We allow three methods for uploading the information from wcvp_distributions:
#' - By choosing a file path to a csv file containing the information using `filepath`.
#' - By loading an R save (file type `.RData`, `.rdata`, `.rda` or `.rds`) via `filepath`.
#' - By using the package  `rWCVPData`, by setting `use_rWCVPData` = TRUE.
#'
#'
#'
#' @param filepath path to wcvp_distributions.csv file.
#' @param wcvp Output from import_wcvp_names
#' @param use_rWCVPdata Flag for whether we use rWCVPdata to get wcvp_distributions
#'
#' @return Desired information from wcvp_distributions.csv combined with wcvp_names.
#' @export
#'
add_wcvp_distributions <- function(filepath = NULL, wcvp, use_rWCVPdata = FALSE){

  ######################################################
  # 1) Load wcvp_distributions (using method as described in the download)
  ######################################################
  cli::cli_h1("Adding wcvp distribution information (wcvp_distributions)")

  cli::cli_h2("Loading wvcp_distribution")
  if(!is.null(filepath)){
    # We have a filepath, try and load depending on file format.
    if(grepl('\\.csv$',filepath)){
      wcvp_distributions = utils::read.table(filepath, sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
    }
    else if(grepl('\\.RData$|\\.rdata$|\\.rda',filepath)){
      load(filepath)
      wcvp_distributions = data.frame(wcvp_distributions) # To reformat to standard data frame not a tibble.
    }
    else if(grepl('\\.rds$',filepath)){
      readRDS(filepath)
      wcvp_distributions = data.frame(wcvp_distributions) # To reformat to standard data frame not a tibble.

    }
    else{
      stop('filepath is not a csv or R data file!')
    }

  }else if(use_rWCVPdata){
    if(system.file(package='rWCVPdata') == ''){
      stop('Package rWCVPdata is not installed please use another method to load wcvp_names or install the package')
    }
    if(requireNamespace('rWCVPdata', quietly = T)){
      wcvp_distributions <- rWCVPdata::wcvp_distributions
    }

  }
  if(!exists('wcvp_distributions')){
    stop('Have not successfully loaded wcvp_distributions, please make sure you provide either: filepath or set use_rWCVPdata = TRUE')
  }
  ######################################################
  # 2) Re-format data.
  ######################################################
  # Select only non-extinct, native, with no location doubt.
  wcvp_distributions$locationcombined = paste0(wcvp_distributions$introduced,wcvp_distributions$extinct, wcvp_distributions$location_doubtful)

  locationOptions = unique(wcvp_distributions$locationcombined)
  plant_ids = sort(unique(wcvp_distributions$plant_name_id),decreasing = FALSE)

  cli::cli_h2('Extracting plant geography for native/introduced, extinct and location doubtful')
    locations = pbapply::pblapply(locationOptions, function(locationOption){
    wcvp_distributions_cur = wcvp_distributions[which(wcvp_distributions$locationcombined == locationOption),]

    wcvp_distributions_plant_id <- wcvp_distributions_cur |>
      dplyr::group_by(.data$plant_name_id) |>
      dplyr::summarise(area = toString(.data$area),
                       region = toString(.data$region),
                       continent = toString(.data$continent),
                       area_code_l3 = toString(.data$area_code_l3)
                       # continent_code_l1 = toString(unique(.data$continent_code_l1)),
                       # region_code_l2 = toString(unique(.data$region_code_l2))
      ) |>
      dplyr::ungroup()

    #Get labels for naive extant with no location doubt.
    if(locationOption == '000'){
      labels = generate_labels(wcvp_distributions_plant_id)
      wcvp_distributions_plant_id = data.frame(wcvp_distributions_plant_id, labels = labels)
    }
    wcvp_distributions_plant_id = wcvp_distributions_plant_id[,-c(2,3,4)]
    # Now perform matching to wcvp using plant_name_id.
    match_across = match(plant_ids, wcvp_distributions_plant_id$plant_name_id)

    dist_info = data.frame(matrix(NA, nrow = length(plant_ids), ncol = ncol(wcvp_distributions_plant_id)))

    names(dist_info) = paste0('Dist_', locationOption, '_',names(wcvp_distributions_plant_id))
    indices = which(!is.na(match_across))
    dist_info[indices,] = wcvp_distributions_plant_id[match_across[!is.na(match_across)],]
    return(dist_info)
  })
  names(locations) = locationOptions

  dist_info = do.call(cbind, locations)

  dist_info = dist_info[-c(4,6,8,10,12)] #Remove excell plant name ids
  names(dist_info) = stringr::str_sub(names(dist_info),5,-1 )
  names(dist_info)[1] = 'plant_name_id'

  #Add dist_info to wcvp$wcvp_names
  wcvp$geography = data.frame(dist_info)
  return(wcvp)
}

