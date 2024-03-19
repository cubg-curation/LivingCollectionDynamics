#' Enrich a collection's database
#'
#' This function enriches plant records in a collection using information from POWO (WCVP), IUCN redlist and BGCI.
#'
#'
#'
#'
#' @param collection A data frame containing a collection.
#' @param wcvp World Checklist of Vascular Plants (WCVP) database, obtained using the function [import_wcvp_names()]. If `NA` WCVP enrichment is not performed.
#' @param iucnRedlist IUCN Red List of Threatened Species database, obtained using the function XXXX.  If `NA` IUCN red list enrichment is not performed.
#' @param BGCI Requires a cleaned BGCI plant search database to obtain how many collections globally a taxon is contained (Not freely available). If `NA` BGCI plant search enrichment is not performed.
#' @param taxon_name_column The name of the column in the `collection` corresponding to taxonomic names.
#' @param taxon_name_full_column The name of the column in the `collection` corresponding to joined taxonomic names and authors.
#' @param taxon_author_column The name of the column in the `collection` corresponding to the authors of the taxonomic names.
#' @param do_is_autonym Flag (TRUE/FALSE) for whether to add the column is_autonym, see [add_is_autonym()] for details.
#' @param do_status_year Flag (TRUE/FALSE) for whether to add the column status_year, see [add_status_year()] for details.
#' @param do_taxon_types Flag (TRUE/FALSE) for whether to add the column taxon_type, see [add_taxon_type()] for details.
#' @param ... Arguments (i.e., attributes) used in the matching algorithm (passed along to [match_collection_to_wcvp()]). Examples include `typo_method`, `do_convert_accepted` and `try_fix_hybrid`.
#' @param wcvp_wanted_info A character vector containing the information we want to extract from WCVP. These are matched to column names in `wcvp$wcvp_names`, `wcvp$geography` and `wcvp$redlist` (if the later two exist).
#' @param redlist_wanted_info A character vector containing the information we want to extract from IUCN red list, corresponding to column names in `iucnRedlist`.
#'
#' @return The `collection` data frame enriched with information dependent on function inputs (new columns).
#' @export
enrich_collection <- function(collection,
                          wcvp = NA,
                          iucnRedlist = NA,
                          BGCI = NA,
                          taxon_name_column = 'TaxonName',
                          taxon_name_full_column = NA,
                          taxon_author_column = NA,
                          do_is_autonym = FALSE,
                          do_status_year = FALSE,
                          do_taxon_types = FALSE,
                          ...,
                          wcvp_wanted_info = c("plant_name_id", "taxon_name", "taxon_authors", "taxon_rank", "taxon_status","powo_id", "family", "genus", "species", "lifeform_description", "climate_description", "geographic_area", "Dist_000_area_code_l3", "Dist_000_labels", "Dist_100_area_code_l3", "Dist_010_area_code_l3", "Dist_001_area_code_l3", "Dist_101_area_code_l3", "Dist_110_area_code_l3", "Dist_011_area_code_l3", "main_common_name", "assessment_date", "category", "criteria", "population_trend"),
                          redlist_wanted_info = c("plant_name_id", "main_common_name", "assessment_date", "category",           "criteria", "population_trend")
                          ){

  ############################################
  # 1) Clean/extract taxon name + taxon author
  ############################################
  cli::cli_h2("Sanitise taxon name and extract author")
  taxon_name_and_author = sanitise_names_authors_report(collection,
                             taxon_name_column = taxon_name_column,
                             taxon_name_full_column = taxon_name_full_column,
                             taxon_author_column = taxon_author_column)


  collection = data.frame(collection,
                               sanitised_taxon = taxon_name_and_author$taxon_name,
                               need_sanitise = taxon_name_and_author$sanitised,
                               extracted_author = taxon_name_and_author$author)
  enriched_report = collection


  ############################################
  # 1) Add is_autonym.
  ############################################
  if(do_is_autonym){
    cli::cli_h2("Adding is autonym")
    enriched_report = add_is_autonym(enriched_report)
  }

  ############################################
  # 2) Add status_year.
  ############################################
  if(do_status_year){
    cli::cli_h2("Adding status year")
    enriched_report = add_status_year(enriched_report)
  }


  ############################################
  # 4) Add information from POWO.
  ############################################
  if(!is.na(wcvp)[1]){
  cli::cli_h2("Adding POWO information")

  # A) find the match between original report and wcvp.
  match_info = match_collection_to_wcvp(collection,
                                      wcvp,
                                      taxon_name_column = 'sanitised_taxon',
                                      taxon_name_full_column = NA,
                                      taxon_author_column = 'extracted_author',
                                      enrich_taxon_name_column = 'taxon_name',
                                      enrich_display_in_message_column = 'powo_id',
                                      enrich_plant_identifier_column = 'plant_name_id',
                                      ...)
  enriched_report = data.frame(enriched_report,
                               POWO_original_authors = match_info$original_authors,
                               POWO_match_taxon_name = match_info$match_taxon_name,
                               POWO_match_authors = match_info$match_authors,
                               POWO_author_check = match_info$author_check,
                               POWO_match_detail = match_info$details,
                               POWO_match_detail_short =match_info$details_short)


  # B) Force plant_name_id to be in wanted info as it's needed for geography and red list info.
  wcvp_wanted_info_orig = wcvp_wanted_info
  wcvp_wanted_info = unique(c('plant_name_id', wcvp_wanted_info))

  # C) Extract the info from wcvp_names.
  wcvp_wanted_columns = wcvp_wanted_info
  wcvp_wanted_columns = wcvp_wanted_columns[wcvp_wanted_columns %in% names(wcvp$wcvp_names)]
  POWO_info = data.frame(matrix(NA, nrow = nrow(collection), ncol = length(wcvp_wanted_columns)))
  names(POWO_info) = paste0('POWO_',wcvp_wanted_columns)
  indices = which(!(is.na(match_info$match) | match_info$match < 0))
  POWO_info[indices,] = wcvp$wcvp_names[match_info$match[indices],match(wcvp_wanted_columns,names(wcvp$wcvp_names))]
  if('plant_name_id' %in% wcvp_wanted_info_orig){
    enriched_report = data.frame(enriched_report, POWO_info)
  }else{
    enriched_report = data.frame(enriched_report, POWO_info[,-1])  #plant_name_id always the first column.
  }

  # D) Add info from wcvp_distributions.
  if('geography' %in% names(wcvp)){
    wcvp_wanted_columns = wcvp_wanted_info
    wcvp_wanted_columns = wcvp_wanted_columns[wcvp_wanted_columns %in% names(wcvp$geography)]
    if(length(wcvp_wanted_columns)>1){
      Geog_info = data.frame(matrix(NA, nrow = nrow(collection), ncol = length(wcvp_wanted_columns)))
      havePOWO_index = which(!is.na(POWO_info$POWO_plant_name_id))
      plant_id_match = match(POWO_info$POWO_plant_name_id[havePOWO_index], wcvp$geography$plant_name_id)

      Geog_info[havePOWO_index,] = wcvp$geography[plant_id_match,match(wcvp_wanted_columns,names(wcvp$geography))]
      names(Geog_info) = paste0('POWO_', wcvp_wanted_columns)
      Geog_info = Geog_info[,-1] # remove plant_name_id column.

      enriched_report = data.frame(enriched_report, Geog_info)

      if(length(wcvp_wanted_columns)==2){
        names(enriched_report)[length(names(enriched_report))] = paste0('POWO_',wcvp_wanted_columns[2])
      }
    }

  }

  # E) Add info from WCVP matched to Redlist.
  if('redList' %in% names(wcvp)){
    wcvp_wanted_columns = wcvp_wanted_info
    wcvp_wanted_columns = wcvp_wanted_columns[wcvp_wanted_columns %in% names(wcvp$redList)]
    if(length(wcvp_wanted_columns)>1){
      Red_info = data.frame(matrix(NA, nrow = nrow(collection), ncol = length(wcvp_wanted_columns)))
      havePOWO_index = which(!is.na(POWO_info$POWO_plant_name_id))
      plant_id_match = match(POWO_info$POWO_plant_name_id[havePOWO_index], wcvp$redList$plant_name_id)

      Red_info[havePOWO_index,] = wcvp$redList[plant_id_match,match(wcvp_wanted_columns,names(wcvp$redList))]
      names(Red_info) = paste0('POWO_Red_',wcvp_wanted_columns)
      Red_info = Red_info[,-1] # remove plant_name_id column.

      enriched_report = data.frame(enriched_report, Red_info)

      if(length(wcvp_wanted_columns)==2){
        names(enriched_report)[length(names(enriched_report))] = paste0('POWO_Red_',wcvp_wanted_columns[2])
      }

    }


  }

  # F) Create POWO web address.
  if('powo_id' %in% wcvp_wanted_info){
    POWO_web_address = rep(NA, nrow(POWO_info))
    indices = !is.na(POWO_info$POWO_powo_id)
    POWO_web_address[indices] = paste0('https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:',POWO_info$POWO_powo_id[indices])

    enriched_report = data.frame(enriched_report, POWO_web_address = POWO_web_address)

  }
  }

  ############################################
  # 4) Add taxon type
  ############################################
  # We run infrageneric_level on the taxonName from POWO unless we didn't find a match then we use the original taxon name.
  if(do_taxon_types){
    cli::cli_h2("Adding taxon types")
    if(!is.na(wcvp)[1]){
    enriched_report = add_taxon_type(enriched_report, POWO_taxon_name_column = 'POWO_taxon_name')
    }else{
      enriched_report = add_taxon_type(enriched_report, POWO_taxon_name_column = NA)

    }
  }

  ############################################
  # 5) Add iucnRedlist information.
  ############################################
  if(!is.na(iucnRedlist)[1]){
    cli::cli_h2("Adding ICUN Red list information")

    if(!exists('enrich_taxon_name_column')){
      enrich_taxon_name_column = 'scientific_name'
    }
    if(!exists('enrich_plant_identifier_column')){
      enrich_plant_identifier_column = 'taxonid'
    }
    if(!exists('enrich_display_in_message_column')){
      enrich_display_in_message_column = 'taxonid'
    }

    match_info = match_collection_to_iucnRedlist(collection,
                                                 iucnRedlist = iucnRedlist,
                                        taxon_name_column = 'sanitised_taxon',
                                        taxon_name_full_column = NA,
                                        taxon_author_column = 'extracted_author',
                                        enrich_taxon_name_column =  enrich_taxon_name_column,
                                        enrich_plant_identifier_column = enrich_plant_identifier_column,
                                        enrich_display_in_message_column = enrich_display_in_message_column,
                                        ...)

    redList_wanted_columns = redlist_wanted_info
    redList_wanted_columns = redList_wanted_columns[redList_wanted_columns %in% names(iucnRedlist)]
    redList_info = data.frame(matrix(NA, nrow = nrow(collection), ncol = length(redList_wanted_columns)))
    names(redList_info) = paste0('redList_',redList_wanted_columns)
    indices = which(!(is.na(match_info$match) | match_info$match < 0))
    redList_info[indices,] = iucnRedlist[match_info$match[indices],match(redList_wanted_columns,names(iucnRedlist))]

    enriched_report = data.frame(enriched_report,
                                 redList_original_authors = match_info$original_authors,
                                 redList_match_taxon_name = match_info$match_taxon_name,
                                 redList_match_authors = match_info$match_authors,
                                 redList_author_check = match_info$author_check,
                                 redList_match_detail = match_info$details,
                                 redList_match_detail_short =match_info$details_short,
                                 redList_info)
  }

  ############################################
  # 5) Add Number of Gardens (BGCI).
  ############################################
  if(!is.na(BGCI)[1]){
    cli::cli_h2("Adding Number of gardens")

    no_gardens = match_original_to_BGCI(collection,
                                        BGCI,
                                        taxon_name_col = 'sanitised_taxon')

    enriched_report = data.frame(enriched_report,
                                 no_gardens = no_gardens$no_gardens)
  }

  #Return the enriched report.
  return(enriched_report)
}
