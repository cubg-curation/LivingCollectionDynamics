#' Extract the taxon type (biological/horticultural characteristics)
#'
#' These functions extract the taxon type from a taxon name, in particular whether an item is: an  indeterminate, a species, a subspecies, a variety, a forma, a cultivar or a hybrid.
#'
#'
#' The aim of these functions to to extract the biological and horticultural characteristics of items via their taxonomic name. For each taxonomic name the functions will return a string containing which taxon types it has. These are:
#' - `"0.indet"`: The item is indeterminate.
#' - `"1.species"`: The item is a species.
#' - `"2.subsp"`: The item is a subspecies.
#' - `"3.var"`: The item is a variety.
#' - `"4.f"`: The item is a forma.
#' - `"5.cultivar"`: The item is a cultivar.
#' - `"6.hybrid"`: The item is a hybrid.
#'
#'  `0-4` are biological characteristics and each item can only have one of these characteristics.
#'
#'  `5-6` are horticultural characteristics and item can have none, one or both of these.
#'
#' `taxon_type()` returns the taxon type given a single taxonomic name.
#'
#' `add_taxon_type()` adjoins a column called `taxon_type` to the original `collection` containing the taxon type of each item. If `is.na(POWO_taxon_name_column) == TRUE` then the taxon type is extracted purely from the taxon names contained in the column `taxon_name_column`. If `POWO_taxon_name_column` is a name of a column in `collection` then the taxon type is extracted from the column corresponding to it, if there are any missing values then the entries corresponding to `taxon_name_column` are used instead.
#'
#' @param taxon_name The taxonomic name of a plant.
#' @param collection A data frame containing a collection.
#' @param taxon_name_column The name of the column containing the (original) taxon name.
#' @param POWO_taxon_name_column The name of the column containing the POWO (WCVP) Taxon name.
#' @param progress_bar Logical flag, if TRUE show progress bar.
#' @export
#'
#'
#'
#' @examples
#' taxon_type("Aridaria sp.")
#' taxon_type("Saxegothaea conspicua")
#' taxon_type("Rhododendron charitopes subsp. tsangpoense")
#' taxon_type("Aquilegia flabellata f. alba 'green'")
#'
#' taxon_names = c("Saintpaulia diplotricha", "Aridaria sp.")
#' POWO_taxon_name = c("Streptocarpus ionanthus var. diplotrichus", NA)
#' collection = data.frame(name = taxon_names, POWO = POWO_taxon_name)
#' add_taxon_type(collection, taxon_name_column = 'name', POWO_taxon_name_column = 'POWO')
taxon_type <- function(taxon_name){

  if(is.na(taxon_name)){return(NA)}
  if(taxon_name ==''){return(NA)}
  # A) Split the taxon name into individual words.
  no_words = stringr::str_count(stringr::str_squish(taxon_name),' ') +1
  if(is.na(no_words)){return(NA)}

  # B) groups stores the infrageneric_levels that are found to describe taxon name
  groups = NULL

  # C) Go through each infrageneric level and add to groups if the correct pattern matches.

  # 0) "0.Indet"
  if(no_words == 1 | grepl(" sp\\.|indet\\.|unkn| cf\\.|aff\\.|spec\\.|Indet\\.|^Indet ",taxon_name)){
    groups = c(groups, '0.indet')
  }

  # 1) "1.species".
  if(no_words == 2){
    groups = c(groups, '1.species')
  }

  # 2) "2.subsp"
  if (grepl('subsp\\.', taxon_name)){
    groups = c(groups, '2.subsp')
  }

  # 3) "3.var"
  if (grepl('var\\.', taxon_name)){
    groups = c(groups, '3.var')
  }

  # 4) "4.f"
  if (grepl(' f\\.', taxon_name)){
    groups = c(groups, '4.f')
  }

  # 5) "5.cultivar"
  if (grepl("cv\\.|'.*?'|CV|cv$|\\[.*?\\]|hort\\.", taxon_name)){
    groups = c(groups, '5.cultivar')
  }

  # 6) "6.hybrid"
  if (grepl('hybrid$|Hybrid|HYBRID|gx\\.| gx | gx|\u00D7', taxon_name)){
    groups = c(groups, '6.hybrid')
  }

  # 6) "6.hybrid" (this checkes the hex code rather than unicode for 'x' not sure if this will be needed once encoding issues are sorted)
  if (stringr::str_detect(taxon_name,'\xd7')){
    groups = c(groups, '6.hybrid')
  }

  groups = unique(groups)
  if(is.null(groups)){
    return(NA)}

  # D) If we have multiple biological characteristics we need to remove excess ones
  # Such as choosing the biological characteristic with the smallest scope.
  # or removing species if we also have cultivar (i.e taxonname length  = 2 but one of the words will be cv. etc)
  if(length(groups) >1){
    # if we have 4.f remove 3.var and 2.subs
    if('4.f' %in% groups){
      groups = groups[!grepl('2|3',groups)]
    }
    # if we have 3.var remove 2.subs
    if('3.var' %in% groups){
      groups = groups[!grepl('2',groups)]
    }
    # if we have 0.indet then remove biological characteristics (i.e "Astragalus sp.", "Mammillaria sp. f. cristata" )
    if('0.indet' %in% groups){
      groups = groups[!grepl('1|2|3|4',groups)]
    }
    # If we have cultivar and species then this is a fake species i.e length = 2 but contains cv or 'XX' so remove species.
    if('5.cultivar' %in% groups){
      groups = groups[!grepl('1',groups)]
    }
    # If we have hybrid and species then this is a fake species i.e length = 2 but contains gx so remove species.
    if('6.hybrid' %in% groups){
      groups = groups[!grepl('1',groups)]
    }
    # If we find [HYBRID] this is not a cultivar so remove.
    if(grepl('\\[HYBRID\\]', taxon_name)){
      groups = groups[!grepl('5',groups)]
    }
    # If we have indet and one hort category return only indet.
    if('0.indet' %in% groups & length(groups) ==2){
      groups = groups[!grepl('5|6',groups)]
    }
  }

  return(paste0(groups,collapse = ', '))
}




#' @rdname taxon_type
#' @export
add_taxon_type <- function(collection, taxon_name_column = 'TaxonName', POWO_taxon_name_column = NA, progress_bar = FALSE){
  ####
  # 1) Get taxon names to use for infrageneric_level
  ####
  #A) we DON'T have POWO_taxon_name_column.
  if(is.na(POWO_taxon_name_column)){
    taxon_names = collection[[taxon_name_column]]
  }
  #B) we DO have POWO_taxon_name_column.
  else{
    taxon_names = collection[[POWO_taxon_name_column]]
    original_taxon_names = collection[[taxon_name_column]]
    taxon_names[is.na(taxon_names)] = original_taxon_names[is.na(taxon_names)]
  }

  ###
  # 2) Get the infrageneric level.
  ###
  if(progress_bar){
    infrageneric_levels = unlist(pbapply::pblapply(taxon_names, taxon_type))
  }
  else{
    infrageneric_levels = unlist(lapply(taxon_names, taxon_type))
  }

  ###
  # 3) Return data with infrageneric_levels added.
  ###
  return(data.frame(collection, taxon_type = infrageneric_levels))
}
