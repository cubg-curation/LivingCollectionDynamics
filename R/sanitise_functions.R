#' Sanitising functions
#'
#' Functions to clean and standardise taxonomic names and authors.
#'
#' `sanitise_name()` returns the sanitised name of a single taxonomic name.
#'
#' `sanitise_authors()` returns the sanitised name of a  taxonomic authors. Where by characters are coerced to Latin-ASCII, thereby removing diacritics (e.g umlauts).
#'
#' `clean_names_authors()` sanitises multiple taxonomic names with or without the corresponding authors, by applying `sanitise_name()` and `sanitise_authors()`. As input a vector of taxonomic names is required (`taxon_names`), in addition a vector of the authors (`taxon_authors`) or joined taxonomic name and authors (`taxon_names_full`) can be provided. If neither `taxon_authors` or `taxon_names_full` are provided the author names are set to `''`. A list is returned where:
#' -  `$taxon_name` a vector of the sanitised taxonomic names,
#' -  `$author`a vector of the sanitised authors,
#' -  `$sanitised` is a logical vector of whether the taxon_name was sanitised.
#'
#' `clean_names_authors_report()` applies `clean_names_authors()` to a collection where the inputs are:
#' -  a data frame of the collection (`collection`),
#' - column name for the taxonomic names (`taxon_name_column`, required),
#' - column name of the authors of the taxonomic names  (`taxon_author_column`, optional),
#' - column name of the combined taxonomic name and author  (`taxon_name_full_column`, optional).
#'
#' @param collection A data frame of a collection.
#' @param taxon_name_column The name of the column in the `collection` corresponding to taxonomic names.
#' @param taxon_name_full_column The name of the column in the `collection` corresponding to joined taxonomic names and authors.
#' @param taxon_author_column The name of the column in the `collection` corresponding to the authors of the taxonomic names.
#' @param taxon_name The taxonomic name of a plant.
#' @param taxon_names Vector of taxonomic names.
#' @param taxon_names_full Vector of joined taxonomic name and author.
#' @param taxon_authors Vector of taxonomic authors.
#' @param console_message Flag (TRUE/FALSE) for showing progress bar in the console.
#'
#' @export
#'
#' @examples
#' sanitise_name('TRIGONELLA afghanica')
#' sanitise_name('Halimium X pauanum')
#' sanitise_name('Aruncus dioicus var acuminatus')
#' sanitise_authors('Stehlé')
#'
#' taxon_names = c('TRIGONELLA afghanica', 'Halimium X pauanum',
#'  'Aruncus dioicus var acuminatus', 'Eupatorium magdalenae')
#' taxon_authors = c('Vassilcz', 'Font Quer', '(Douglas ex Hook.) H.Hara', 'Stehlé')
#'
#' sanitise_names_authors(taxon_names, taxon_authors)
#'
#' collection = data.frame(names = taxon_names, full = paste0(taxon_names, ' ', taxon_authors))
#' sanitise_names_authors_report(collection, taxon_name_column = 'names',
#'  taxon_name_full_column = 'full')
sanitise_name <- function(taxon_name){

  if(is.na(taxon_name)){
    return('NA')
  }
  ### 1) fix x/X/h/H to \u00D7.
  if(grepl(' [xXhH] |^[xXhH] | [xXhH]$',taxon_name)){
    # Change the symbol
    taxon_name = stringr::str_replace(taxon_name,'^[xX] ','\u00D7 ')
    taxon_name = stringr::str_replace(taxon_name,' [xX]$',' \u00D7')
    taxon_name = stringr::str_replace(taxon_name,' [xX] ',' \u00D7 ')
  }

  ### 2) Make sure there are spaces surrounding the hybrid marker
  if(grepl('\u00D7|\\+', taxon_name)){
    # Make sure the symbol is preceded by and followed by a space if not at the start or end of the taxon name.
    length_taxon = stringr::str_length(taxon_name)
    locations = as.numeric(stringr::str_locate_all(taxon_name, '\u00D7|\\+')[[1]][,1])
    locations = locations[!locations %in% c(1,length_taxon)]
    if(length(locations) > 0){
      for(i in 1:length(locations)){
        # Position before a space?
        if(stringr::str_sub(taxon_name, locations[i]-1, locations[i]-1) != ' '){
          taxon_name = paste0(stringr::str_sub(taxon_name, 1, locations[i]-1),
                              ' ',
                              stringr::str_sub(taxon_name, locations[i], length_taxon))
          locations[i] = locations[i] + 1
          length_taxon = length_taxon + 1
        }
        # Position after a space?
        if(stringr::str_sub(taxon_name, locations[i]+1, locations[i]+1) != ' '){
          taxon_name = paste0(stringr::str_sub(taxon_name, 1, locations[i]),
                              ' ',
                              stringr::str_sub(taxon_name, locations[i]+1, length_taxon))
          length_taxon = length_taxon + 1

        }
      }
    }
  }

  ### 3) Fix casing only first letter of Genus Capital rest lower
  # If the starting letter means hybrid (x,+,\u00D7)
  if(grepl('^[+\u00D7]',taxon_name)){
    taxon_part = stringr::str_sub(taxon_name,3,-1)
    taxon_part = gsub("(\\D)(\\D+)", "\\U\\1\\L\\2", taxon_part, perl = TRUE)
    taxon_name = paste0(stringr::str_sub(taxon_name,1,2), taxon_part, collapse ='')
  }
  else{
    # First letter capital the rest lower.
    taxon_name = gsub("(\\D)(\\D+)", "\\U\\1\\L\\2", taxon_name, perl = TRUE)
  }

  ### 4) fix f or var to f. and var.
  if(grepl(' f | var | subsp | v ',taxon_name)){
    taxon_name = stringr::str_replace(taxon_name,' f ',' f\\. ')
    taxon_name = stringr::str_replace(taxon_name,' var | v ',' var\\. ')
    taxon_name = stringr::str_replace(taxon_name,' subsp ',' subsp\\. ')
    taxon_name = stringr::str_replace(taxon_name,' nothosubsp ',' nothosubsp\\. ')

  }

  ### 5) Remove excess whitespace.
  taxon_name = stringr::str_squish(taxon_name)

  ### 6) Simplify all special characters (except hybrid marker)
  words = stringr::str_split(taxon_name, ' ')[[1]]
  taxon_name = paste0(unlist(lapply(words, function(word){
    if(!grepl('\u00D7|\\+', word)){
      word = stringi::stri_trans_general(word, id = "Latin-ASCII")
    }
    word
  })), collapse = ' ')
  return(taxon_name)
}

#' @rdname sanitise_name
#' @export
sanitise_authors <- function(taxon_authors){
  stringi::stri_trans_general(taxon_authors, id = "Latin-ASCII") # simplify characters, i.e remove upstroph, tilde.
}


#' @rdname sanitise_name
#' @export
sanitise_names_authors <- function(taxon_names,
                                taxon_authors = NA,
                                taxon_names_full = NA,
                                console_message = FALSE){
  # A) Sanitise the taxon names.
  if(console_message){
    clean_taxon_name = unlist(pbapply::pblapply(taxon_names, sanitise_name))
  }else{
    clean_taxon_name = unlist(lapply(taxon_names, sanitise_name))

  }

  # B) Extract author if needed.
  # i) Both taxon_authors and taxon_names_full = NA
  if(length(taxon_authors) == 1  & length(taxon_names_full) == 1 & all(is.na(taxon_authors)) & all(is.na(taxon_names_full))){
    author = rep('',length(taxon_names))
  }
  # ii) Both taxon_authors is NA  and taxon_names_full is not
  else if(length(taxon_authors) == 1  & length(taxon_names_full) > 1 & all(is.na(taxon_authors))){
    author = author_from_taxon_name_full(taxon_names, taxon_names_full)
    author = sanitise_authors(author)
  }
  #iii) Taxon authors in not NA.
  else{
    author = taxon_authors
    author = sanitise_authors(author)
  }


  #C) Was sanitising needed.
  sanitised = rep(F,length(taxon_names))
  sanitised[taxon_names != clean_taxon_name] = T


  return(list(taxon_name = clean_taxon_name, author = author, sanitised = sanitised))
}

#' @rdname sanitise_name
#' @export
sanitise_names_authors_report <- function(collection,
                                       taxon_name_column = 'TaxonName',
                                       taxon_name_full_column = NA,
                                       taxon_author_column = NA,
                                       console_message = FALSE){
  # Get the values out of original report.
  taxon_names = collection[[taxon_name_column]]
  if(is.na(taxon_name_full_column)){
    taxon_name_full = NA
  }else{
    taxon_name_full = collection[[taxon_name_full_column]]
  }

  if(is.na(taxon_author_column)){
    taxon_authors = NA
  }else{
    taxon_authors = collection[[taxon_author_column]]
  }

  return(sanitise_names_authors(taxon_names = taxon_names, taxon_names_full = taxon_name_full, taxon_authors = taxon_authors, console_message = console_message))
}
