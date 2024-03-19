#' Separate author from combined taxonomic name and author
#'
#' Extract the taxonomic author from combined taxonomic name and author combined.
#'
#' `taxon_names` and `taxon_names_full` must have the same number of elements where the ith element of each correspond to the same name. To obtain the author of ith record and word/s found in `taxon_names[i]` are removed from  `taxon_names_full[i]`.
#'
#' @param taxon_names Vector of taxonomic names.
#' @param taxon_names_full Vector of joined taxonomic name and author.
#'
#' @return Vector of the author/s of the taxonomic names.
#' @export
#'
#' @examples
#' taxon_names = c('Achatocarpus praecox', 'Anthurium affine')
#' taxon_names_full = c('Achatocarpus praecox Griseb.', 'Anthurium affine Schott')
#' author_from_taxon_name_full(taxon_names, taxon_names_full)
author_from_taxon_name_full <- function(taxon_names, taxon_names_full){
  # Remove all backslashes
  taxon_names = stringr::str_remove_all(taxon_names, pattern = '\\\\')
  taxon_names_full = stringr::str_remove_all(taxon_names_full, pattern = '\\\\')

  #As we're using grepl need to add escape for special characters.
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\.', '\\\\.')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\[', '\\\\[')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\]', '\\\\]')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\+', '\\\\+')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\?', '\\\\?')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\(', '\\\\(')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\)', '\\\\)')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\*', '\\\\*')
  taxon_names = stringr::str_replace_all(taxon_names,pattern = '\\*', '\\\\*')

  # Convert taxon name to a grepl statement, where AA BB goes to AA|BB.
  taxon_name_words_grepl = unlist(lapply(taxon_names, function(x){
    words = stringr::str_split(x, ' ')[[1]]
    return(paste0(words,collapse='|'))
  }))

  #Loop over all taxon name full removing any word that is also in taxon name.
  authors = rep(NA,length(taxon_names_full))
  for(i in 1:length(authors)){
    if(is.na(taxon_names_full[i])){
      authors[i] = ''
    }
    else if(taxon_names_full[i] != '' | taxon_name_words_grepl[i] != ''){
      auth_cur = stringr::str_replace_all(taxon_names_full[i],taxon_name_words_grepl[i],'')
      authors[i] = stringr::str_squish(auth_cur)
    }
    else{
      authors[i] = ''
    }

  }
  # authors = stringi::stri_trans_general(authors, id = "Latin-ASCII") # simplify characters, i.e remove upstroph, tilde.
  return(authors)
}




#' Compare two authors
#'
#' @details
#' The result from the comparison is:
#' - `"Identical"` if the two author strings are equal (`original_author == proposed_author`),
#' - `"Partial"` if any of the words in either author is found in the other. Author words are found using [author_words()].
#' - `"Different"` otherwise.
#'
#' @param original_author,proposed_author The authors to be compared.
#'
#' @return Either `Identical', 'Partial' or 'Different'.
#' @export
#'
#' @examples
#' author_check("Schott", "Schott")
#' author_check("Scott", "Schott")
#' author_check("(Jacq.) Schott", "Schott")
author_check <- function(original_author, proposed_author){
  if(is.na(original_author) || is.na(proposed_author)){
    return('Different')
  }
  # Check if the authors are identical.
  if(original_author == proposed_author){
    return('Identical')
  }
  # strip each author into words.
  original_words = author_words(original_author)
  proposed_words = author_words(proposed_author)
  orig_words_in_proposed = FALSE
  prop_words_in_original = FALSE
  if(length(original_words) > 0){
    orig_words_in_proposed = unlist(lapply(original_words, function(x){grepl(x,proposed_author)}))
  }
  if(length(proposed_words) > 0){
    prop_words_in_original= unlist(lapply(proposed_words, function(x){grepl(x,original_author)}))
  }
  if(any(c(orig_words_in_proposed, prop_words_in_original))){
    return('Partial')
  }

  return('Different')
}

#' Extract words from author
#'
#' @details
#' From the string `author` we extract any words that begin with a capital and is followed by at least two lower case letters or hyphens. Or the special cases `DC.`, `Sm.` and `Br.`. In the special cases `\\` will be added before `.` to allow regular expression searching of the authors such as `grepl("DC\\.", "Some Author")`. We have not included `L.` since this could be mistaken for an initial of another author rather than Carl Linnaeus.
#'
#' This function is used when comparing authors in [author_check()].
#'
#' @param author Author that wants splitting into words.
#'
#' @return Words contained in the author.
#' @export
#' @examples
#' author_words('(Jacq.) Schott')
#' author_words('Villarroel & J.R.I.Wood')
#' author_words('(DC.) F.Muell.')
author_words <- function(author){
  author_wordsA = unlist(stringr::str_extract_all(author,'[A-Z]{1}[a-z-]{2,}'))
  author_wordsB = unlist(stringr::str_extract_all(author,'DC\\.|Sm\\.|Br\\.'))
  author_wordsB = stringr::str_replace_all(author_wordsB,'\\.','\\\\.')
  return(c(author_wordsA,author_wordsB))
}
