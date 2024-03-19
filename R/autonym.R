#' Find if a taxonomic name is an autonym
#'
#' These function find when a  taxonomic name is an autonym.
#'
#'To obtain whether or not taxonomic names are autonyms we compare the words either side of a 'splitter'. The splitters are defined to be any of `subsp.`, `var.`, `f.`, `ssp.` or `nothosubsp.`. If the words either side are identical we return `TRUE` otherwise we return `FALSE.`
#'
#' @param taxon_name The taxonomic name of a plant.
#' @param collection A data frame containing a collection.
#' @param taxon_name_column The name of the column containing the (original) Taxon name.
#' @param progress_bar Logical flag, if TRUE show progress bar.

#' @return TRUE if the taxonomic name is an autonym, otherwise FALSE.
#'
#' @examples
#' is_autonym("Codiaeum variegatum var. variegatum")
#' is_autonym("Crinum pedunculatum f. purple")
#'
#' taxon_names = c("Saintpaulia diplotricha", "Codiaeum variegatum var. variegatum")
#' collection = data.frame(ID = 1:2, name = taxon_names)
#' add_is_autonym(collection, taxon_name_column = 'name')
#' @export
is_autonym <- function(taxon_name){
  # 1) Is the plant a hybrid, i.e contains '×' (unicode \u00D7)
  if(grepl('\u00D7',taxon_name)){
    return(FALSE)
  }

  # 2) Split word by level i.e var. f., etc.
  # And 'squish' the two resultant parts (i.e remove excess whitespace)
  split_name = stringr::str_split(taxon_name,' var\\. | subsp\\. | f\\. | ssp\\. | nothosubsp\\. ')[[1]]
  split_name = unlist(lapply(split_name, stringr::str_squish))

  # 3) If there is only one chunk return no (i.e no var., f., etc)
  if(length(split_name) < 2){return(FALSE)}


  # 4) Split each chunk into words.
  split_parts = stringr::str_split(split_name,' ')

  # 5) Loop over each pair of chunks and check for autonym.
  for(i in 1:(length(split_name)-1)){
    # If the length of the second part is more than one word then cultivars (i.e will most likely contain '' or []).
    if(length(split_parts[[i+1]]) > 1){ return(FALSE)}

    #Get the word either side of the split (var./f./subsp.)
    pre_word = split_parts[[i]][length(split_parts[[i]])]
    post_word = split_parts[[i+1]][1]

    # If the words don't math return false.
    if(pre_word != post_word){
      return(FALSE)
    }

  }
  # Each chunk matches so return true.
  return(TRUE)

}

#' @rdname is_autonym
#' @export
add_is_autonym <- function(collection, taxon_name_column = 'TaxonName', progress_bar = FALSE){
  if(progress_bar){
    autonyms = unlist(pbapply::pblapply(collection[,taxon_name_column], is_autonym))
  }
  else{
    autonyms = unlist(lapply(collection[,taxon_name_column], is_autonym))
  }
  data_new = data.frame(collection, is_autonym = autonyms)
  return(data_new)
}
