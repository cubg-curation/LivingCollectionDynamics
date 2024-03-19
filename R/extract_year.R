#' Extract the year from different date formats
#'
#' A function that extracts the year from dates where the date can be in such formats as: `DD/MM/YYYY`, `YYYY-MM-DD`, `MM/YYYY`, `YYYY-MM`, `YYYY`, `DD Month YYYY`.
#'
#' The method works by extracting any patterns of 4 numbers in a row.
#' Thus, the year must be in the format YYYY. I.e. cannot have formats such as DD/MM/YY.
#'
#' @param dates A vector of dates.
#'
#' @return A vector of corresponding years.
#' @export
#' @examples
#' extract_year(c('03/05/2021', '3rd of May 2022', '2023-05-03'))
extract_year <- function(dates){
  return(as.numeric(unlist(stringr::str_extract(dates, '[0-9]{4}'))))
}


#' Add status year column to collection
#'
#' A function to adjoin the year from the item status date to the collection database.
#'
#' @param collection A data frame containing a collection.
#' @param item_status_date_column The name of the column containing the item status date.
#'
#' @return The `collection` data frame with `status_year` adjoined.
#' @export
#'
#' @examples
#' taxon_names = c('Trigonella afghanica', 'Eupatorium magdalenae')
#' item_status = c('2021-04-12', '12th August 1991')
#' collection = data.frame(name = taxon_names, item_status = item_status)
#' add_status_year(collection, item_status_date_column ='item_status')
add_status_year <- function(collection, item_status_date_column = 'ItemStatusDate'){
  ItemStatusdate = collection[,match(item_status_date_column, names(collection))]
  status_year = extract_year(ItemStatusdate)
  return(data.frame(collection, status_year = status_year))
}
