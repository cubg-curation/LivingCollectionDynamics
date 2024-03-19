
#' Generate data frame of records that exist each year
#'
#' This function calculates which items in a collection are existing at a given date (or dates).
#'
#'  Accession date is inputted through `AccessionYear` which accepts only the year, for example `c(2014,2017,2015)`. The accession date is then set to the 1st of January in the given year.
#'
#'  Item death is calculated via the inputs `ItemStatusDate` and `ItemStatusType`. If the item is still alive (`ItemStatusType == 'Existing`) then the death date is set to today (`Sys.Date()`). Otherwise the item is dead (`ItemStatusType == 'NotExisting`) and we use the corresponding date given in `ItemStatusDate`. If only the year is available in `ItemStatusDate` we set the date of death to the 31st of December that year. If only the month and year is available we set the day of the item death to be the 28th.
#'
#'  `AccessionYear`, `ItemStatusDate` and `ItemStatusType` inputs must have the same length and the ith value in each corresponds to the ith item in the collection.
#'
#' @param date dates in the format `YYYY-MM-DD` or `YYYY/MM/DD`.
#' @param AccessionYear The accession year of the records.
#' @param ItemStatusDate The data when the records were last updated.
#' @param ItemStatusType The status of the record when last updated. Either `Existing` or `NotExisting`.
#' @param post_date Set the date that existing plants in the LC are classed as alive until . Default is the present day.

#'
#' @return data frame where each column corresponds to a date, each row an item in the collection. The (i,j)th value is a logical (TRUE/FALSE) corresponding whether the ith item existed at the jth date.
#' @export
#' @examples
#' exist_at_date(date = c('2003-06-25', '2010-03-07'),
#' AccessionYear = c('1956', '1988', '2005', '2018'),
#' ItemStatusDate = c('2004-01', '2010-03-05', '2022-04-08', '2022-08-19'),
#' ItemStatusType = c('NotExisting','NotExisting','Existing','Existing'))
exist_at_date <- function(date, AccessionYear, ItemStatusDate, ItemStatusType, post_date = as.character(Sys.Date())){

  date = as.Date(date)
  # Accession date.
  pre_date = rep(NA,length(AccessionYear)) ; pre_date_char = rep(NA,length(AccessionYear))
  wanted_index = which(AccessionYear > 1650 & AccessionYear <= as.numeric(format(Sys.Date(),'%Y')))
  pre_date[wanted_index] =as.Date(paste(AccessionYear[wanted_index], '01', '01', sep = "-"), "%Y-%m-%d")
  pre_date_char[wanted_index] = paste(AccessionYear[wanted_index], '01', '01', sep = "-")

  # use the current day, unless not existing then use the date of that entry.
  post_date = rep(post_date,length(AccessionYear))
  post_date[ItemStatusType == 'NotExisting'] = ItemStatusDate[ItemStatusType == 'NotExisting']
  # If only a year is given assume it occurs on the 31st of Dec
  index_string_length_4 = which(stringr::str_length(post_date) == 4)
  post_date[index_string_length_4] = paste( post_date[index_string_length_4], '12', '31', sep = "-")
  # If only a year and month is given assume the day is the 28th.
  index_string_length_7 = which(stringr::str_length(post_date) == 7)

  post_date[index_string_length_7] = paste( post_date[index_string_length_7], '28', sep = "-")
  dates = data.frame(pre = pre_date_char, post = post_date, pre_date = pre_date, post_date = as.Date(post_date, "%Y-%m-%d"))

  # Vector of whether the plant is existing on the date.
  out = matrix(NA, nrow = length(AccessionYear), ncol = length(date))
  for(i in 1:length(date)){
    existing_on_date = rep(FALSE, length(AccessionYear))
    existing_on_date[which(date[i] >= dates$pre_date & date[i] < dates$post_date)] = TRUE

    out[,i] = existing_on_date
  }

  out = data.frame(out)
  names(out) = date

  return(out)
}


#' Add transparency to colours
#'
#' @param col A vector of colours.
#' @param alpha numeric in `(0,1]`, the transparency wanted.
#'
#' @return a vector of colours with added transparency
#' @export
#' @examples
#' add_alpha(col = c("blue", rgb(100,20,30, max = 255), '#0A461E', viridis::viridis(1)),
#' alpha = 0.5)
add_alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, grDevices::col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
