# Requires: 
library(lubridate)
library(chron)
library(dplyr)

format_dates_and_age_verbatim <- function(metadata) {
  ##Format date of birth and calculate ages
  ##-- caluclate age 
  head(metadata$DOB)
  head(mdy(metadata$DOB))
  ##yeas before 1950 are in the wrong centory
  metadata$DOB_formatted<-mdy(metadata$DOB)
  #some of the DOBs when formatted  to YYYY from YY had assigned the wrong century
  as.Date(chron(format(as.Date(metadata$DOB_formatted, "%m/%d/%y"), "%m/%d/%y"))) 
  options(chron.year.expand =   #expanding the cutoff year
            function (y, cut.off = 12, century = c(1900, 2000), ...) {
              chron:::year.expand(y, cut.off = cut.off, century = century, ...)
            })
  metadata$DOB_formatted <- as.Date(chron(format(as.Date(metadata$DOB_formatted, "%m/%d/%y"), "%m/%d/%y"))) #re-formatting the years with the right century
  head(metadata$DOB_formatted)
  class(metadata$DOB_formatted)     
  
  ## Format the date collected to date format
  head(metadata$Date)
  metadata$Date_formatted <- mdy(metadata$Date) 
  head(metadata$Date_formatted, 4)
  
  ## calculate the age of the individual at the time of collection
  
  metadata <-metadata %>%
    mutate(
      age = year(Date_formatted) - year(DOB_formatted))
  
  head(metadata$age, 4)
  metadata
}
