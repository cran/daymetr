% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_daymet_batch.r
\name{download_daymet_batch}
\alias{download_daymet_batch}
\title{This function downloads 'Daymet' data for several single pixel
location, as specified by a batch file.}
\usage{
download_daymet_batch(
  file_location = NULL,
  start = 1980,
  end = as.numeric(format(Sys.time(), "\%Y")) - 1,
  internal = TRUE,
  force = FALSE,
  silent = FALSE,
  path = tempdir(),
  simplify = FALSE
)
}
\arguments{
\item{file_location}{file with several site locations and coordinates
in a comma delimited format: site, latitude, longitude}

\item{start}{start of the range of years over which to download data}

\item{end}{end of the range of years over which to download data}

\item{internal}{assign or FALSE, load data into workspace or save to disc}

\item{force}{\code{TRUE} or \code{FALSE} (default),
override the conservative end year setting}

\item{silent}{suppress the verbose output (default = FALSE)}

\item{path}{set path where to save the data
if internal = FALSE (default = tempdir())}

\item{simplify}{output data to a tibble, logical \code{FALSE}
or \code{TRUE} (default = \code{TRUE})}
}
\value{
Daymet data for point locations as a nested list or
data written to csv files
}
\description{
This function downloads 'Daymet' data for several single pixel
location, as specified by a batch file.
}
\examples{

\dontrun{
# The download_daymet_batch() routine is a wrapper around
# the download_daymet() function. It queries a file with
# coordinates to easily download a large batch of daymet
# pixel locations. When internal = TRUE, the data is stored
# in a structured list in an R variable. If FALSE, the data
# is written to disk.

# create demo locations (two sites)
locations <- data.frame(site = c("site1", "site2"),
                      lat = rep(36.0133, 2),
                      lon = rep(-84.2625, 2))

# write data to csv file
write.table(locations, paste0(tempdir(),"/locations.csv"),
           sep = ",",
           col.names = TRUE,
           row.names = FALSE,
           quote = FALSE)

# download data, will return nested list of daymet data
df_batch <- download_daymet_batch(file_location = paste0(tempdir(),
                                                        "/locations.csv"),
                                    start = 1980,
                                    end = 1980,
                                    internal = TRUE,
                                    silent = TRUE)

# For other practical examples consult the included
# vignette. 
}
}
