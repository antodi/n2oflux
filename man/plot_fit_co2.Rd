% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fit_co2.R
\name{plot_fit_co2}
\alias{plot_fit_co2}
\title{Show results of linear and non-linear fits for CO2 data, along with deadband and outliers when used.}
\usage{
plot_fit_co2(id, data, deadband_c = 30, stop_time_ag = 120, offset_k = "json")
}
\arguments{
\item{id}{The ID (written as date.LABEL) of measurements to plot with CO2 data. Several measurements can be processed at the same time.}

\item{data}{The data containing the measurements seconds by seconds, given by the process_json_files function.}

\item{deadband_c}{The desired value (sec) for the deadband. CO2 observations before that time will be discarded.}

\item{stop_time_ag}{The duration of observation (sec).}

\item{offset_k}{The value of the offset (cm) for the collar. If set to "json", the offset will be retrieved from the json file.}
}
\value{
A list containing the plots for all IDs
}
\description{
Show results of linear and non-linear fits for CO2 data, along with deadband and outliers when used.
}
