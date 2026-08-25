# Packages
library(terra)
library(modifiedmk)


# Load short-season NDVI stack
short_ndvi_stack <- rast(
  "~/SeasonalTrendsEdited/data/ndvi_short_season_stack_1997_2024.tif"
)

# Pixel-wise Hamed–Rao modified Mann–Kendall function
mk_fun <- function(x) {

  # Return NA if the entire pixel time series is missing
  if (all(is.na(x))) {
    return(c(
      slope = NA_real_,
      pvalue = NA_real_
    ))
  }

  # Remove missing observations
  x <- x[!is.na(x)]

  # Require a minimum number of valid observations
  if (length(x) < 10) {
    return(c(
      slope = NA_real_,
      pvalue = NA_real_
    ))
  }


  # Apply Hamed–Rao modified Mann–Kendall test. mmkh3lag accounts for serial autocorrelation. Here using the first three lags.
  mk <- try(
    suppressWarnings(
      modifiedmk::mmkh3lag(x)
    ),
    silent = TRUE
  )


  # Return NA if the test fails for a pixel
  if (inherits(mk, "try-error")) {
    return(c(
      slope = NA_real_,
      pvalue = NA_real_
    ))
  }


  # Extract outputs from mmkh3lag
  # [[7]] = Sen's slope
  # [[2]] = corrected p-value
  slope <- as.numeric(
    mk[[7]]
  )

  pvalue <- as.numeric(
    mk[[2]]
  )


  # Convert non-finite outputs such as NaN/Inf to NA
  if (!is.finite(slope)) {
    slope <- NA_real_
  }

  if (!is.finite(pvalue)) {
    pvalue <- NA_real_
  }


  # Return Sen's slope and corrected p-value
  return(c(
    slope = slope,
    pvalue = pvalue
  ))
}


# Apply modified Mann–Kendall to every pixel
short_mk <- app(
  short_ndvi_stack,
  fun = mk_fun,
  filename = "~/SeasonalTrendsEdited/output/short_season_MK_results.tif",
  overwrite = TRUE
)


# Assign descriptive layer names
names(short_mk) <- c(
  "Sen_slope",
  "p_value"
)