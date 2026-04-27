#' Download and Cache Stirling State Data
#'
#' Downloads the pre-computed long-double state blocks from GitHub and
#' saves them to the user's local data directory. Once downloaded,
#' \code{logStirling2} will automatically detect and use these states
#' for accelerated calculations.
#'
#' @param force Logical; if \code{TRUE}, re-downloads the data even if it
#'   already exists locally.
#'
#' @return Invisible \code{TRUE} on success.
#' @export
get_state_data <- function(force = FALSE) {
  # 1. Define paths
  data_dir <- tools::R_user_dir("logStirling2", "data")
  dest_file <- file.path(data_dir, "logStirling2_cache_v01.rds")

  # 2. Check if already exists
  if (file.exists(dest_file) && !force) {
    message("State data already exists at: ", dest_file)
    return(invisible(TRUE))
  }

  # 3. Create directory if needed
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }

  # 4. Download from your GitHub Repository
  # Replace with your actual GitHub username/repo
  url <- "https://github.com/jblood94/logStirling2/releases/download/v0.1.0/logStirling2_cache_v01.rds"

  message("Downloading state data (approx 10MB)...")
  tryCatch({
    utils::download.file(url, destfile = dest_file, mode = "wb")
  }, error = function(e) {
    stop("Download failed. Please check your internet connection or the URL.")
  })

  # 5. Hot-load into the session cache (.state_env)
  # This ensures the current session sees the new data immediately
  if (.Machine$sizeof.longdouble == 16) {
    temp_states <- readRDS(dest_file)
    .state_env$logS_states <- temp_states
    message("Data successfully downloaded and loaded into memory.")
  } else {
    warning("Data downloaded, but your architecture does not support 16-byte
            long doubles. Calculation speed will not be affected.")
  }

  invisible(TRUE)
}
