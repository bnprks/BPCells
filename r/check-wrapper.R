
args <- commandArgs(TRUE)

# Disable timestamp check completely
Sys.setenv("_R_CHECK_FUTURE_FILE_TIMESTAMPS_"="false")
Sys.setenv("_R_CHECK_SYSTEM_CLOCK_"="0")

# Call R CMD check as normal
tools:::.check_packages(args)