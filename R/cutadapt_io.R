#' Takes the 'summary' text from a cutadapt log-file and extracts the number of
#' readpairs / basepairs that were input and which pass/fail various filters.
#'
#' Not exported
#'
#' @param        x             Text from the summary section of a cutadapt log.
#'
#' @importFrom   tibble        data_frame
#'

parse_cutadapt_summary <- function(x) {
  # Should x be newline-stripped / bareline stripped /
  #   have the "=== SUMMARY ===" line removed prior to calling this?
  if(missing(x)){
    stop("`x` should be defined in `parse_cutadapt_summary`")
  }

  fieldnames <- tibble::data_frame(
    expected = c(
      "Total read pairs processed", "Read 1 with adapter",
      "Read 2 with adapter",
      "Pairs that were too short", "Pairs written (passing filters)",
      "Total basepairs processed", "Read 1", "Read 2",
      "Total written (filtered)", "Read 1", "Read 2"
    ),
    output = c(
      "rp_input", "r1_with_adapter",
      "r2_with_adapter",
      "rp_too_short", "rp_output",
      "bp_input", "bp_input_r1", "bp_input_r2",
      "bp_output", "bp_output_r1", "bp_output_r2"
    )
  )

  # parse_numeric_fields(x, fieldnames)
  tibble::data_frame(
    rp_input = 114728,    r1_with_adapter = 2923, r2_with_adapter = 3833,
    rp_too_short = 57,    rp_output = 114671,
    bp_input = 17254556,  bp_input_r1 = 8628890,  bp_input_r2 = 8625666,
    bp_output = 17217852, bp_output_r1 = 8612728, bp_output_r2 = 8605124
  )
}
