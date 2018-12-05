.is_single_string <- function(.x) {
  is.character(.x) && length(.x) == 1
}
