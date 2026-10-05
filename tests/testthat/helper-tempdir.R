# A fresh directory under the session temp dir (removed when R exits).
test_tempdir <- function() {
  path <- tempfile("rsemflow-test-")
  dir.create(path)
  path
}
