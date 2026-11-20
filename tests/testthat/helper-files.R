# A new folder that is removed when the calling test finishes.
temp_dir <- function(env = parent.frame()) {
  d <- tempfile("joinpoint-test")
  dir.create(d)
  do.call(on.exit, list(bquote(unlink(.(d), recursive = TRUE)), add = TRUE),
          envir = env)
  d
}
