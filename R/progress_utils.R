# =============================================================================
# progress_utils.R
# Console progress reporting: a weighted multi-step tracker with a live text
# progress bar, elapsed time, and an adaptive ETA for the whole pipeline.
#
# Two levels:
#   - progress_next()  -- call once per major pipeline step (parsing,
#     generating metabolites, building the workbook, etc). Prints a new
#     line: bar, step name, elapsed so far, ETA for the remaining steps.
#   - progress_tick()  -- for fine-grained progress *within* one expensive
#     step (e.g. metabolite-by-metabolite while building the Charge
#     Envelopes sheet). Overwrites the same console line via \r so it
#     doesn't spam the console with one line per metabolite.
#
# The ETA is adaptive, not a fixed prediction: step weights below are rough
# starting guesses (the workbook step -- specifically its Charge Envelopes
# sheet -- dominates real runtime, per earlier profiling), and the tracker
# recomputes "time per unit of weight" from actual elapsed time after each
# completed step, so the estimate self-corrects as the run progresses
# rather than trusting the initial guesses throughout.
# =============================================================================

.fmt_duration <- function(seconds) {
  if (is.na(seconds) || seconds < 0) return("--")
  if (seconds < 60) return(sprintf("%.0fs", seconds))
  m <- floor(seconds / 60); s <- round(seconds - m * 60)
  if (m < 60) return(sprintf("%dm %ds", m, s))
  h <- floor(m / 60); m <- m - h * 60
  sprintf("%dh %dm", h, m)
}

.progress_bar_str <- function(frac, width = 24) {
  frac <- max(0, min(1, frac))
  filled <- round(frac * width)
  paste0("[", strrep("#", filled), strrep("-", width - filled), "]")
}

# steps: named numeric vector of relative weights, e.g.
#   c("Parsing input" = 1, "Building workbook" = 15, ...)
# Order matters -- steps are expected to complete in the order given.
progress_tracker <- function(steps) {
  env <- new.env(parent = emptyenv())
  env$names <- names(steps)
  env$weights <- unname(as.numeric(steps))
  env$total_weight <- sum(env$weights)
  env$done_weight <- 0
  env$idx <- 0
  env$t0 <- Sys.time()
  env$step_t0 <- Sys.time()
  env
}

# Advance to the next named step; prints the step-level status line.
# Call this once per step, in order -- it automatically marks the
# *previous* step's weight as complete (if any) before reporting progress,
# so callers don't need a separate "mark done" call in between.
progress_next <- function(tracker) {
  if (tracker$idx >= 1 && tracker$idx <= length(tracker$weights)) {
    tracker$done_weight <- tracker$done_weight + tracker$weights[tracker$idx]
  }

  now <- Sys.time()
  elapsed_total <- as.numeric(now - tracker$t0, units = "secs")
  frac_done <- if (tracker$total_weight > 0) tracker$done_weight / tracker$total_weight else 0

  # Adaptive ETA: time-per-unit-weight from actual progress so far,
  # applied to the weight still remaining. First step has no data yet.
  eta_str <- "estimating..."
  if (tracker$done_weight > 0) {
    rate <- elapsed_total / tracker$done_weight
    remaining <- (tracker$total_weight - tracker$done_weight) * rate
    eta_str <- .fmt_duration(remaining)
  }

  tracker$idx <- tracker$idx + 1
  step_name <- if (tracker$idx <= length(tracker$names)) tracker$names[tracker$idx] else "Finishing"
  n_steps <- length(tracker$names)

  cat(sprintf("\n%s %3.0f%%  Step %d/%d: %s\n",
              .progress_bar_str(frac_done), frac_done * 100,
              tracker$idx, n_steps, step_name))
  cat(sprintf("  elapsed %s | ETA (remaining) %s\n",
              .fmt_duration(elapsed_total), eta_str))

  tracker$step_t0 <- now
  invisible(tracker)
}

# Fine-grained in-step progress (e.g. metabolite N/M while building a
# sheet). Overwrites the same console line. Call progress_tick_end() once
# to move to a fresh line when that sub-loop finishes.
progress_tick <- function(tracker, current, total, label = "") {
  now <- Sys.time()
  step_elapsed <- as.numeric(now - tracker$step_t0, units = "secs")
  frac <- if (total > 0) current / total else 1
  eta_step <- if (current > 0) .fmt_duration(step_elapsed / current * (total - current)) else "--"
  msg <- sprintf("  %s %3.0f%%  %s (%d/%d)  step elapsed %s  step ETA %s",
                 .progress_bar_str(frac, 16), frac * 100, label,
                 current, total, .fmt_duration(step_elapsed), eta_step)
  # Pad to overwrite any longer previous line, then carriage-return.
  cat("\r", msg, strrep(" ", max(0, 100 - nchar(msg))), sep = "")
  utils::flush.console()
  invisible(tracker)
}

progress_tick_end <- function() {
  cat("\n")
  invisible(NULL)
}

# Final summary line -- call once, after the last step's work completes
# (i.e. after the final progress_next() call and that step's actual work).
progress_finish <- function(tracker) {
  if (tracker$idx >= 1 && tracker$idx <= length(tracker$weights)) {
    tracker$done_weight <- tracker$done_weight + tracker$weights[tracker$idx]
    tracker$idx <- tracker$idx + 1  # prevent double-counting on repeat calls
  }
  elapsed_total <- as.numeric(Sys.time() - tracker$t0, units = "secs")
  cat(sprintf("\n%s 100%%  Done\n", .progress_bar_str(1)))
  cat(sprintf("  total elapsed: %s\n\n", .fmt_duration(elapsed_total)))
  invisible(tracker)
}
