#' Prepare inspectable phenology cross-validation splits
#'
#' Split complete season entries before fitting, keeping their hourly weather
#' intact. Years are inspection metadata, not unique observation identifiers.
#' Reserve independent validation entries by index or proportion. Preparation
#' does not fit or predict any model.
#'
#' @param seasons A seasonlist for single/stage models, or nested seasonlists
#'   for combined models. Each season is a complete hourly weather data frame.
#' @param observed Observations in the format used by phenology_rss(),
#'   phenology_rss_combined() or phenology_rss_stages(). NULL slots mark missing
#'   observations; they are preserved when subsetting.
#' @param years Phenological year per season, or one vector per cultivar for
#'   combined data. NULL uses the maximum Year in each weather frame. Supply
#'   explicit years when this does not identify your phenological season.
#'   Repeated years at different locations are allowed and are split separately.
#' @param layout One of "single", "combined", "stages"; describes the data,
#'   independently of the model to be fitted.
#' @param v Number of CV folds, at least two and at most the number of remaining
#'   entry groups (one per season entry by default).
#' @param repeats Number of independent random fold assignments. Optimizer
#'   restarts are a separate argument of fit_phenology_cv().
#' @param validation_indices Optional one-based indices of season entries to
#'   withhold. For single/stage data, a numeric vector indexing seasons. For
#'   combined data, a list of numeric index vectors, one per member in input
#'   order (use integer() for members with no validation entries). All stages of
#'   a selected common season are withheld together. Indices select exactly
#'   these entries, without also withholding other entries from the same year.
#' @param validation_prop Alternatively, randomly reserve this proportion of
#'   observed season entries, rounded up; 0.2 means 20 percent. NULL-only entries
#'   are not sampled. A common season with at least one observed stage counts
#'   once for stage data; combined entries are counted across all members.
#'   Cannot be combined with validation_indices. NULL reserves no entries.
#' @param groups Optional CV grouping labels per season, or one vector per member
#'   for combined data. Entries with the same label stay in the same CV fold.
#'   NULL splits entries independently, except that stages of a common season
#'   always stay together. Use location/year labels to group related observations,
#'   or years for whole-year CV. Groups affect CV folds only; validation_indices
#'   and validation_prop select entries independently of these labels.
#' @param seed Non-negative integer seed, or NULL to use the caller's RNG.
#'   A supplied seed preserves the caller's RNG state.
#' @return A phenology_cv_data list with three elements: data (original seasons,
#'   observed dates and a cases table of entry identifiers), assignments (case_id,
#'   member, season, year, set, repeat_id and fold), and settings (preparation
#'   options). Dates and groups are not repeated in assignments. Training and
#'   assessment indices are derived from this table when fitting, not stored
#'   separately. Printing shows counts and six preview rows. Inspect the full
#'   assignments table before passing this object to fit_phenology_cv().
#' @examples
#' weather <- data.frame(Temp = rep(8, 480), Year = 2008,
#'                       JDay = rep(50:69, each = 24))
#' seasons <- rep(list(weather), 6)
#' prepared <- prepare_phenology_cv(seasons, rep(60, 6), years = 2001:2006,
#'                                  v = 2, repeats = 2, validation_indices = 6)
#' prepared$assignments
#' subset(prepared$assignments, set == "cv" & repeat_id == 1 & fold == 1)
#' @md
#' @export
prepare_phenology_cv <- function(seasons, observed, years = NULL,
                                 layout = c("single", "combined", "stages"),
                                 v = 5, repeats = 1, validation_indices = NULL,
                                 validation_prop = NULL, seed = 12345, groups = NULL) {
  layout <- match.arg(layout)
  .fit_positive_scalar(v, "v", integer = TRUE)
  if (v < 2) stop("v must be at least two.", call. = FALSE)
  .fit_positive_scalar(repeats, "repeats", integer = TRUE)
  if (!is.null(seed)) .fit_positive_scalar(seed, "seed", integer = TRUE, zero = TRUE)
  if (layout == "combined") {
    if (!is.list(seasons) || is.data.frame(seasons) || !length(seasons) ||
        !all(vapply(seasons, .is_pheno_seasonlist, logical(1))))
      stop("combined seasons must contain a non-empty seasonlist per member.", call. = FALSE)
    weather <- seasons
  } else {
    if (!.is_pheno_seasonlist(seasons))
      stop("seasons must be a non-empty seasonlist of weather data frames.", call. = FALSE)
    weather <- if (layout == "single") list(seasons) else {
      if (!is.list(observed) || is.data.frame(observed) || !length(observed))
        stop("stages observed must contain one observation vector/list per stage.", call. = FALSE)
      rep(list(seasons), length(observed))
    }
  }
  observations <- if (layout == "single") list(observed) else observed
  if (!is.list(observations) || is.data.frame(observations) ||
      length(observations) != length(weather))
    stop("observed must contain one observation vector/list per member.", call. = FALSE)
  year_vectors <- if (is.null(years)) lapply(weather, function(w) {
    vapply(w, function(s) {
      if (!is.numeric(s$Year) || !length(s$Year) || any(!is.finite(s$Year)))
        stop("Supply years explicitly or include finite Year values in each season.", call. = FALSE)
      max(s$Year)
    }, numeric(1))
  }) else if (layout == "combined") years else rep(list(years), length(weather))
  if (!is.list(year_vectors) || length(year_vectors) != length(weather))
    stop("years must match the seasonlist layout.", call. = FALSE)
  labels <- if (layout == "single") "model" else
    if (layout == "stages") names(observed) else names(seasons)
  explicit_member_names <- !is.null(labels)
  if (is.null(labels)) labels <- paste0("member", seq_along(weather))
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
    stop("Member names must be unique and non-empty.", call. = FALSE)
  cases <- do.call(rbind, lapply(seq_along(weather), function(i) {
    n <- length(weather[[i]])
    y <- year_vectors[[i]]
    if (!is.numeric(y) || !is.null(dim(y)) || length(y) != n ||
        any(!is.finite(y)) || any(y != floor(y)))
      stop("years must contain one finite integer year per season.", call. = FALSE)
    dates <- .phenology_observations(observations[[i]], n)
    season_names <- names(weather[[i]])
    if (is.null(season_names)) season_names <- as.character(seq_len(n))
    data.frame(member = i, member_name = labels[i], season = seq_len(n),
               season_name = season_names, year = unname(y), observed = unname(dates),
               stringsAsFactors = FALSE)
  }))
  rownames(cases) <- NULL
  cases$case_id <- seq_len(nrow(cases))
  # A stage season is one sampling unit, even though it contains several dates.
  cases$entry_id <- if (layout == "stages") cases$season else cases$case_id
  if (is.null(groups)) {
    cases$group <- paste0("entry", cases$entry_id)
  } else {
    group_vectors <- if (layout == "combined") groups else rep(list(groups), length(weather))
    if (!is.list(group_vectors) || length(group_vectors) != length(weather))
      stop("groups must match the seasonlist layout.", call. = FALSE)
    cases$group <- unlist(lapply(seq_along(weather), function(i) {
      g <- group_vectors[[i]]
      if ((!is.character(g) && !is.numeric(g) && !is.factor(g)) ||
          !is.null(dim(g)) || length(g) != length(weather[[i]]) || anyNA(g) ||
          (is.numeric(g) && any(!is.finite(g))) || any(!nzchar(as.character(g))))
        stop("groups must contain one non-missing, non-empty label per season.", call. = FALSE)
      as.character(g)
    }), use.names = FALSE)
  }
  if (!is.null(validation_indices) && !is.null(validation_prop))
    stop("Choose validation_indices or validation_prop, not both.", call. = FALSE)
  reserved_cases <- .cv_validation_cases(validation_indices, cases, layout, lengths(weather), labels)
  if (!is.null(validation_prop) &&
      (!is.numeric(validation_prop) || length(validation_prop) != 1L ||
       !is.null(dim(validation_prop)) ||
       !is.finite(validation_prop) || validation_prop <= 0 || validation_prop >= 1))
    stop("validation_prop must be a number strictly between zero and one.", call. = FALSE)
  result <- .cv_with_seed(seed, function() {
    if (!is.null(validation_prop)) {
      eligible <- unique(cases$entry_id[!is.na(cases$observed)])
      if (!length(eligible)) stop("There must be observed entries to reserve validation data.", call. = FALSE)
      reserved <- eligible[sample.int(length(eligible), ceiling(length(eligible) * validation_prop))]
      reserved_cases <- which(cases$entry_id %in% reserved)
    }
    cv_indices <- setdiff(cases$case_id, reserved_cases)
    cv_groups <- unique(cases$group[cv_indices])
    if (length(cv_groups) < v)
      stop("There must be at least v CV entry groups after reserving validation entries.", call. = FALSE)
    assignments <- list()
    for (r in seq_len(repeats)) {
      fold <- rep(seq_len(v), length.out = length(cv_groups))
      fold <- fold[sample.int(length(fold))]
      assignments[[r]] <- cbind(cases[cv_indices, c("case_id", "member", "season", "year"), drop = FALSE], set = "cv",
        repeat_id = r, fold = fold[match(cases$group[cv_indices], cv_groups)])
      for (f in seq_len(v)) {
        assessment <- cv_indices[cases$group[cv_indices] %in% cv_groups[fold == f]]
        train <- cv_indices[cases$group[cv_indices] %in% cv_groups[fold != f]]
        if (!any(!is.na(cases$observed[train])) ||
            !any(!is.na(cases$observed[assessment])))
          stop("Every CV train/assessment split must contain an observed date.", call. = FALSE)
      }
    }
    if (length(reserved_cases)) assignments[[length(assignments) + 1L]] <-
      cbind(cases[reserved_cases, c("case_id", "member", "season", "year"), drop = FALSE], set = "validation",
            repeat_id = NA_integer_, fold = NA_integer_)
    assignments <- do.call(rbind, assignments)
    rownames(assignments) <- NULL
    assignments
  })
  structure(list(data = list(seasons = seasons, observed = observed,
                             cases = cases[, setdiff(names(cases), "observed"), drop = FALSE]),
                 assignments = result,
                 settings = list(layout = layout, v = v, repeats = repeats, seed = seed,
                   validation_indices = validation_indices, validation_prop = validation_prop,
                   grouped = !is.null(groups), member_names = labels,
                   explicit_member_names = explicit_member_names)),
            class = "phenology_cv_data")
}

.cv_validation_cases <- function(indices, cases, layout, counts, labels) {
  if (is.null(indices)) return(integer())
  valid <- function(x, n) is.numeric(x) && is.null(dim(x)) &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (layout != "combined") {
    if (!valid(indices, counts[1]))
      stop("validation_indices must contain distinct one-based season indices within range.", call. = FALSE)
    return(which(cases$season %in% indices))
  }
  if (!is.list(indices) || is.data.frame(indices) || length(indices) != length(counts) ||
      (!is.null(names(indices)) && !identical(names(indices), labels)))
    stop("combined validation_indices must be a list in member order, one index vector per member.", call. = FALSE)
  for (i in seq_along(counts)) {
    if (!valid(indices[[i]], counts[i]))
      stop("validation_indices must contain distinct one-based season indices within range for each member.", call. = FALSE)
  }
  which(vapply(seq_len(nrow(cases)), function(k)
    cases$season[k] %in% indices[[cases$member[k]]], logical(1)))
}

.cv_with_seed <- function(seed, fn) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    })
    set.seed(seed)
  }
  fn()
}

.cv_subset <- function(data, indices) {
  layout <- data$settings$layout
  cases <- data$data$cases
  if (layout == "single") {
    k <- cases$season[indices]
    return(list(seasons = data$data$seasons[k], observed = data$data$observed[k]))
  }
  n_members <- length(data$settings$member_names)
  by_member <- lapply(seq_len(n_members), function(i)
    cases$season[indices[cases$member[indices] == i]])
  observed <- lapply(seq_len(n_members), function(i) data$data$observed[[i]][by_member[[i]]])
  names(observed) <- names(data$data$observed)
  if (layout == "stages")
    return(list(seasons = data$data$seasons[by_member[[1]]], observed = observed))
  seasons <- lapply(seq_len(n_members), function(i) data$data$seasons[[i]][by_member[[i]]])
  names(seasons) <- names(data$data$seasons)
  list(seasons = seasons, observed = observed)
}

.cv_cases <- function(data) {
  cases <- data$data$cases
  dates <- if (data$settings$layout == "single") list(data$data$observed) else data$data$observed
  cases$observed <- unlist(lapply(seq_along(dates), function(i)
    .phenology_observations(dates[[i]], sum(cases$member == i))), use.names = FALSE)
  cases
}

.cv_validation_indices <- function(data) {
  unique(data$assignments$case_id[data$assignments$set == "validation"])
}

.cv_splits <- function(data) {
  assignments <- data$assignments
  cv <- setdiff(data$data$cases$case_id, .cv_validation_indices(data))
  splits <- list()
  for (r in seq_len(data$settings$repeats)) {
    rows <- assignments[assignments$set == "cv" & !is.na(assignments$repeat_id) &
                          assignments$repeat_id == r, , drop = FALSE]
    for (f in seq_len(data$settings$v)) {
      id <- paste0("repeat", r, "_fold", f)
      assessment <- cv[cv %in% rows$case_id[rows$fold == f]]
      splits[[id]] <- list(id = id, repeat_id = r, fold = f,
                           train = setdiff(cv, assessment), assessment = assessment)
    }
  }
  splits
}

#' @export
print.phenology_cv_data <- function(x, ...) {
  cat("Prepared phenology CV:", length(unique(x$data$cases$entry_id)), "entries;",
      x$settings$v, "folds x", x$settings$repeats, "repeat(s);",
      length(unique(x$data$cases$entry_id[.cv_validation_indices(x)])), "independent validation entries\n")
  print(utils::head(x$assignments, 6L), row.names = FALSE, ...)
  if (nrow(x$assignments) > 6L)
    cat("...", nrow(x$assignments) - 6L, "more rows. Inspect $assignments for the full table.\n")
  invisible(x)
}
