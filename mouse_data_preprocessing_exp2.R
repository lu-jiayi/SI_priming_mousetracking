if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  source_path <- rstudioapi::getSourceEditorContext()$path
  if (nzchar(source_path)) {
    setwd(dirname(source_path))
  }
}

library(tidyverse)

first_value <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

min_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

max_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

decode_mouse_trace <- function(trace) {
  empty_trace <- tibble(
    sample_index = integer(),
    time_ms = double(),
    x_px = double(),
    y_px = double(),
    trace_width_px = double(),
    trace_height_px = double()
  )
  
  if (is.na(trace) || trace == "") {
    return(empty_trace)
  }
  
  header_pattern <- "^x(-?\\d+)y(-?\\d+)w(\\d+)h(\\d+)"
  header <- str_match(trace, header_pattern)
  
  if (any(is.na(header[1, 2:5]))) {
    warning("Could not decode a MouseTracker trace: ", str_trunc(trace, 80))
    return(empty_trace)
  }
  
  x_start <- as.double(header[1, 2])
  y_start <- as.double(header[1, 3])
  trace_width <- as.double(header[1, 4])
  trace_height <- as.double(header[1, 5])
  
  changes_text <- str_remove(trace, header_pattern)
  changes <- str_match_all(
    changes_text,
    "t(\\d+)([+-]\\d+)([+-]\\d+)"
  )[[1]]
  
  if (nrow(changes) == 0) {
    return(tibble(
      sample_index = 0L,
      time_ms = 0,
      x_px = x_start,
      y_px = y_start,
      trace_width_px = trace_width,
      trace_height_px = trace_height
    ))
  }
  
  time_change <- as.double(changes[, 2])
  x_change <- as.double(changes[, 3])
  y_change <- as.double(changes[, 4])
  
  tibble(
    sample_index = seq.int(0L, nrow(changes)),
    time_ms = c(0, cumsum(time_change)),
    x_px = c(x_start, x_start + cumsum(x_change)),
    y_px = c(y_start, y_start + cumsum(y_change)),
    trace_width_px = trace_width,
    trace_height_px = trace_height
  )
}

# -----------------------------------------------------------------------------
# Input and output locations
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "raw_data_exp2.csv"
output_dir <- if (length(args) >= 2) args[[2]] else "cleaned_data"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read the mixed-width PCIbex file
# -----------------------------------------------------------------------------

# The follow-up experiment has:
#   * ordinary PCIbex rows with 13 fields; and
#   * experimental rows with 28 fields.
#
# We remove only full comment lines before parsing. Using comment.char = "#"
# would be unsafe here because the new partner-color fields contain hex values
# beginning with # (e.g., #0072B2).
raw_lines <- readLines(input_file, warn = FALSE, encoding = "UTF-8")
if (length(raw_lines) > 0) {
  raw_lines[[1]] <- sub("^\ufeff", "", raw_lines[[1]])
}
data_lines <- raw_lines[!str_detect(raw_lines, "^\\s*#")]

if (length(data_lines) == 0) {
  stop("No PCIbex data rows were found in ", input_file, ".")
}

pcibex_names <- c(
  "reception_time", "participant_hash", "controller", "item_order",
  "inner_element", "label", "latin_square_group", "element_type",
  "element_name", "parameter", "value", "event_time",
  paste0("field_", 13:28)
)

raw_connection <- textConnection(data_lines)
raw_import <- utils::read.csv(
  raw_connection,
  header = FALSE,
  col.names = pcibex_names,
  colClasses = "character",
  comment.char = "",
  na.strings = c("", "NULL", "NA"),
  quote = "\"",
  fill = TRUE,
  strip.white = FALSE,
  check.names = FALSE
) |>
  as_tibble()
close(raw_connection)

if (ncol(raw_import) != 28) {
  stop("Expected 28 columns after padding, but read ", ncol(raw_import), ".")
}

# In 28-field rows, fields 13--27 contain the follow-up trial metadata and
# field 28 contains comments. In ordinary 13-field rows, field 13 is comments.
raw_events <- raw_import |>
  mutate(
    extended_layout = !is.na(field_14),
    condition = if_else(extended_layout, field_13, NA_character_),
    trial_type = if_else(extended_layout, field_14, NA_character_),
    picture_l = if_else(extended_layout, field_15, NA_character_),
    picture_r = if_else(extended_layout, field_16, NA_character_),
    sentence = if_else(extended_layout, field_17, NA_character_),
    partner = if_else(extended_layout, field_18, NA_character_),
    partner_color = if_else(extended_layout, field_19, NA_character_),
    strong_prime_speaker = if_else(extended_layout, field_20, NA_character_),
    weak_prime_speaker = if_else(extended_layout, field_21, NA_character_),
    other_condition_speaker = if_else(extended_layout, field_22, NA_character_),
    strong_prime_color = if_else(extended_layout, field_23, NA_character_),
    weak_prime_color = if_else(extended_layout, field_24, NA_character_),
    other_condition_color = if_else(extended_layout, field_25, NA_character_),
    window_width = if_else(extended_layout, field_26, NA_character_),
    window_height = if_else(extended_layout, field_27, NA_character_),
    comments = if_else(extended_layout, field_28, field_13)
  ) |>
  transmute(
    reception_time = parse_double(
      reception_time,
      na = c("", "NULL", "NA", "Never")
    ),
    participant_hash,
    controller,
    item_order = parse_integer(item_order),
    inner_element = parse_integer(inner_element),
    label,
    latin_square_group = na_if(latin_square_group, "NULL"),
    element_type,
    element_name,
    parameter,
    value,
    event_time = parse_double(
      event_time,
      na = c("", "NULL", "NA", "Never")
    ),
    condition,
    trial_type,
    picture_l,
    picture_r,
    sentence,
    partner,
    partner_color,
    strong_prime_speaker,
    weak_prime_speaker,
    other_condition_speaker,
    strong_prime_color,
    weak_prime_color,
    other_condition_color,
    window_width = parse_double(window_width),
    window_height = parse_double(window_height),
    comments
  )

# -----------------------------------------------------------------------------
# Participant-level information
# -----------------------------------------------------------------------------

demographic_events <- raw_events |>
  filter(label == "intro-1", element_type == "Html")

participants <- demographic_events |>
  group_by(participant_hash) |>
  summarise(
    prolific_id = first_value(value[parameter == "ID"]),
    consent = first_value(value[parameter == "consent"]),
    age = first_value(
      value[
        parameter == "age" &
          str_detect(coalesce(value, ""), "^\\s*\\d+\\s*$")
      ]
    ),
    gender = first_value(
      value[
        parameter == "age" &
          !str_detect(coalesce(value, ""), "^\\s*\\d+\\s*$")
      ]
    ),
    native_language = first_value(value[parameter == "natlang"]),
    dominant_language = first_value(value[parameter == "domlang"]),
    parent_language = first_value(value[parameter == "parentlang"]),
    other_languages = first_value(value[parameter == "otherlang"]),
    state = first_value(value[parameter == "state"]),
    .groups = "drop"
  ) |>
  mutate(
    participant_id = coalesce(prolific_id, participant_hash),
    age = parse_integer(age),
    gender = str_to_lower(str_squish(gender)),
    across(
      c(native_language, dominant_language, parent_language),
      ~ str_to_lower(str_squish(.x))
    ),
    across(c(other_languages, state), str_squish)
  ) |>
  select(
    participant_id, prolific_id, participant_hash, consent, age, gender,
    native_language, dominant_language, parent_language, other_languages, state
  ) |>
  arrange(participant_id)

# -----------------------------------------------------------------------------
# Trial-level data
# -----------------------------------------------------------------------------

speaker_fields <- c(
  "partner", "partner_color", "strong_prime_speaker", "weak_prime_speaker",
  "other_condition_speaker", "strong_prime_color", "weak_prime_color",
  "other_condition_color"
)

experiment_events <- raw_events |>
  filter(label == "actual") |>
  left_join(
    participants |> select(participant_hash, participant_id, prolific_id),
    by = "participant_hash"
  ) |>
  mutate(
    participant_id = coalesce(participant_id, participant_hash),
    received_at_utc = as.POSIXct(
      reception_time,
      origin = "1970-01-01",
      tz = "UTC"
    )
  ) |>
  select(
    participant_id, prolific_id, participant_hash, received_at_utc,
    item_order, inner_element, condition, trial_type, picture_l, picture_r,
    sentence, all_of(speaker_fields), window_width, window_height,
    element_type, element_name, parameter, value, event_time, comments
  ) |>
  arrange(participant_id, item_order, event_time)

trials <- experiment_events |>
  group_by(participant_id, prolific_id, participant_hash, item_order) |>
  summarise(
    condition = first_value(condition),
    trial_type = first_value(trial_type),
    picture_l = first_value(picture_l),
    picture_r = first_value(picture_r),
    sentence = first_value(sentence),
    partner = first_value(partner),
    partner_color = first_value(partner_color),
    strong_prime_speaker = first_value(strong_prime_speaker),
    weak_prime_speaker = first_value(weak_prime_speaker),
    other_condition_speaker = first_value(other_condition_speaker),
    strong_prime_color = first_value(strong_prime_color),
    weak_prime_color = first_value(weak_prime_color),
    other_condition_color = first_value(other_condition_color),
    mouse_trajectory = first_value(
      value[element_type == "MouseTracker" & parameter == "Move"]
    ),
    window_width = as.double(first_value(as.character(window_width))),
    window_height = as.double(first_value(as.character(window_height))),
    trial_start_time = min_value(
      event_time[
        element_type == "PennController" & parameter == "_Trial_" &
          value == "Start"
      ]
    ),
    tracker_start_time = min_value(
      event_time[element_type == "MouseTracker" & parameter == "Move"]
    ),
    selection_time = min_value(
      event_time[element_type == "Selector" & parameter == "Selection"]
    ),
    trial_end_time = max_value(
      event_time[
        element_type == "PennController" & parameter == "_Trial_" &
          value == "End"
      ]
    ),
    choice = first_value(
      value[element_type == "Selector" & parameter == "Selection"]
    ),
    selector_order = first_value(
      comments[
        element_type == "Selector" & parameter == "Selection"
      ]
    ),
    too_slow = as.integer(first_value(
      value[element_type == "Var" & element_name == "too_slow"]
    )),
    n_clicks = sum(
      element_type == "MouseTracker" & parameter == "Click",
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  group_by(participant_id, participant_hash) |>
  arrange(item_order, .by_group = TRUE) |>
  mutate(trial_number = row_number()) |>
  ungroup() |>
  mutate(
    left_choice = str_extract(selector_order, "^[^;]+"),
    right_choice = str_extract(selector_order, "[^;]+$"),
    
    # Physical screen location, used to align mouse trajectories.
    selected_side = case_when(
      choice == left_choice ~ "left",
      choice == right_choice ~ "right",
      TRUE ~ NA_character_
    ),
    
    # Picture identity, used for response coding. choiceL always means
    # picture_l and choiceR always means picture_r, regardless of screen order.
    selected_picture = case_when(
      choice == "choiceL" ~ picture_l,
      choice == "choiceR" ~ picture_r,
      TRUE ~ NA_character_
    ),
    
    trial_duration_ms = trial_end_time - trial_start_time,
    response_time_from_trial_start_ms = selection_time - trial_start_time,
    response_time_from_tracker_start_ms = selection_time - tracker_start_time,
    too_slow = coalesce(too_slow, 0L),
    valid_response =
      too_slow == 0L &
      !is.na(selected_side) &
      !is.na(selected_picture) &
      !is.na(tracker_start_time) &
      !is.na(selection_time)
  ) |>
  select(
    participant_id, prolific_id, participant_hash, trial_number, item_order,
    condition, trial_type, sentence, picture_l, picture_r,
    all_of(speaker_fields), choice, selector_order, selected_side,
    selected_picture, mouse_trajectory, window_width, window_height,
    trial_start_time, tracker_start_time, selection_time, trial_end_time,
    trial_duration_ms, response_time_from_trial_start_ms,
    response_time_from_tracker_start_ms, n_clicks, too_slow, valid_response
  ) |>
  arrange(participant_id, trial_number)

# -----------------------------------------------------------------------------
# Decode mouse movements into one row per sample
# -----------------------------------------------------------------------------

movement_rows <- experiment_events |>
  filter(element_type == "MouseTracker", parameter == "Move") |>
  transmute(
    participant_id, prolific_id, participant_hash, item_order,
    tracker_start_time = event_time,
    trace = value,
    decoded = map(trace, decode_mouse_trace)
  )

mouse_trajectories <- movement_rows |>
  select(-trace) |>
  unnest(decoded) |>
  left_join(
    trials |>
      select(
        participant_hash, item_order, trial_number, condition, trial_type,
        sentence, picture_l, picture_r, all_of(speaker_fields), selected_side,
        selected_picture, trial_start_time, selection_time, valid_response
      ),
    by = c("participant_hash", "item_order")
  ) |>
  mutate(
    time_from_trial_start_ms =
      tracker_start_time + time_ms - trial_start_time,
    time_before_response_ms =
      tracker_start_time + time_ms - selection_time,
    x_prop = x_px / trace_width_px,
    y_prop = y_px / trace_height_px,
    x_centered = 2 * (x_prop - 0.5),
    y_centered = 2 * (y_prop - 0.5),
    # Align horizontal movement so the selected physical side is positive.
    x_toward_choice = case_when(
      selected_side == "right" ~ x_centered,
      selected_side == "left" ~ -x_centered,
      TRUE ~ NA_real_
    )
  ) |>
  select(
    participant_id, prolific_id, participant_hash, trial_number, item_order,
    condition, trial_type, sentence, picture_l, picture_r,
    all_of(speaker_fields), selected_side, selected_picture, valid_response,
    sample_index, time_ms, time_from_trial_start_ms, time_before_response_ms,
    x_px, y_px, trace_width_px, trace_height_px, x_prop, y_prop, x_centered,
    y_centered, x_toward_choice
  ) |>
  arrange(participant_id, trial_number, sample_index)

# -----------------------------------------------------------------------------
# Keep mouse clicks in a separate long-format table
# -----------------------------------------------------------------------------

clicks <- experiment_events |>
  filter(element_type == "MouseTracker", parameter == "Click") |>
  transmute(
    participant_id, prolific_id, participant_hash, item_order,
    click_time = event_time,
    click_location = value
  ) |>
  extract(
    click_location,
    into = c("click_x_px", "click_y_px"),
    regex = "^(-?\\d+):(-?\\d+)$",
    remove = TRUE,
    convert = TRUE
  ) |>
  left_join(
    trials |>
      select(
        participant_hash, item_order, trial_number, condition, trial_type,
        all_of(speaker_fields), trial_start_time, tracker_start_time,
        selection_time, valid_response
      ),
    by = c("participant_hash", "item_order")
  ) |>
  group_by(participant_id, item_order) |>
  arrange(click_time, .by_group = TRUE) |>
  mutate(
    click_index = row_number(),
    click_time_from_trial_start_ms = click_time - trial_start_time,
    click_time_from_tracker_start_ms = click_time - tracker_start_time,
    click_time_before_response_ms = click_time - selection_time
  ) |>
  ungroup() |>
  select(
    participant_id, prolific_id, participant_hash, trial_number, item_order,
    condition, trial_type, all_of(speaker_fields), valid_response, click_index,
    click_time, click_time_from_trial_start_ms,
    click_time_from_tracker_start_ms, click_time_before_response_ms,
    click_x_px, click_y_px
  ) |>
  arrange(participant_id, trial_number, click_index)

# A convenient analysis subset. No rows are removed from trials.csv; instead,
# this object contains only trials with a timely selection and a mouse trace.
trials_analysis <- trials |>
  filter(valid_response)

# -----------------------------------------------------------------------------
# Save outputs and report basic quality checks
# -----------------------------------------------------------------------------

write_csv(participants, file.path(output_dir, "participants.csv"), na = "")
write_csv(
  experiment_events,
  file.path(output_dir, "experiment_events.csv"),
  na = ""
)
write_csv(trials, file.path(output_dir, "trials.csv"), na = "")
write_csv(
  mouse_trajectories,
  file.path(output_dir, "mouse_trajectories.csv"),
  na = ""
)
write_csv(clicks, file.path(output_dir, "clicks.csv"), na = "")

message("Participants: ", nrow(participants))
message("Experimental trials: ", nrow(trials))
message("Valid responses: ", sum(trials$valid_response, na.rm = TRUE))
message("Decoded mouse samples: ", nrow(mouse_trajectories))
message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))

