args <- commandArgs(trailingOnly = TRUE)

# Portable project root: this runner lives in 00_project_setup/.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(gt)
})

run_date_tag <- format(Sys.Date(), "%Y%m%d")

default_day03_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_patient_day0_day3_integrated_longitudinal_itt460_v.xlsx"
default_day60_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_patient_day0_day60_integrated_longitudinal_itt460_.xlsx"
default_severity_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_primary_severity_L0_L3_itt460.xlsx"
default_itt_rds_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/itt.dataset.RDS"
default_out_dir <- file.path("/Users/bing/Documents/中介效应/outputs", paste0("regardvap_manuscript_pipeline_", run_date_tag))

day03_path <- if (length(args) >= 1) args[1] else default_day03_path
day60_path <- if (length(args) >= 2) args[2] else default_day60_path
severity_path <- if (length(args) >= 3) args[3] else default_severity_path
itt_rds_path <- if (length(args) >= 4) args[4] else default_itt_rds_path
out_dir <- if (length(args) >= 5) args[5] else default_out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

section_dirs <- list(
  sec31 = file.path(out_dir, "3.1_study_population_and_baseline"),
  sec32 = file.path(out_dir, "3.2_longitudinal_trajectories"),
  sec33 = file.path(out_dir, "3.3_longitudinal_mediation_analysis"),
  sec34 = file.path(out_dir, "3.4_heterogeneity_of_treatment_effect"),
  scratch = file.path(out_dir, "_scratch")
)
invisible(lapply(section_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

analysis_core_dir <- file.path(section_dirs[["scratch"]], "longitudinal_core")
hte_dir <- file.path(section_dirs[["scratch"]], "hte_interaction")

message_block <- function(...) {
  cat(paste0(..., collapse = ""), "\n")
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

copy_required <- function(src, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) stop("Failed to copy file: ", src, " -> ", dest, call. = FALSE)
  invisible(dest)
}

run_r_script <- function(script_path, script_args) {
  assert_file_exists(script_path, "Script")
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  cmd_args <- shQuote(c(script_path, script_args))
  out <- system2(rscript_bin, args = cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (length(out) > 0) {
    message_block(paste(out, collapse = "\n"))
  }
  if (status != 0L) {
    stop("Child script failed: ", script_path, call. = FALSE)
  }
  invisible(out)
}

safe_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[is.nan(out)] <- NA_real_
  out
}

make_figure4 <- function(table3_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  df <- read.csv(table3_path, stringsAsFactors = FALSE, check.names = FALSE)
  risk_df <- data.frame(
    panel = "Panel A. Counterfactual risks",
    label = c("R(1, G1)", "R(1, G0)", "R(0, G0)"),
    estimate = c(df$R_1_G1[1], df$R_1_G0[1], df$R_0_G0[1]),
    stringsAsFactors = FALSE
  )
  effect_df <- data.frame(
    panel = "Panel B. TE decomposition",
    label = c("TE", "IDE", "IIE"),
    estimate = c(df$TE[1], df$IDE[1], df$IIE[1]),
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(risk_df, effect_df)
  plot_df$label <- factor(plot_df$label, levels = c("R(1, G1)", "R(1, G0)", "R(0, G0)", "TE", "IDE", "IIE"))

  fig <- ggplot(plot_df, aes(x = label, y = estimate, fill = label)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "gray45", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", estimate)), vjust = ifelse(plot_df$estimate >= 0, -0.4, 1.2), size = 4) +
    facet_wrap(~panel, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c(
      "R(1, G1)" = "#1F4E79",
      "R(1, G0)" = "#5B8DB8",
      "R(0, G0)" = "#9AA5B1",
      "TE" = "#203A57",
      "IDE" = "#0E7C86",
      "IIE" = "#D98E04"
    )) +
    labs(
      title = "Figure 4. Counterfactual risks and interventional effect decomposition",
      x = NULL,
      y = "Risk or risk difference"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )

  png_path <- file.path(out_dir, "figure4_te_decomposition.png")
  csv_path <- file.path(out_dir, "figure4_te_decomposition_data.csv")
  ggsave(png_path, fig, width = 11, height = 5.5, dpi = 300)
  write.csv(plot_df, csv_path, row.names = FALSE)
  list(png = png_path, csv = csv_path)
}

make_figure1_svg <- function(out_path) {
  lines <- c(
    '<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1120" viewBox="0 0 1600 1120">',
    '  <defs>',
    '    <style>',
    '      .title { font-family: Georgia, "Times New Roman", serif; font-size: 42px; font-weight: 700; fill: #111; }',
    '      .label { font-family: Georgia, "Times New Roman", serif; font-size: 28px; fill: #222; }',
    '      .small { font-family: Georgia, "Times New Roman", serif; font-size: 22px; fill: #222; }',
    '      .tablehead { font-family: Georgia, "Times New Roman", serif; font-size: 28px; font-weight: 700; fill: #222; }',
    '      .cellhead { font-family: Georgia, "Times New Roman", serif; font-size: 23px; font-weight: 400; fill: #222; }',
    '      .cell { font-family: Georgia, "Times New Roman", serif; font-size: 22px; fill: #222; }',
    '      .foot { font-family: Georgia, "Times New Roman", serif; font-size: 18px; fill: #222; }',
    '      .line { stroke: #222; stroke-width: 2.3; fill: none; }',
    '      .thin { stroke: #333; stroke-width: 1.6; fill: none; }',
    '      .arrow { marker-end: url(#arrow); }',
    '    </style>',
    '    <marker id="arrow" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">',
    '      <path d="M0,0 L10,4 L0,8 z" fill="#222"/>',
    '    </marker>',
    '  </defs>',
    '',
    '  <rect x="0" y="0" width="1600" height="1120" fill="#fff"/>',
    '  <text x="110" y="95" class="title">Timeline of measurements</text>',
    '  <text x="165" y="165" class="label" font-style="italic">C, A, L0</text>',
    '  <text x="455" y="165" class="label" font-style="italic">L0, M0</text>',
    '  <text x="690" y="165" class="label" font-style="italic">L1, M1</text>',
    '  <text x="980" y="165" class="label" font-style="italic">L2, M2</text>',
    '  <text x="1240" y="165" class="label" font-style="italic">L3, M3</text>',
    '  <text x="1460" y="165" class="label" font-style="italic">Y</text>',
    '  <line x1="205" y1="190" x2="205" y2="465" class="line"/>',
    '  <line x1="470" y1="190" x2="470" y2="465" class="line"/>',
    '  <line x1="715" y1="190" x2="715" y2="465" class="line"/>',
    '  <line x1="1005" y1="190" x2="1005" y2="465" class="line"/>',
    '  <line x1="1270" y1="190" x2="1270" y2="465" class="line"/>',
    '  <line x1="1500" y1="190" x2="1500" y2="465" class="line"/>',
    '  <line x1="470" y1="300" x2="1270" y2="300" class="line arrow"/>',
    '  <text x="775" y="250" text-anchor="middle" class="label">Daily clinical state and appropriate-treatment process</text>',
    '  <text x="775" y="286" text-anchor="middle" class="label">measured during the Day 0-3 mediation window</text>',
    '  <line x1="205" y1="390" x2="1500" y2="390" class="line arrow"/>',
    '  <text x="840" y="352" text-anchor="middle" class="label">Day 0-3 longitudinal mediation window</text>',
    '  <text x="840" y="440" text-anchor="middle" class="label">Follow-up for 60-day mortality</text>',
    '  <line x1="165" y1="450" x2="1530" y2="450" class="line arrow"/>',
    '  <text x="160" y="508" class="label">Baseline</text>',
    '  <text x="435" y="508" class="label">Day 0</text>',
    '  <text x="684" y="508" class="label">Day 1</text>',
    '  <text x="973" y="508" class="label">Day 2</text>',
    '  <text x="1238" y="508" class="label">Day 3</text>',
    '  <text x="1454" y="508" class="label">Day 60</text>',
    '  <text x="105" y="610" class="tablehead">Variables measured at each time point</text>',
    '  <rect x="105" y="630" width="1430" height="392" class="line" fill="none"/>',
    '  <line x1="105" y1="710" x2="1535" y2="710" class="thin"/>',
    '  <line x1="105" y1="762" x2="1535" y2="762" class="thin"/>',
    '  <line x1="105" y1="814" x2="1535" y2="814" class="thin"/>',
    '  <line x1="105" y1="866" x2="1535" y2="866" class="thin"/>',
    '  <line x1="105" y1="918" x2="1535" y2="918" class="thin"/>',
    '  <line x1="105" y1="970" x2="1535" y2="970" class="thin"/>',
    '  <line x1="240" y1="630" x2="240" y2="1022" class="thin"/>',
    '  <line x1="500" y1="630" x2="500" y2="1022" class="thin"/>',
    '  <line x1="770" y1="630" x2="770" y2="1022" class="thin"/>',
    '  <line x1="1360" y1="630" x2="1360" y2="1022" class="thin"/>',
    '  <text x="172" y="684" text-anchor="middle" class="cellhead">Time point</text>',
    '  <text x="370" y="684" text-anchor="middle" class="cellhead">Exposure</text>',
    '  <text x="635" y="684" text-anchor="middle" class="cellhead">Mediator</text>',
    '  <text x="1065" y="684" text-anchor="middle" class="cellhead">Time-varying clinical state</text>',
    '  <text x="1448" y="684" text-anchor="middle" class="cellhead">Outcome</text>',
    '  <text x="172" y="744" text-anchor="middle" class="cell">Baseline</text>',
    '  <text x="370" y="744" text-anchor="middle" class="cell" font-style="italic">A</text>',
    '  <text x="635" y="744" text-anchor="middle" class="cell">-</text>',
    '  <text x="1065" y="744" text-anchor="middle" class="cell" font-style="italic">L0</text>',
    '  <text x="1448" y="744" text-anchor="middle" class="cell">-</text>',
    '  <text x="172" y="796" text-anchor="middle" class="cell">Day 0</text>',
    '  <text x="370" y="796" text-anchor="middle" class="cell" font-style="italic">A</text>',
    '  <text x="635" y="796" text-anchor="middle" class="cell" font-style="italic">M0</text>',
    '  <text x="1065" y="796" text-anchor="middle" class="cell" font-style="italic">L0</text>',
    '  <text x="1448" y="796" text-anchor="middle" class="cell">-</text>',
    '  <text x="172" y="848" text-anchor="middle" class="cell">Day 1</text>',
    '  <text x="370" y="848" text-anchor="middle" class="cell" font-style="italic">A</text>',
    '  <text x="635" y="848" text-anchor="middle" class="cell" font-style="italic">M1</text>',
    '  <text x="1065" y="848" text-anchor="middle" class="cell" font-style="italic">L1</text>',
    '  <text x="1448" y="848" text-anchor="middle" class="cell">-</text>',
    '  <text x="172" y="900" text-anchor="middle" class="cell">Day 2</text>',
    '  <text x="370" y="900" text-anchor="middle" class="cell" font-style="italic">A</text>',
    '  <text x="635" y="900" text-anchor="middle" class="cell" font-style="italic">M2</text>',
    '  <text x="1065" y="900" text-anchor="middle" class="cell" font-style="italic">L2</text>',
    '  <text x="1448" y="900" text-anchor="middle" class="cell">-</text>',
    '  <text x="172" y="952" text-anchor="middle" class="cell">Day 3</text>',
    '  <text x="370" y="952" text-anchor="middle" class="cell" font-style="italic">A</text>',
    '  <text x="635" y="952" text-anchor="middle" class="cell" font-style="italic">M3</text>',
    '  <text x="1065" y="952" text-anchor="middle" class="cell" font-style="italic">L3</text>',
    '  <text x="1448" y="952" text-anchor="middle" class="cell">-</text>',
    '  <text x="172" y="1004" text-anchor="middle" class="cell">Day 60</text>',
    '  <text x="370" y="1004" text-anchor="middle" class="cell" font-style="italic">A</text>',
    '  <text x="635" y="1004" text-anchor="middle" class="cell">-</text>',
    '  <text x="1065" y="1004" text-anchor="middle" class="cell">-</text>',
    '  <text x="1448" y="1004" text-anchor="middle" class="cell" font-style="italic">Y</text>',
    '  <text x="105" y="1052" class="foot">V = baseline covariates; A = baseline carbapenem resistance; L0-L3 = daily severity states; M0-M3 = daily appropriate treatment; Y = 60-day mortality.</text>',
    '  <text x="105" y="1077" class="foot">The Day 0-3 mediation window captures the early treatment process after baseline microbiologic resistance status is defined.</text>',
    '  <text x="105" y="1100" class="foot">This figure is generated directly from code in the manuscript pipeline and is not imported from a static asset.</text>',
    '</svg>'
  )
  writeLines(lines, out_path)
}

make_figure2_svg <- function(out_path) {
  lines <- c(
    '<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="980" viewBox="0 0 1800 980">',
    '  <defs>',
    '    <style>',
    '      .paper { font-family: "Times New Roman", Times, serif; fill: #1f1f1f; }',
    '      .node { font-size: 38px; font-style: italic; }',
    '      .smallnode { font-size: 34px; font-style: italic; }',
    '      .caplead { font: 700 29px Arial, sans-serif; fill: #111; }',
    '      .captext { font: 400 27px Arial, sans-serif; fill: #111; }',
    '      .edge { stroke: #202020; stroke-width: 3.1; fill: none; marker-end: url(#arrow); }',
    '    </style>',
    '    <marker id="arrow" markerWidth="11" markerHeight="11" refX="8.5" refY="3.8" orient="auto" markerUnits="strokeWidth">',
    '      <path d="M0,0 L9,3.8 L0,7.6 z" fill="#202020"/>',
    '    </marker>',
    '  </defs>',
    '  <rect width="1800" height="980" fill="#ffffff"/>',
    '  <g class="paper">',
    '    <text x="300" y="430" class="node">A</text>',
    '    <text x="520" y="300" class="smallnode">M0</text>',
    '    <text x="760" y="300" class="smallnode">M1</text>',
    '    <text x="1000" y="300" class="smallnode">M2</text>',
    '    <text x="1240" y="300" class="smallnode">M3</text>',
    '    <text x="520" y="500" class="smallnode">L0</text>',
    '    <text x="760" y="500" class="smallnode">L1</text>',
    '    <text x="1000" y="500" class="smallnode">L2</text>',
    '    <text x="1240" y="500" class="smallnode">L3</text>',
    '    <text x="1485" y="400" class="smallnode">Y</text>',
    '  </g>',
    '  <line x1="342" y1="415" x2="500" y2="290" class="edge"/>',
    '  <line x1="342" y1="420" x2="500" y2="490" class="edge"/>',
    '  <path d="M342 404 C500 325, 650 305, 740 292" class="edge"/>',
    '  <path d="M342 412 C530 390, 675 465, 740 492" class="edge"/>',
    '  <path d="M342 396 C635 250, 945 250, 980 292" class="edge"/>',
    '  <path d="M342 422 C610 455, 920 495, 980 492" class="edge"/>',
    '  <path d="M342 388 C740 205, 1215 210, 1220 292" class="edge"/>',
    '  <path d="M342 430 C760 510, 1180 518, 1220 492" class="edge"/>',
    '  <path d="M342 384 C825 160, 1385 190, 1470 384" class="edge"/>',
    '  <line x1="555" y1="300" x2="735" y2="300" class="edge"/>',
    '  <line x1="795" y1="300" x2="975" y2="300" class="edge"/>',
    '  <line x1="1035" y1="300" x2="1215" y2="300" class="edge"/>',
    '  <line x1="555" y1="500" x2="735" y2="500" class="edge"/>',
    '  <line x1="795" y1="500" x2="975" y2="500" class="edge"/>',
    '  <line x1="1035" y1="500" x2="1215" y2="500" class="edge"/>',
    '  <line x1="540" y1="320" x2="540" y2="475" class="edge"/>',
    '  <line x1="780" y1="320" x2="780" y2="475" class="edge"/>',
    '  <line x1="1020" y1="320" x2="1020" y2="475" class="edge"/>',
    '  <line x1="1260" y1="320" x2="1260" y2="475" class="edge"/>',
    '  <path d="M550 500 C625 465, 690 355, 750 320" class="edge"/>',
    '  <path d="M790 500 C865 465, 930 355, 990 320" class="edge"/>',
    '  <path d="M1030 500 C1105 465, 1170 355, 1230 320" class="edge"/>',
    '  <path d="M555 300 C635 360, 700 455, 750 485" class="edge"/>',
    '  <path d="M795 300 C875 360, 940 455, 990 485" class="edge"/>',
    '  <path d="M1035 300 C1115 360, 1180 455, 1230 485" class="edge"/>',
    '  <path d="M550 500 C735 500, 905 500, 1230 500" class="edge"/>',
    '  <path d="M550 500 C720 470, 900 435, 1470 390" class="edge"/>',
    '  <path d="M790 500 C930 470, 1080 435, 1470 392" class="edge"/>',
    '  <path d="M1030 500 C1160 470, 1280 435, 1470 394" class="edge"/>',
    '  <path d="M1270 500 C1360 470, 1430 430, 1470 400" class="edge"/>',
    '  <path d="M555 300 C760 205, 1115 205, 1470 385" class="edge"/>',
    '  <path d="M795 300 C970 235, 1200 235, 1470 388" class="edge"/>',
    '  <path d="M1035 300 C1175 255, 1315 270, 1470 392" class="edge"/>',
    '  <line x1="1275" y1="300" x2="1470" y2="396" class="edge"/>',
    '  <text x="120" y="835" class="caplead">Figure 2.</text>',
    '  <text x="280" y="835" class="captext">Causal diagram for the longitudinal mediation analysis in the REGARD-VAP study.</text>',
    '  <text x="120" y="880" class="captext">V is omitted from the diagram but is assumed to affect A, L0-L3, M0-M3, and Y.</text>',
    '  <text x="120" y="925" class="captext">A denotes baseline carbapenem resistance; L0-L3 denote daily severity states; M0-M3 denote daily appropriate treatment; Y denotes 60-day mortality.</text>',
    '  <text x="120" y="970" class="captext">This figure is generated directly from code in the manuscript pipeline and is not imported from a static asset.</text>',
    '</svg>'
  )
  writeLines(lines, out_path)
}

write_section_readme <- function(path, lines) {
  writeLines(lines, path)
}

gt_save_local <- function(tbl, html_path) {
  old_cache <- Sys.getenv("R_SASS_CACHE_DIR", unset = "")
  on.exit(Sys.setenv(R_SASS_CACHE_DIR = old_cache), add = TRUE)
  Sys.setenv(R_SASS_CACHE_DIR = file.path(tempdir(), "r-sass-cache"))
  gt::gtsave(tbl, html_path)
}

make_pretty_table2 <- function(input_csv, out_dir) {
  df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
  df$group <- ifelse(df$group == "Non-resistant", "Non-resistant", "Resistant")
  df$day_label <- paste0("Day ", df$day)

  display_df <- data.frame(
    Group = df$group,
    Day = df$day_label,
    `Patient-days` = df$n_rows,
    `Treatment observed` = df$n_treatment_observed,
    `Severity observed` = df$n_severity_observed,
    `Micro info observed` = df$n_micro_info_observed,
    `Appropriate treatment, n (%)` = df$appropriate_treatment,
    `Daily severity (primary), mean (SD)` = df$daily_severity_main,
    `Daily severity (alternative), mean (SD)` = df$daily_severity_alt,
    `Microbiology info available, n (%)` = df$microbiology_info_available,
    `Prior carbapenem-R info, n (%)` = df$prior_carbapenem_r_info,
    `Cumulative appropriate days, median [IQR]` = df$cumulative_appropriate_days,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tbl <- display_df |>
    gt::gt(rowname_col = "Day", groupname_col = "Group") |>
    gt::tab_header(
      title = gt::md("**Table 2. Day 0-3 treatment and severity summary**"),
      subtitle = "Observed longitudinal processes stratified by baseline carbapenem resistance."
    ) |>
    gt::tab_spanner(
      label = "Observed data support",
      columns = c(`Patient-days`, `Treatment observed`, `Severity observed`, `Micro info observed`)
    ) |>
    gt::tab_spanner(
      label = "Treatment and clinical state summary",
      columns = c(
        `Appropriate treatment, n (%)`,
        `Daily severity (primary), mean (SD)`,
        `Daily severity (alternative), mean (SD)`,
        `Microbiology info available, n (%)`,
        `Prior carbapenem-R info, n (%)`,
        `Cumulative appropriate days, median [IQR]`
      )
    ) |>
    gt::cols_align(align = "center", columns = everything()) |>
    gt::tab_source_note(
      source_note = gt::md(
        "Primary severity uses the main `L_t` definition; alternative severity uses the sensitivity-analysis `L_t` definition."
      )
    ) |>
    gt::tab_source_note(
      source_note = gt::md(
        "`Patient-days` is the total number of rows in that group/day stratum; the three `observed` columns show the non-missing denominators used for the corresponding summaries."
      )
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#1E3A5F"),
        gt::cell_text(color = "white", weight = "bold")
      ),
      locations = gt::cells_column_labels(everything())
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#D9E7F5"),
        gt::cell_text(weight = "bold", color = "#17324D")
      ),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = "#F7FAFC"),
      locations = gt::cells_body(rows = seq(2, nrow(display_df), by = 2))
    ) |>
    gt::tab_options(
      table.font.names = c("Arial", "Helvetica", "sans-serif"),
      table.font.size = 13,
      heading.title.font.size = 18,
      heading.subtitle.font.size = 12,
      data_row.padding = gt::px(7),
      row_group.font.size = 13,
      source_notes.font.size = 11,
      column_labels.background.color = "#1E3A5F",
      table.border.top.color = "#1E3A5F",
      table.border.bottom.color = "#1E3A5F",
      row_group.background.color = "#D9E7F5",
      table.width = gt::pct(100)
    )

  write.csv(display_df, file.path(out_dir, "table2_day0_day3_treatment_severity_summary_pretty.csv"), row.names = FALSE)
  gt_save_local(tbl, file.path(out_dir, "table2_day0_day3_treatment_severity_summary_pretty.html"))
}

make_pretty_table3 <- function(input_csv, out_dir) {
  df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

  display_df <- data.frame(
    `Complete cases` = df$n_complete,
    `Risk under R(1, G1)` = sprintf("%.3f", df$R_1_G1),
    `Risk under R(1, G0)` = sprintf("%.3f", df$R_1_G0),
    `Risk under R(0, G0)` = sprintf("%.3f", df$R_0_G0),
    `Total effect (TE)` = sprintf("%.3f", df$TE),
    `Indirect effect (IIE)` = sprintf("%.3f", df$IIE),
    `Direct effect (IDE)` = sprintf("%.3f", df$IDE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tbl <- display_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Table 3. G-formula estimates of counterfactual risks and effect decomposition**"),
      subtitle = "Primary longitudinal mediation analysis on the 60-day mortality risk-difference scale."
    ) |>
    gt::tab_spanner(
      label = "Counterfactual risks",
      columns = c(`Risk under R(1, G1)`, `Risk under R(1, G0)`, `Risk under R(0, G0)`)
    ) |>
    gt::tab_spanner(
      label = "Effect decomposition",
      columns = c(`Total effect (TE)`, `Indirect effect (IIE)`, `Direct effect (IDE)`)
    ) |>
    gt::cols_align(align = "center", columns = everything()) |>
    gt::tab_source_note(
      source_note = gt::md(
        "`R(1, G1)` = risk under resistant exposure with the resistant-group treatment trajectory law; `R(1, G0)` = risk under resistant exposure with the reference treatment law; `R(0, G0)` = risk under reference exposure with the reference treatment law."
      )
    ) |>
    gt::tab_source_note(
      source_note = gt::md(
        "Effect decomposition is defined on the risk-difference scale: `TE = IIE + IDE`."
      )
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#1E3A5F"),
        gt::cell_text(color = "white", weight = "bold")
      ),
      locations = gt::cells_column_labels(everything())
    ) |>
    gt::tab_options(
      table.font.names = c("Arial", "Helvetica", "sans-serif"),
      table.font.size = 13,
      heading.title.font.size = 18,
      heading.subtitle.font.size = 12,
      data_row.padding = gt::px(8),
      source_notes.font.size = 11,
      column_labels.background.color = "#1E3A5F",
      table.border.top.color = "#1E3A5F",
      table.border.bottom.color = "#1E3A5F",
      table.width = gt::pct(100)
    )

  write.csv(display_df, file.path(out_dir, "table3_counterfactual_risks_te_ide_iie_pretty.csv"), row.names = FALSE)
  gt_save_local(tbl, file.path(out_dir, "table3_counterfactual_risks_te_ide_iie_pretty.html"))
}

make_pretty_table4 <- function(input_csv, out_dir) {
  df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

  df$Specification <- ifelse(
    df$include_ot,
    "Main L_t + O_t",
    ifelse(df$lt_variant == "alt", "Alternative L_t", ifelse(df$mediator_history == "none", "No lagged mediator history", "Main specification"))
  )

  display_df <- data.frame(
    Specification = df$Specification,
    `Complete cases` = df$n_complete,
    `Active covariates` = df$active_covars,
    `Dropped covariates` = ifelse(df$dropped_covars == "", "-", df$dropped_covars),
    `Risk under R(1, G1)` = sprintf("%.3f", df$R_1_G1),
    `Risk under R(1, G0)` = sprintf("%.3f", df$R_1_G0),
    `Risk under R(0, G0)` = sprintf("%.3f", df$R_0_G0),
    `Total effect (TE)` = sprintf("%.3f", df$TE),
    `Indirect effect (IIE)` = sprintf("%.3f", df$IIE),
    `Direct effect (IDE)` = sprintf("%.3f", df$IDE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tbl <- display_df |>
    gt::gt(rowname_col = "Specification") |>
    gt::tab_header(
      title = gt::md("**Table 4. Sensitivity analyses for the longitudinal g-formula**"),
      subtitle = "Alternative severity, mediator-history, and O_t specifications."
    ) |>
    gt::tab_spanner(
      label = "Model diagnostics",
      columns = c(`Complete cases`, `Active covariates`, `Dropped covariates`)
    ) |>
    gt::tab_spanner(
      label = "Counterfactual risks",
      columns = c(`Risk under R(1, G1)`, `Risk under R(1, G0)`, `Risk under R(0, G0)`)
    ) |>
    gt::tab_spanner(
      label = "Effect decomposition",
      columns = c(`Total effect (TE)`, `Indirect effect (IIE)`, `Direct effect (IDE)`)
    ) |>
    gt::cols_align(align = "center", columns = everything()) |>
    gt::tab_source_note(
      source_note = gt::md(
        "`Main L_t + O_t` includes the microbiology-information term; `Alternative L_t` uses the alternative severity definition; `No lagged mediator history` removes the lagged treatment history from the mediator model."
      )
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#1E3A5F"),
        gt::cell_text(color = "white", weight = "bold")
      ),
      locations = gt::cells_column_labels(everything())
    ) |>
    gt::tab_options(
      table.font.names = c("Arial", "Helvetica", "sans-serif"),
      table.font.size = 12,
      heading.title.font.size = 18,
      heading.subtitle.font.size = 12,
      data_row.padding = gt::px(8),
      source_notes.font.size = 11,
      column_labels.background.color = "#1E3A5F",
      table.border.top.color = "#1E3A5F",
      table.border.bottom.color = "#1E3A5F",
      table.width = gt::pct(100)
    )

  write.csv(display_df, file.path(out_dir, "table4_sensitivity_analyses_pretty.csv"), row.names = FALSE)
  gt_save_local(tbl, file.path(out_dir, "table4_sensitivity_analyses_pretty.html"))
}

assert_file_exists(day03_path, "Day 0-3 file")
assert_file_exists(day60_path, "Day 0-60 file")
assert_file_exists(severity_path, "Severity file")
assert_file_exists(itt_rds_path, "ITT RDS file")

message_block("Running longitudinal mediation core pipeline...")
run_r_script(
  file.path(project_root, "04_longitudinal_mediation", "01_run_regardvap_longitudinal_mediation_v3_no_wgs.R"),
  c(day03_path, day60_path, severity_path, itt_rds_path, analysis_core_dir)
)

analysis_wide_path <- file.path(analysis_core_dir, "step00_data_prep", "analysis_panel_wide_for_gformula.csv")
assert_file_exists(analysis_wide_path, "Longitudinal analysis panel")

message_block("Running interaction-model HTE pipeline...")
run_r_script(
  file.path(project_root, "06_heterogeneity_analyses", "01_run_hte_interaction_models.R"),
  c(day03_path, day60_path, hte_dir)
)

message_block("Collecting manuscript-structured outputs...")

# 3.1 需要找danny看一下内容
#write_section_readme(
#  file.path(section_dirs[["sec31"]], "README_3.1.md"),
# c(
#   "# 3.1 Study population and baseline characteristics",
#    "",
#    "- Baseline section intentionally left blank in this version.",
#    "- You asked to fill the baseline part later after we review Danny's code.",
#    "- No baseline table is exported by the current manuscript pipeline run."
#  )
#)

# 3.2
make_figure1_svg(file.path(section_dirs[["sec32"]], "figure1_timeline.svg"))
make_figure2_svg(file.path(section_dirs[["sec32"]], "figure2_dag.svg"))
copy_required(
  file.path(analysis_core_dir, "step02_longitudinal_summary", "table2_longitudinal_day0_day3_summary.csv"),
  file.path(section_dirs[["sec32"]], "table2_day0_day3_treatment_severity_summary.csv")
)
copy_required(
  file.path(analysis_core_dir, "step02_longitudinal_summary", "table2_observed_trajectory_distribution.csv"),
  file.path(section_dirs[["sec32"]], "table2_observed_trajectory_distribution.csv")
)
copy_required(
  file.path(analysis_core_dir, "step02_longitudinal_summary", "figure3_observed_processes.png"),
  file.path(section_dirs[["sec32"]], "figure3_observed_trajectories.png")
)
figure3_svg_src <- file.path(analysis_core_dir, "step02_longitudinal_summary", "figure3_observed_processes.svg")
if (file.exists(figure3_svg_src)) {
  copy_required(figure3_svg_src, file.path(section_dirs[["sec32"]], "figure3_observed_trajectories.svg"))
}
copy_required(
  file.path(analysis_core_dir, "step02_longitudinal_summary", "figure3_daywise_summary.csv"),
  file.path(section_dirs[["sec32"]], "figure3_daywise_summary.csv")
)
make_pretty_table2(
  file.path(section_dirs[["sec32"]], "table2_day0_day3_treatment_severity_summary.csv"),
  section_dirs[["sec32"]]
)
write_section_readme(
  file.path(section_dirs[["sec32"]], "README_3.2.md"),
  c(
    "# 3.2 Longitudinal trajectories",
    "",
    "- Figure 1: figure1_timeline.svg",
    "- Figure 2: figure2_dag.svg",
    "- Table 2: table2_day0_day3_treatment_severity_summary.csv",
    "- Pretty Table 2 HTML: table2_day0_day3_treatment_severity_summary_pretty.html",
    "- Figure 3: figure3_observed_trajectories.png",
    "- Supporting files: figure3_daywise_summary.csv, table2_observed_trajectory_distribution.csv"
  )
)

# 3.3
table3_source <- file.path(analysis_core_dir, "step03_main_gformula", "table3_main_gformula_estimates.csv")
copy_required(table3_source, file.path(section_dirs[["sec33"]], "table3_counterfactual_risks_te_ide_iie.csv"))
copy_required(
  file.path(analysis_core_dir, "step04_sensitivity", "table4_sensitivity_analyses.csv"),
  file.path(section_dirs[["sec33"]], "table4_sensitivity_analyses.csv")
)
figure4_files <- make_figure4(table3_source, section_dirs[["sec33"]])
make_pretty_table3(
  file.path(section_dirs[["sec33"]], "table3_counterfactual_risks_te_ide_iie.csv"),
  section_dirs[["sec33"]]
)
make_pretty_table4(
  file.path(section_dirs[["sec33"]], "table4_sensitivity_analyses.csv"),
  section_dirs[["sec33"]]
)
write_section_readme(
  file.path(section_dirs[["sec33"]], "README_3.3.md"),
  c(
    "# 3.3 Longitudinal mediation analysis",
    "",
    "- Table 3: table3_counterfactual_risks_te_ide_iie.csv",
    "- Pretty Table 3 HTML: table3_counterfactual_risks_te_ide_iie_pretty.html",
    "- Figure 4: figure4_te_decomposition.png",
    "- Table 4: table4_sensitivity_analyses.csv",
    "- Pretty Table 4 HTML: table4_sensitivity_analyses_pretty.html",
    "- Supporting plotting data: figure4_te_decomposition_data.csv"
  )
)

# 3.4
copy_required(
  file.path(analysis_core_dir, "step05_hte", "table5_hte_gformula_estimates.csv"),
  file.path(section_dirs[["sec34"]], "table5_hte_gformula_estimates.csv")
)
hte_copy_targets <- c(
  "table_hte_interaction_global_tests.csv",
  "table_hte_interaction_sofa_curve.csv",
  "table_hte_interaction_bacteria.csv",
  "table_hte_interaction_country.csv",
  "table_hte_interaction_balance_counts.csv",
  "figure_hte_interaction_baseline_sofa.png",
  "figure_hte_interaction_bacteria.png",
  "figure_hte_interaction_country.png"
)
for (nm in hte_copy_targets) {
  src <- file.path(hte_dir, nm)
  if (file.exists(src)) {
    copy_required(src, file.path(section_dirs[["sec34"]], nm))
  }
}
write_section_readme(
  file.path(section_dirs[["sec34"]], "README_3.4.md"),
  c(
    "# 3.4 Heterogeneity of treatment effect",
    "",
    "- Longitudinal g-formula subgroup table: table5_hte_gformula_estimates.csv",
    "- Interaction-model outputs are copied here as sensitivity / supporting HTE analyses."
  )
)

write_section_readme(
  file.path(out_dir, "README_manuscript_pipeline.md"),
  c(
    "# REGARD-VAP manuscript pipeline v1",
    "",
    "This run organizes outputs into the manuscript-facing structure:",
    "- 3.1 Study population and baseline characteristics",
    "- 3.2 Longitudinal trajectories",
    "- 3.3 Longitudinal mediation analysis",
    "- 3.4 Heterogeneity of treatment effect",
    "",
    "Core analysis directories used internally:",
    paste0("- Longitudinal core: ", analysis_core_dir),
    paste0("- HTE interaction: ", hte_dir),
    "",
    "Input files:",
    paste0("- Day 0-3: ", day03_path),
    paste0("- Day 0-60: ", day60_path),
    paste0("- Severity: ", severity_path),
    paste0("- ITT RDS: ", itt_rds_path),
    "",
    "All manuscript-facing figures in this pipeline are generated by code.",
    "Baseline outputs are intentionally omitted in this version pending the separate baseline review."
  )
)

message_block("Manuscript pipeline finished.")
message_block("Output directory:")
message_block(out_dir)
