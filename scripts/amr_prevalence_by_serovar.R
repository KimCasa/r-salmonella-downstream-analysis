library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)

#year parsing if present, otherwise falls back to collection.date
#coalesces possible identifier fields
df <- merged %>%
  mutate(
    Year = suppressWarnings(as.integer(ifelse(!is.na(Year), Year, Collection.date))),
    isolate_id = coalesce(as.character(Isolate), as.character(Isolate.identifiers),
                          as.character(BioSample), as.character(Strain), as.character(Run)),
    Serovar = str_squish(str_to_title(Serovar)),
    Serovar = ifelse(str_detect(Serovar, regex("^typh", TRUE)), "Typhimurium", Serovar)
    #
  )

#Pick a usable AMR-gene column automatically
candidates <- c("Element.symbol","Element.name","AMR_Gene","AMR_gene_plot")
candidates <- candidates[candidates %in% names(df)]

score_col <- function(v) {
  v <- as.character(v)
  v <- str_squish(v)
  non_empty <- sum(!is.na(v) & v != "")
  plusminus <- sum(v %in% c("+","-"), na.rm=TRUE)
  uniq <- length(unique(na.omit(v[v!="" & !v %in% c("+","-")])))
  data.frame(non_empty=non_empty, plusminus=plusminus, uniq=uniq)
}

diag <- lapply(candidates, function(nm) {
  s <- score_col(df[[nm]]); s$col <- nm; s
}) %>% bind_rows()

# Prefer a column with many non-empty, many uniques, and NOT mostly +/-.
pick <- diag %>%
  mutate(score = non_empty + 10*uniq - 100*plusminus) %>%
  arrange(desc(score)) %>% slice(1)

gene_col <- if (nrow(pick)) pick$col[1] else NA_character_
cat("Chosen gene column:", gene_col, "\n")
if (is.na(gene_col)) gene_col <- "Element.symbol"  # default guess

df$gene <- df[[gene_col]] %>% as.character() %>% str_squish()
# If the chosen column is actually junk (mostly +/-), set to NA to trigger fallback:
if (mean(df$gene %in% c("+","-") | is.na(df$gene) | df$gene=="") > 0.9) {
  df$gene <- NA_character_
}


#keep likely AMR alleles; drop virulence/efflux/housekeeping that swamp the ranking
amr_regex  <- regex(
  "^(bla|ampC|kpc|ctx|tem|shv|oxa|ndm|vim|imp|ges|tet|sul|aac|aad|aph|ant|strA|strB|cat|cml|erm|mph|mef|msr|qnr|oqx|fosA|dfr|folP|floR|mcr|van|lnu|optrA|cfr|rmt)",
  ignore_case = TRUE
)
drop_regex <- regex(
  "^(mdsA|mdsB|acr.?|tolC|invA|sinH|sspH\\d*|sseK\\d*|lpf[B-Z]|avrA|iro[BC]|gol[ST])$",
  ignore_case = TRUE
)

p_newport      <- make_serovar_plot(plot_base, "Newport",      top_n = 8, end_year = 2024)
p_typhimurium  <- make_serovar_plot(plot_base, "Typhimurium",  top_n = 8, end_year = 2024)

plot_base <- plot_base %>%
  filter(str_detect(gene, amr_regex)) %>%      # keep AMR alleles
  filter(!str_detect(gene, drop_regex))        # remove near-ubiquitous non-AMR hits

#Pick serovars to show (prefer Dublin & Typhimurium)
targets <- c("Cerro","Typhimurium", "Newport")
present_targets <- intersect(targets, unique(plot_base$Serovar))
if (!length(present_targets)) {
  present_targets <- names(sort(table(plot_base$Serovar), decreasing=TRUE))[1:min(2, length(unique(plot_base$Serovar)))]
  message("Using most abundant serovars: ", paste(present_targets, collapse=", "))
}
plot_base <- plot_base %>% filter(Serovar %in% present_targets)

stopifnot(nrow(plot_base) > 0)

#Prevalence calculations
denom <- plot_base %>% group_by(Serovar, Year) %>%
  summarise(total_isolates = n_distinct(isolate_id), .groups="drop")

num <- plot_base %>% group_by(Serovar, Year, gene) %>%
  summarise(isolates_with_gene = n_distinct(isolate_id), .groups="drop")

prev <- left_join(num, denom, by=c("Serovar","Year")) %>%
  mutate(prevalence = isolates_with_gene/total_isolates)

#Top-N genes per serovar 
TOP_N <- 8
top_genes <- prev %>% group_by(Serovar, gene) %>%
  summarise(overall_isolates = sum(isolates_with_gene), .groups="drop_last") %>%
  filter(overall_isolates>0) %>%
  group_by(Serovar) %>% slice_max(overall_isolates, n=TOP_N, with_ties=FALSE) %>% ungroup()

prev_top <- inner_join(prev, top_genes %>% select(Serovar, gene), by=c("Serovar","gene"))

# complete years so lines are continuous
yr <- sort(unique(prev_top$Year))
prev_complete <- prev_top %>% group_by(Serovar, gene) %>%
  complete(Year = yr, fill=list(isolates_with_gene=0, total_isolates=NA, prevalence=0)) %>%
  ungroup()

p <- ggplot(prev_complete, aes(Year, prevalence, color = gene, group = gene)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  facet_wrap(~ Serovar, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Yearly Prevalence within Each Serovar",
    x = "Year",
    y = "Prevalence of isolates with allele",
    color = "AMR allele"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p)