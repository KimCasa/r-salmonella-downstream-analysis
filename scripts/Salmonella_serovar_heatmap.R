#written by Kimberly Casares 
#the purpose of this code is to create visualizations of bovine isolates
#created a heatmap 

library(tidyverse)
library(lubridate)
library(scales)  # For better axis formatting

# Convert to a tibble if it's not already (makes working with it easier)
bovine_isolates <- as_tibble(bovine_isolates)

# Standardize serotype names first- before any data analysis
bovine_isolates <- bovine_isolates %>%
  mutate(
    # Create a standardized serotype column
    StandardSerovar = case_when(
      # Extract just the serovar name if it's in the full format
      grepl("serovar", !!sym(servar_col), ignore.case = TRUE) ~ 
        str_extract(!!sym(servar_col), "(?<=serovar\\s)\\w+"),
      # Otherwise use the existing name
      TRUE ~ !!sym(servar_col)
    )
  )

# Update the servar_col variable to use the standardized column
servar_col <- "StandardSerovar"

# Fix the date conversion issues using the detected date column
bovine_isolates <- bovine_isolates %>%
  mutate(
    # Try multiple date formats to ensure correct parsing
    CleanDate = case_when(
      # Try ISO format first (YYYY-MM-DD)
      !is.na(as.Date(!!sym(date_col), format = "%Y-%m-%d")) ~ as.Date(!!sym(date_col), format = "%Y-%m-%d"),
      # Try MM/DD/YYYY format
      !is.na(as.Date(!!sym(date_col), format = "%m/%d/%Y")) ~ as.Date(!!sym(date_col), format = "%m/%d/%Y"),
      # Try DD/MM/YYYY format
      !is.na(as.Date(!!sym(date_col), format = "%d/%m/%Y")) ~ as.Date(!!sym(date_col), format = "%d/%m/%Y"),
      # Add more formats if needed
      TRUE ~ as.Date(NA)
    ),
    Year = year(CleanDate)
  )


#identify top serovars

#Get the top 15 serovars for all visualizations
top_serovars <- bovine_isolates %>%
  filter(!is.na(!!sym(servar_col))) %>%
  count(!!sym(servar_col)) %>%
  arrange(desc(n)) %>%
  slice_head(n = 15) %>%
  pull(!!sym(servar_col))

# For heatmap
#gets the complete set of years in data
all_years <- sort(unique(bovine_isolates$Year))

# Create a complete grid of all serovar-year combinations
complete_grid <- expand.grid(
  Year = all_years,
  Serovar = top_serovars
)
heatmap_data_complete <- bovine_isolates %>%
  filter(!is.na(Year)) %>%
  filter(!!sym(servar_col) %in% top_serovars) %>%
  group_by(Year, !!sym(servar_col)) %>%
  summarise(Count = n(), .groups = "drop") %>%
  rename(Serovar = !!sym(servar_col)) %>%
  right_join(complete_grid, by = c("Year", "Serovar")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count))

# Order serovars by total abundance
serovar_order <- heatmap_data_complete %>%
  group_by(Serovar) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  arrange(desc(TotalCount)) %>%
  pull(Serovar)

heatmap_data_complete$Serovar <- factor(
  heatmap_data_complete$Serovar,
  levels = serovar_order
)

# Heatmap
max_count <- max(heatmap_data_complete$Count)

zero_color <- "#CCCCCC"
low_color <- "#EDF8FB"
high_color <- "#08306B"

heatmap_plot <- ggplot(
  heatmap_data_complete,
  aes(x = Year, y = Serovar, fill = Count)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = c(
      zero_color,
      low_color,
      "#C7E9F3",
      "#73BFE2",
      "#3182BD",
      high_color
    ),
    values = scales::rescale(
      c(0, 0.1, 5, 20, 50, 100, max_count)
    ),
    limits = c(0, max_count),
    name = "Count"
  ) +
  geom_text(
    data = subset(heatmap_data_complete, Count == 0),
    aes(label = "✗"),
    color = "#666666",
    size = 2.5
  ) +
  theme_minimal() +
  labs(
    title = "Heatmap of Salmonella Serovar Distribution by Year",
    subtitle = "Top 15 serovars among bovine isolates",
    x = "Year",
    y = "Serovar"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  ) +
  scale_x_continuous(breaks = all_years)

print(heatmap_plot)