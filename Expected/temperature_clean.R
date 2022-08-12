library(tidyverse)

# Load in data
temp_data <- read_csv("../WilliamData/multinomial_model/temperatures.csv")

temp_data <- temp_data %>%
  mutate(iso3 = str_sub(ISO3, start = -3),
         month = rep(1:12, times = nrow(temp_data) / 12),
         temperature = `Temperature - (Celsius)`,
         year = Year) %>%
  select(iso3, year, month, temperature) %>% 
  arrange(iso3, year, month)

# Niue is the only country without monthly historical acm data
# without temperatures in the file William uploaded
# Downloaded data from https://climateknowledgeportal.worldbank.org/download-data
niue_data <- read_csv("Niue_temperature_data.csv", skip = 1)
niue_temp <- data.frame(iso3 = "NIU",
                        year = 2020,
                        month = 1:12,
                        temperature = unname(unlist(niue_data[1, 3:14])))
temp_data <- rbind(temp_data, niue_temp) %>% 
  arrange(iso3, year, month)

save(temp_data, file = "temperature_data_clean.RData")
