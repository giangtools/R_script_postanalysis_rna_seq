library(tidyverse)
library(ontologyIndex)
go_obo_file <- "D:/Python/exp/annotation_handling/script/go.obo"
go_ontology <- get_ontology(go_obo_file)
# Extract GO term names
go_term_names <- data.frame(
  GO_terms = names(go_ontology$name),
  name = unname(go_ontology$name)
)
file_path <- "D:/Backup_data_lab/3192_analysis/salmon_results_with_go_terms.tsv"
df <- read_tsv(file_path)

# Convert GO_terms column to a list
df <- df %>%
  mutate(GO_terms = map(GO_terms, ~str_remove_all(., "[\\[\\]' ]") %>% str_split(",") %>% unlist()))
# Calculate GO term frequency
go_freq <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(freq = n(), .groups = 'drop')

# Calculate GO term size
go_size <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(size = sum(mean) / n(), .groups = 'drop')

# Combine frequency and size
go_data <- go_freq %>%
  inner_join(go_size, by = "GO_terms")
# Merge GO term names
go_data <- go_data %>%
  left_join(go_term_names, by = "GO_terms")
top_20_go <- go_data %>%
  arrange(desc(size)) %>%
  head(20)

ggplot(top_20_go, aes(x = reorder(name, -size), y = size, size = freq)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Top 20 GO Term Sizes", x = "GO Term Names", y = "Average TPM", size = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(top_20_go, aes(x = reorder(name, -size), y = size)) +
  geom_point(aes(size = size), alpha = 0.7) +
  theme_minimal() +
  labs(title = "Top 20 GO Term Sizes", x = "GO Term Names", y = "Average TPM", size = "Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Vẽ biểu đồ dựa trên namespace: 
# Tạo biến cho mỗi namespace
# Tạo biến cho mỗi namespace
go_data <- go_data %>%
  mutate(namespace = factor(go_term_namespace, levels = c("biological_process", "molecular_function", "cellular_component")))

# Sắp xếp lại go_data theo namespace và size
go_data <- go_data %>%
  arrange(namespace, desc(size))

# Vẽ biểu đồ
ggplot(go_data, aes(x = reorder(name, -size), y = size, color = namespace, size = freq)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Top 20 GO Term Sizes", x = "GO Term Names", y = "Average TPM", size = "Frequency", color = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("biological_process" = "blue", "molecular_function" = "red", "cellular_component" = "green"))
