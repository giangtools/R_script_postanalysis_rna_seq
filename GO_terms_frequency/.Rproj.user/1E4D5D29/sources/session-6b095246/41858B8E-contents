# Cài đặt các thư viện cần thiết
library(tidyverse)

# Đọc dữ liệu từ file TSV
file_path <- "D:/Backup_data_lab/3192_analysis/salmon_results_with_go_terms.tsv"
df <- read_tsv(file_path)

# Chuyển đổi cột GO_terms từ chuỗi ký tự thành danh sách
df <- df %>%
  mutate(GO_terms = str_remove_all(GO_terms, "[\\[\\]' ]") %>%
           str_split(","))

# Tính tần số xuất hiện của các GO terms
go_freq <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(freq = n(), .groups = 'drop')

# Tính độ lớn của mỗi GO term
go_size <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(size = sum(mean, na.rm = TRUE) / n(), .groups = 'drop')
print(go_size)
# Kết hợp tần số và độ lớn của các GO terms
go_data <- go_freq %>%
  inner_join(go_size, by = "GO_terms")

# Đọc file TSV có cột namespace
namespace_file_path <- "d:/R/Script/GO_terms_frequency/go_term_counts_namespace.tsv"
namespace_df <- read_tsv(namespace_file_path)

# Kết hợp dữ liệu tần số, độ lớn và namespace
go_data <- go_data %>%
  left_join(namespace_df, by = c("GO_terms" = "GO_term"))

print(go_data)

# Vẽ biểu đồ chấm tròn, gom nhóm theo namespace
ggplot(go_data, aes(x = GO_terms, y = size, size = freq)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "GO Term Sizes", x = "GO Terms", y = "Average TPM", size = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ namespace, scales = "free_x")

# Lọc 20 GO terms có độ lớn lớn nhất
top_20_go <- go_data %>%
  arrange(desc(size)) %>%
  head(20)

print(top_20_go)

# Vẽ biểu đồ chấm tròn cho top 20 GO terms, gom nhóm theo namespace
ggplot(top_20_go, aes(x = reorder(GO_terms, -size), y = size, size = freq)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Top 20 GO Term Sizes", x = "GO Terms", y = "Average TPM", size = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ namespace, scales = "free_x")

# Lưu biểu đồ dưới dạng file
ggsave("go_term_sizes.png")
