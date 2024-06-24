# Cài đặt các thư viện cần thiết
library(tidyverse)

# Đọc dữ liệu từ file TSV
file_path <- "D:/Backup_data_lab/3192_analysis/salmon_results_with_go_terms.tsv"
df <- read_tsv(file_path)

file_path_ns <- "D:/R/Script/GO_terms_frequency/go_term_counts_namespace_filtered_1.tsv"
df_ns <- read_tsv(file_path_ns)

# Chuyển đổi cột GO_terms từ chuỗi ký tự thành danh sách
df <- df %>%
  mutate(GO_terms = map(GO_terms, ~str_remove_all(., "[\\[\\]' ]") %>% str_split(",") %>% unlist))

# Tính tần số xuất hiện của các GO terms
go_freq <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(freq = n())

# Tính độ lớn của mỗi GO term
go_size <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(size = sum(mean) / n())

# Kết hợp tần số và độ lớn của các GO terms
go_data <- go_freq %>%
  inner_join(go_size, by = "GO_terms")
print(go_data)
print(df_ns)

# Kết hợp hai dataframe lại với nhau
combined_df <- go_data %>%
  inner_join(df_ns, by = "GO_terms")
print(combined_df)

# Lọc dữ liệu để lấy 20 GO_terms có size nhỏ nhất cho mỗi namespace
bottom_20_go_terms <- combined_df %>%
  group_by(namespace) %>%
  arrange(namespace, size) %>%
  slice_head(n = 20) %>%
  ungroup()

# Vẽ biểu đồ cột sử dụng ggplot2
# Vẽ biểu đồ cột sử dụng ggplot2 và scale lại giá trị trục tung
ggplot(bottom_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ namespace, scales = "free_x") +
  theme_minimal() +
  labs(title = "Top 20 GO Terms by Size in Each Namespace",
       x = "Name",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(limits = c(0, 30), labels = scales::comma)
print(bottom_20_go_terms, n=100)

size_1 <- bottom_20_go_terms$GO_terms
name_sp <- bottom_20_go_terms$namespace
nam <- bottom_20_go_terms$name
new_df <- data.frame(GO_terms = size_1, name = nam)
print(new_df)
library(ggplot2)
# Vẽ biểu đồ cột
ggplot(bottom_20_go_terms, aes(x = reorder(GO_terms, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +  # Vẽ cột
  facet_wrap(~ namespace, scales = "free_x") +
  theme_minimal() +  # Sử dụng giao diện minimal
  labs(title = "Top 20 GO Terms by Size in Each Namespace",  # Đặt tiêu đề
       x = "GO Terms",  # Đặt nhãn trục hoành
       y = "Size",  # Đặt nhãn trục tung
       fill = "Namespace") +  # Đặt nhãn fill
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Xoay nhãn trục hoành

print(bottom_20_go_terms)

# Lưu biểu đồ với kích thước lớn hơn
ggsave("bottom_20_go_terms_plot.png", width = 10, height = 6)


print(bottom_20_go_terms)
write_tsv(new_df, "C:/Users/ADMIN/TGiang/Powerpoint_mollab/bottom_gene.tsv")
# Lưu biểu đồ dưới dạng file
ggsave("bottom_20_go_term_sizes.png")
