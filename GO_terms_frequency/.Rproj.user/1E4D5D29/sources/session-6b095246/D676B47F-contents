# Cài đặt các thư viện cần thiết
library(tidyverse)

# Đọc dữ liệu từ file TSV
#file_path <- "D:/Backup_data_lab/3192_analysis/salmon_results_with_go_terms.tsv"
file_path <- "C:/Users/ADMIN/TGiang/GD_63_Postanalysis/Test/salmon_results_with_go_terms.tsv"
df <- read_tsv(file_path)

# Chuyển đổi cột GO_terms từ chuỗi ký tự thành danh sách
df <- df %>%
  mutate(GO_terms = map(GO_terms, ~str_remove_all(., "[\\[\\]' ]") %>% str_split(",") %>% unlist))
print(df)
print(df)
# Đếm số lượng các gene theo các mức biểu hiện
expression_counts <- df %>%
  group_by(expression_level) %>%
  summarise(count = n())
# Hiển thị kết quả
print(expression_counts)

ggplot(expression_counts, aes(x = expression_level, y = count)) +
     geom_bar(stat = "identity", fill = "skyblue") +
     labs(title = "Gene Expression Levels",
          x = "Expression Level",
          y = "Count") +
  theme_minimal()

# Tính tần số xuất hiện của các GO terms
go_freq <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(freq = n())
print(go_freq)

go_freq_gene <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(freq = n(), genes = list(unique(gene_id))) %>%
  mutate(genes = map_chr(genes, ~paste(., collapse = ",")))

print(go_freq_gene)
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

#Kết hợp hai dataframe lại với nhau
combined_df <- go_data %>%
  inner_join(df_ns, by = "GO_terms")
print(combined_df)
# Lọc ra 20 go_term có tầng số cao nhất
top_20_freq <- combined_df %>%
  arrange(desc(freq.x)) %>%
  head(20)
print(top_20_freq)

# Vẽ biểu đồ cột 20 CỤM LỚN NHẤT sử dụng ggplot2 với tên 
# Vẽ biểu đồ cột 20 CỤM LỚN NHẤT sử dụng ggplot2 với tên 
# Vẽ biểu đồ cột 20 CỤM LỚN NHẤT sử dụng ggplot2 với tên 
# Vẽ biểu đồ cột 20 CỤM LỚN NHẤT sử dụng ggplot2 với tên 
# Lọc dữ liệu để lấy 20 GO_terms có size lớn nhất cho mỗi namespace
top_20_go_terms <- combined_df %>%
  group_by(namespace) %>%
  arrange(namespace, desc(size)) %>%
  slice_head(n = 20) %>%
  ungroup()
print(top_20_go_terms)
ggplot(top_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ namespace, scales = "free_x") +
  geom_text(aes(label = freq.x), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal() +
  labs(title = "Top 20 GO Terms by Size in Each Namespace",
       x = "GO Terms",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




# Vẽ biểu đồ cột 20 CỤM NHỎ NHẤT sử dụng ggplot2 với tên 
# Vẽ biểu đồ cột 20 CỤM NHỎ NHẤT sử dụng ggplot2 với tên 
# Vẽ biểu đồ cột 20 CỤM NHỎ NHẤT sử dụng ggplot2 với tên 
bottom_20_go_terms <- combined_df %>%
  group_by(namespace) %>%
  arrange(namespace, size) %>%
  slice_head(n = 20) %>%
  ungroup()

# Vẽ biểu đồ cột sử dụng ggplot2
# Giới hạn chiều dài của các tên để dễ hiển thị trên trục x
bottom_20_go_terms <- bottom_20_go_terms %>%
  mutate(name = ifelse(nchar(name) > 30, paste0(substr(name, 1, 27), "..."), name))
print(n=60, bottom_20_go_terms)
ggplot(bottom_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ namespace, scales = "free_x") +
  geom_text(aes(label = freq.x), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal() +
  labs(title = "Bottom 20 GO Terms by Size in Each Namespace",
       x = "Name",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
print(bottom_20_go_terms)



