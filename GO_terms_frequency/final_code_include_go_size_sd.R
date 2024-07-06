# Cài đặt các thư viện cần thiết
library(tidyverse)
library(dplyr)
# Đọc dữ liệu từ file TSV

file_path <- "C:/Users/ADMIN/TGiang/GD_63_Postanalysis/Res/go_terms_python_merged/salmon_results_with_go_terms.tsv"
df <- read_tsv(file_path)

# Chuyển đổi cột GO_terms từ chuỗi ký tự thành danh sách
df <- df %>%
  mutate(GO_terms = map(GO_terms, ~str_remove_all(., "[\\[\\]' ]") %>% str_split(",") %>% unlist()))
print(df)


# Giả định rằng df đã được tạo từ trước
df <- df %>%
  mutate(expression_level = case_when(
    mean < 0.5 ~ "below cutoff",
    mean >= 0.5 & mean <= 10 ~ "low",
    mean > 10 & mean <= 1000 ~ "medium",
    mean > 1000 ~ "high"
  ))

# Kiểm tra kết quả
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
write_tsv(go_freq, "D:/R/Script/GO_terms_frequency/go_freq.tsv")
#Chạy thêm Đoạn code filter.py
file_path_ns <- "go_term_counts_namespace_name.tsv"
df_ns <- read_tsv(file_path_ns)
print(df_ns)

go_freq_gene <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(freq = n(), genes = list(unique(gene_id))) %>%
  mutate(genes = map_chr(genes, ~paste(., collapse = ",")))
print(go_freq_gene)
##########################################
# Tạo một dataframe mới chứa thông tin gene_id và product từ df
gene_product_df <- df %>%
  select(gene_id, product)

# Tạo cột mới trong go_freq_gene chứa "genes:product"
go_freq_gene <- go_freq_gene %>%
  mutate(genes_products = map_chr(genes, function(gene_list) {
    gene_ids <- unlist(str_split(gene_list, ","))
    products <- gene_product_df %>%
      filter(gene_id %in% gene_ids) %>%
      pull(product)
    paste(paste(gene_ids, products, sep = ":"), collapse = ",")
  }))

# Hiển thị dataframe go_freq_gene sau khi thêm cột mới
print(go_freq_gene)
write_csv(go_freq_gene,"C:/Users/ADMIN/TGiang/GD_63_Postanalysis/Test/go_freq_gene.csv")
##########################################

# Tính trung bình và độ lệch chuẩn của các giá trị gene cho mỗi GO_terms
go_stats <- df %>%
  unnest(GO_terms) %>%
  filter(GO_terms != "") %>%
  group_by(GO_terms) %>%
  summarise(size = mean(mean), sd = sd(mean))
print(go_stats)

# Kết hợp tần số và độ lớn của các GO terms
go_data <- go_freq %>%
  inner_join(go_stats, by = "GO_terms")
print(go_data)

# Giả sử df_ns là một dataframe có thông tin namespace của mỗi GO_terms
# Kết hợp với dữ liệu namespace
combined_df <- go_data %>%
  inner_join(df_ns, by = "GO_terms")
##Lưu ý bổ sunng 
combined_df <- combined_df %>%
  inner_join(go_freq_gene %>% select(GO_terms, genes, genes_products), by = "GO_terms")
print(combined_df)

# Lọc ra 20 go_term có tầng số cao nhất
top_20_freq <- combined_df %>%
  arrange(desc(freq.x)) %>%
  head(20)
print(top_20_freq)

# Vẽ biểu đồ cột với thanh lỗi (error bars)
top_20_go_terms <- combined_df %>%
  group_by(namespace) %>%
  arrange(namespace, desc(size)) %>%
  slice_head(n = 10) %>%
  ungroup()
######
top_20_go_terms <- combined_df %>%
  filter(!is.na(namespace)) %>%
  group_by(namespace) %>%
  arrange(namespace, desc(size)) %>%
  slice_max(size, n = 20) %>%
  ungroup()

print(top_20_go_terms)
######
ggplot(top_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = size - sd, ymax = size + sd), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ namespace, scales = "free_x") +
  geom_text(aes(label = freq.x), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal() +
  labs(title = "Top 20 GO Terms by Size in Each Namespace",
       x = "GO Terms",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


top_20_go_terms <- top_20_go_terms %>%
  left_join(select(combined_df, GO_terms, genes, gene), by = "GO_terms")
print(top_20_go_terms)
write_csv(top_20_go_terms,"C:/Users/ADMIN/TGiang/GD_63_Postanalysis/Test/top_20_go_terms.csv")


# Vẽ biểu đồ cột 20 CỤM NHỎ NHẤT sử dụng ggplot2
bottom_20_go_terms <- combined_df %>%
  group_by(namespace) %>%
  arrange(namespace, size) %>%
  slice_head(n = 20) %>%
  ungroup()
#####
bottom_20_go_terms <- combined_df %>%
  group_by(namespace) %>%
  arrange(namespace, size) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  filter(!is.na(namespace)) %>%  # Loại bỏ các dòng có giá trị NA trong namespace
  mutate(name = ifelse(nchar(name) > 30, paste0(substr(name, 1, 27), "..."), name))
#####
# Giới hạn chiều dài của các tên để dễ hiển thị trên trục x
bottom_20_go_terms <- bottom_20_go_terms %>%
  mutate(name = ifelse(nchar(name) > 30, paste0(substr(name, 1, 27), "..."), name))
print(n=60, bottom_20_go_terms)

ggplot(bottom_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = size - sd, ymax = size + sd), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ namespace, scales = "free_x") +
  geom_text(aes(label = freq.x), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal() +
  labs(title = "Bottom 20 GO Terms by Size in Each Namespace",
       x = "Name",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(n=60,bottom_20_go_terms)
bottom_20_go_terms <- bottom_20_go_terms %>%
  left_join(select(combined_df, GO_terms, genes), by = "GO_terms")
write_csv(bottom_20_go_terms,"C:/Users/ADMIN/TGiang/GD_63_Postanalysis/Test/bottom_20_go_terms.csv")


# Thêm cột namespace và name vào dataframe go_freq_gene
go_freq_gene_name <- go_freq_gene %>%
  left_join(select(combined_df, GO_terms, namespace, name, genes), by = "GO_terms")
print(go_freq_gene_name)
write_csv(go_freq_gene_name,"C:/Users/ADMIN/TGiang/GD_63_Postanalysis/Test/go_freq_gene_name.csv")



library(patchwork)

# Vẽ biểu đồ top_20_go_terms
p1 <- ggplot(top_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = size - sd, ymax = size + sd), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ namespace, scales = "free_x") +
  geom_text(aes(label = freq.x), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal() +
  labs(title = "Top 20 GO Terms by Size in Each Namespace",
       x = "GO Terms",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Vẽ biểu đồ bottom_20_go_terms
p2 <- ggplot(bottom_20_go_terms, aes(x = reorder(name, -size), y = size, fill = namespace)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = size - sd, ymax = size + sd), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ namespace, scales = "free_x") +
  geom_text(aes(label = freq.x), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal() +
  labs(title = "Bottom 20 GO Terms by Size in Each Namespace",
       x = "GO Terms",
       y = "Size",
       fill = "Namespace") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Hiển thị hai biểu đồ cùng nhau với cùng tỉ lệ scale
p1 + p2 + plot_layout(ncol = 2)

