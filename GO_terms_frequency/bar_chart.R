# Cài đặt các gói cần thiết
install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")

# Nạp các gói
library(dplyr)
library(ggplot2)
library(stringr)

# Đọc dữ liệu từ file CSV
file_path <- "D:/R/Script/GO_terms_frequency/id_goterms.tsv"  # Thay 'path/to/your/file.csv' bằng đường dẫn tới file CSV của bạn
df <- read.csv(file_path, sep = "\t", stringsAsFactors = FALSE)

# Chuyển đổi cột GO_terms từ chuỗi thành danh sách
df$GO_terms <- strsplit(df$GO_terms, ", ")

# Tính tần số xuất hiện của mỗi mã GO
all_go_terms <- unlist(df$GO_terms)
go_term_counts <- as.data.frame(table(all_go_terms))
colnames(go_term_counts) <- c("GO_term", "Frequency")

# Vẽ biểu đồ cột
ggplot(go_term_counts, aes(x = reorder(GO_term, -Frequency), y = Frequency)) +
  geom_bar(stat = "identity") +
  xlab("GO Term") +
  ylab("Frequency") +
  ggtitle("Frequency of GO Terms") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()
