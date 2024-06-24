# Cài đặt các thư viện cần thiết
library(ggplot2)
library(dplyr)

# Đọc dữ liệu từ file TSV
df_cv <- read.csv("d:/Backup_data_lab/3192_analysis/salmon_results_with_mean_sd_sv.tsv", sep = '\t', header = TRUE)

# Kiểm tra dữ liệu
print(head(df_cv))

# Vẽ biểu đồ chấm tròn với mean và SV
ggplot(df_cv, aes(x = mean, y = SV)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Scatter Plot of Mean vs SV",
       x = "Mean",
       y = "SV") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Lưu biểu đồ chấm tròn dưới dạng file
ggsave("mean_sv_dot_plot.png")
