suppressWarnings({
  suppressPackageStartupMessages({
    library(bnlearn)
    library(tidyverse)
    library(here)
    library(parallel)
  })
})

# --- Функция уведомлений ---
send_telegram <- function(message) {
  token   <- "8601597077:AAElWfcPY8iE2OWVokIxvgnH2yOxVyUeDtY"
  chat_id <- "8601597077" # Убедись, что пробелов в ID нет
  url <- paste0("https://api.telegram.org/bot", token, "/sendMessage")
  tryCatch({
    system(paste0('curl -s -X POST "', url, '" -d chat_id="', chat_id, '" -d text="', message, '"'), 
           ignore.stdout = TRUE)
  }, error = function(e) cat("Telegram failed\n"))
}

cat("--- Script 03b: Robust BN bootstrap ---\n")

# 1. Загрузка данных
full_counts <- readRDS(here("results", "full_count_matrix.rds"))
sig_genes   <- read.csv(here("results", "significant_genes.csv"), row.names = 1)
metadata    <- readRDS(here("results", "cleaned_metadata.rds"))

# 2. Подготовка матрицы (как в твоем коде)
deg_genes <- rownames(sig_genes)
available <- intersect(deg_genes, rownames(full_counts))
expr_mat  <- t(full_counts[available, ])
expr_mat  <- expr_mat[, apply(expr_mat, 2, var) > 0]
bn_data   <- as.data.frame(expr_mat)

# 3. Настройка бутстрепа
R_total <- 50             # Сколько всего хотим итераций
batch_size <- 5           # Сохраняем результат каждые 5 итераций
n_batches <- R_total / batch_size
n_cores <- max(1, detectCores() - 1)

dir.create(here("results", "bn_checkpoints"), showWarnings = FALSE)

send_telegram(sprintf("🚀 Старт BN: %d генов, R=%d, Cores=%d", ncol(bn_data), R_total, n_cores))

# Цикл по батчам (для надежности)
all_boot_list <- list()

for (b in 1:n_batches) {
  batch_file <- here("results", "bn_checkpoints", paste0("batch_", b, ".rds"))
  
  if (file.exists(batch_file)) {
    all_boot_list[[b]] <- readRDS(batch_file)
    cat("Батч", b, "загружен из файла\n")
  } else {
    cat("Считаем батч", b, "...\n")
    cl <- makeCluster(n_cores)
    clusterExport(cl, "bn_data")
    
    tryCatch({
      # Считаем маленькую порцию бутстрепа
      batch_res <- boot.strength(
        data = bn_data, R = batch_size, algorithm = "hc",
        algorithm.args = list(score = "bic-g", max.iter = Inf), # Снимаем лимит итераций
        cluster = cl
      )
      saveRDS(batch_res, batch_file)
      all_boot_list[[b]] <- batch_res
      send_telegram(sprintf("📊 Прогресс: выполнено %d/%d итераций", b * batch_size, R_total))
    }, error = function(e) {
      send_telegram(paste("❌ Ошибка в батче", b, ":", e$message))
    })
    
    stopCluster(cl)
  }
}

# 4. Объединяем и усредняем
cat("Собираем результаты...\n")
# Правильное объединение результатов boot.strength
final_boot_res <- all_boot_list[[1]]
if(length(all_boot_list) > 1) {
  for(i in 2:length(all_boot_list)) {
    final_boot_res$strength <- final_boot_res$strength + all_boot_list[[i]]$strength
    final_boot_res$direction <- final_boot_res$direction + all_boot_list[[i]]$direction
  }
  final_boot_res$strength <- final_boot_res$strength / length(all_boot_list)
  final_boot_res$direction <- final_boot_res$direction / length(all_boot_list)
}

# Итоговая сеть
cons_bn <- averaged.network(final_boot_res, threshold = 0.5)

# 5. Сохранение
saveRDS(cons_bn, here("results", "BN_FINAL_structure.rds"))
saveRDS(final_boot_res, here("results", "BN_FINAL_boot_res.rds"))

send_telegram(sprintf("✅ BN готова! Узлов: %d, Ребер: %d", 
                      length(nodes(cons_bn)), nrow(arcs(cons_bn))))