# Вставь в начало любого скрипта
send_telegram <- function(message) {
  token   <- "ТВОЙ_ТОКЕН"
  chat_id <- "5095424200"
  
  url <- paste0("https://api.telegram.org/bot", token, "/sendMessage")
  
  cmd <- paste(
    "curl -s -X POST",
    shQuote(url),
    "-d", shQuote(paste0("chat_id=", chat_id)),
    "--data-urlencode", shQuote(paste0("text=", message))
  )
  
  result <- system(cmd, intern = TRUE)
  cat(result, sep = "\n")
}

# Перед запуском
send_telegram("🚀 Bayes_26: BN bootstrap started (R=50, 1116 genes)")

time_start <- proc.time()

boot_res <- boot.strength(...)

time_end <- proc.time()
elapsed  <- round((time_end - time_start)["elapsed"] / 60, 1)

# После завершения
n_edges <- nrow(boot_res[boot_res$strength > 0.5, ])
send_telegram(paste0(
  "✅ Bayes_26: BN готова!\n",
  "⏱ Время: ", elapsed, " мин\n",
  "🔗 Рёбер (strength>0.5): ", n_edges
))

# При ошибке — обернуть в tryCatch
tryCatch({
  boot_res <- boot.strength(...)
  send_telegram("✅ BN completed!")
}, error = function(e) {
  send_telegram(paste0("❌ BN crashed: ", e$message))
})