# PoLCA-R-case
# The present case describes the experience of R-usage for proceeding Latent Class Analysis during sociological research and further appliance of the results in SPSS analysis.

Анализ данных проводился в среде статистических вычислений R: A Language and Environment for Statistical Computing v. 4.5.3 с использованием пакета poLCA: Polytomous Variable Latent Class Analysis v. 1.6.0.2.

Стояла задача присвоить классы респондентам, принявшим участие в социологическом исследовании, посвященном их участию в реализации проекта "Демография" (N = 428). Распределение проводилось по вопросу "Почему Вы принимаете участие в проекте "Содействие занятости" в качестве преподавателя?". Вопрос предполагал возможность множественного выбора следующих вариантов ответа:
1. Финансовое вознаграждение
2. Интерес к участию в новом проекте
3. Перспективы для профессионального саморазвития
4. Желание помочь команде / руководителю
5. Перспективы развития для своей организации
6. Поиск потенциальных сотрудников
7. Желание реализовать свой потенциал / применить свои компетенции
8. Не было выбора

Данные были изначально представлены в формате .sav для их обработке в программе IBM SPSS Statistics v. 25.0. В этом виде они были перенесены в среду R для проведения латентного классового анализа. После его проведения и выявления оптимального количества классов (2 класса), данные классы были присвоены респондентам; итоговый массив был обратно транспонирован в .sav для проведения дальнейших расчетов по классам (линейный регрессионный анализ, сравнение средних по параметрам, корреляционный анализ и др.). Исходя из полученных результатов было решено номинировать полученные классы как "Респонденты, ориентированные на организацию" и "Респонденты, ориентированные на себя".

Проведенный анализ решил ряд исследовательских задач, поставленных в ходе подготовки научной статьи и показал статистически значимые различия между выявленными классами по таким параметрам как уровень лояльности, степень вовлеченности в проект, новообретенные компетенции и др.

Полный ход работы в среде R отражен ниже:

```install.packages("haven")
install.packages("poLCA")
install.packages("tidyverse")

library(haven)
library(poLCA)
library(tidyverse)

cat("\n=== Uploading data ===\n")
cat("Choose the file .sav\n")
data_full <- read_sav(file.choose())
cat("Data is uploaded. Size:", dim(data_full), "\n")

indicator_vars <- c("Q10_1", "Q10_2", "Q10_3", "Q10_4", "Q10_5", "Q10_6", "Q10_7", "Q10_8")

cat("\n=== Transforming to factors ===\n")
for (var in indicator_vars) {
  data_full[[var]] <- as.factor(data_full[[var]])
  cat("Variable", var, "Transformed to factors. Levels:", 
      length(levels(data_full[[var]])), "\n")
}

cat("\n=== Missing values processing ===\n")
n_before <- nrow(data_full)
data_clean <- na.omit(data_full[, indicator_vars])
n_after <- nrow(data_clean)

cat("Initial number of data points:", n_before, "\n")
cat("After removing missing data points:", n_after, "\n")
cat("Data points removed:", n_before - n_after, "\n")

if (n_after < 50) {
  stop("Too little data points after removing missing data points (less than 50). 
       Analysis is imposible.")
}

lca_formula <- as.formula(paste("cbind(", 
                                paste(indicator_vars, collapse = ", "), 
                                ") ~ 1"))

cat("\n=== Model formula ===\n")
cat(deparse(lca_formula), "\n")

cat("\n=== Finding optimal number of classes ===\n")

max_classes <- 6
results <- data.frame()

for (k in 1:max_classes) {
  cat("Model estimation to", k, "classes...")
  
    model <- tryCatch({
    poLCA(lca_formula, 
          data = data_clean, 
          nclass = k, 
          nrep = 10,
          maxiter = 5000,
          graphs = FALSE,
          verbose = FALSE)
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
    return(NULL)
  })
  
   if (!is.null(model)) {
       posterior <- model$posterior
    n <- nrow(posterior)
    C <- model$nclass
    entropy <- 1 - (sum(-posterior * log(posterior + 1e-10)) / (n * log(C)))
    
        results <- rbind(results, data.frame(
      Classes = k,
      LogLik = model$llik,
      AIC = model$aic,
      BIC = model$bic,
      Entropy = entropy,
      MinClassSize = min(model$P) * 100
    ))
    
    cat(" BIC =", round(model$bic, 1), 
        "Entropy =", round(entropy, 3),
        "MinClass =", round(min(model$P) * 100, 1), "%\n")
  } else {
    results <- rbind(results, data.frame(
      Classes = k,
      LogLik = NA, AIC = NA, BIC = NA, Entropy = NA, MinClassSize = NA
    ))
  }
}

cat("\n=== Table of model comparison ===\n")
print(results)

optimal_classes <- 2

cat("\n=== Choice of classes number ===\n")
cat("Classes are chosen:", optimal_classes, "\n")

cat("\n=== Final model estimation ===\n")
final_model <- poLCA(lca_formula, 
                     data = data_clean, 
                     nclass = optimal_classes, 
                     nrep = 20, 
                     maxiter = 5000, 
                     graphs = TRUE,
                     verbose = TRUE)

cat("\n=== Final model ===\n")
cat("Log-likelihood function:", final_model$llik, "\n")
cat("AIC:", final_model$aic, "\n")
cat("BIC:", final_model$bic, "\n")

cat("\nClasses sizes:\n")
for (i in 1:optimal_classes) {
  cat(sprintf("  Класс %d: %.1f%% (%d respondents)\n", 
              i, final_model$P[i] * 100, round(final_model$P[i] * final_model$Nobs)))
}

cat("\n=== Classes assumption ===\n")

classification_df <- data.frame(
  row_id = rownames(data_clean),
  assigned_class = final_model$predclass,
  stringsAsFactors = FALSE
)

posterior_probs <- final_model$posterior
for (c in 1:optimal_classes) {
  classification_df[[paste0("prob_class_", c)]] <- posterior_probs[, c]
}

classification_df$max_probability <- apply(posterior_probs, 1, max)

cat("Dataframe of classification is created", nrow(classification_df), "lines\n")
cat("Classes distribution in clear data:\n")
print(table(classification_df$assigned_class))

cat("\n=== Merging with initial data ===\n")

data_full$row_id <- as.character(1:nrow(data_full))
classification_df$row_id <- as.character(classification_df$row_id)

data_final <- left_join(data_full, classification_df, by = "row_id")

data_final$row_id <- NULL

cat("Size of final data:", dim(data_final), "\n")
cat("Number of data points with class assumption:", 
    sum(!is.na(data_final$assigned_class)), "\n")
cat("Number of missing data points (not included to LCA):", 
    sum(is.na(data_final$assigned_class)), "\n")

cat("\n=== Saving results ===\n")

write_sav(data_final, "complete_data_with_classes.sav")
cat("SPSS file is saved: complete_data_with_classes.sav\n")

write.csv(results, "model_comparison.csv", row.names = FALSE)
cat("Table of models comparison is saved: model_comparison.csv\n")

save(final_model, file = "lca_model.RData")
cat("Model is saved: lca_model.RData\n")
