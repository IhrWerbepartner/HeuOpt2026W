#file_name <- "/home/johannes/repos/HeuOpt2026W/ex1/instances/test/10000/n10000_k3000_m3_s20011.txt"
#file_name <- "/home/johannes/repos/HeuOpt2026W/ex1/instances/competition/n5000_k1500_m2_s30003.txt"
# file_name <- "/home/johannes/repos/HeuOpt2026W/ex1/instances/competition/n5000_k1500_m2_s30002.txt"
file_name <- "/home/johannes/repos/HeuOpt2026W/ex1/instances/competition/n5000_k1500_m2_s30001.txt"
#file_name <- "/home/johannes/repos/HeuOpt2026W/ex1/instances/test/10000/n10000_k3000_m3_s20020.txt"


all_values <- integer()

con <- file(file_name, open = "r")
on.exit(close(con), add = TRUE)

for (line in readLines(con)) {
  if (!startsWith(line, "#")) {
    words <- strsplit(trimws(line), "\\s+")[[1]]
    words <- words[nzchar(words)]
    all_values <- c(all_values, as.integer(words))
  }
}

if (length(all_values) < 5) {
  stop(sprintf("invalid instance file: %s", file_name))
}

n <- all_values[1]
edge_count <- all_values[2]
k <- all_values[3]
m <- all_values[4]

expected_values <- 4 + (n * m) + (edge_count * (3 + m))

if (length(all_values) != expected_values) {
  stop(sprintf(
    "invalid number of values in RCkST instance: %d, expected: %d",
    length(all_values), expected_values
  ))
}

resources_end_index <- 4 + n * m

resources <- t(matrix(
  all_values[5:resources_end_index],
  nrow = m,
  ncol = n
))

edges <- t(matrix(
  all_values[(resources_end_index + 1):length(all_values)],
  nrow = 3 + m,
  ncol = edge_count
))

# vertices
par(mfrow = c(m, 1))
for (i in 1:m) {
  hist.data =  hist(resources[, i], breaks = "FD", main = sprintf("Resource availability %d", i))
}
par(mfrow = c(1, 1))

# edge weights
hist.data =  hist(edges[, 3], breaks = "FD", main = "Edge Weights")

# edges
par(mfrow = c(m, 1))
for (i in 1:m) {
  hist.data =  hist(edges[, i + 3], breaks = "FD", main = sprintf("Resources consumption %d", i))
}
par(mfrow = c(1, 1))

plot(edges[, 3], edges[, 4])
plot(edges[, 3], edges[, 5])
#plot(edges[, 3], edges[, 6])

vertex_degrees <- table(edges[, 1:2])
hist(vertex_degrees, breaks = "FD")

# vertex degree vs resource availabilty
plot(vertex_degrees[, 1], resources[, 1])

#hist(edges[, 6], breaks = "FD")
library(robustbase)
weight_value <- edges[, c(3, 5)]
weight_value.rob_mcd <- covMcd(weight_value)

best <- rep(2, length(weight_value[,1]))
best[weight_value.rob_mcd$best] <- 3

plot(weight_value, col = best, xlab = "Weight", ylab = "Resource Cost", main = "Weight vs Cost with outlier detection")
legend(x = "topleft", legend = c("Discarded set",
                                 "Optimal Subset"),
       pch = c(1, 1),
       col = c("red", "green"))

