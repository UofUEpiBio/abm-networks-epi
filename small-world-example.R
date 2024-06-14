

library(epiworldR)
library(data.table)

n      <- 2e3
p_r    <- 1/3
C_rate <- 7
R0     <- 1.5
p_t <- p_r/(C_rate/R0 + p_r - 1)

x <- ModelSIR(
  name              = "as", 
  prevalence        = 10/n,
  transmission_rate = p_t,
  recovery_rate     = p_r
  )

agents_smallworld(x, n, k = 14, p = 0.1, d = FALSE)

run_multiple(x, ndays = 50, nsims = 500, seed = 111, nthreads = 6, saver = make_saver("reproductive"))

ans <- run_multiple_get_results(x)
rt <- plot(ans$reproductive)

rt <- data.table(rt)

rt_daily <- rt[, mean(avg, na.rm = TRUE), by = date]

library(ggplot2)
ggplot(rt, aes(x = date, y = avg)) +
  geom_boxplot(aes(group = date), fill = "#63c588", colour = "#63c588") +
  geom_hline(yintercept = 1, linetype = 2, colour = "black", lwd = 1) +
  geom_line(data = rt_daily, aes(x = date, y = V1), colour = "#e72020", lwd = 1) +
  annotate("text", x = 2, y = 1.1, label = "R0 = 1", size = 6) +
  # Adding a legend of the colors used in geom_line
  annotate("text", x = 45, y = .3, label = "Mean R", size = 6, colour = "#a51919", ) +
  scale_y_log10() +
  labs(x = "Day", y = "Reproductive number\n(log10)")

# Saving the plot in a 4 x 4 inches png file
ggsave("data/small-world-example.png", width = 6, height = 6, units = "in", dpi = 300)

