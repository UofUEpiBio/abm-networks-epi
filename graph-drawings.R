library(igraph)
library(netplot)
library(ergm)

set.seed(2133)

n <- 50
x <- sample.int(2, n, TRUE)

myplot <- function(g) {
  nplot(
    g,
    edge.line.breaks = 1,
    edge.color = "darkgray",
    vertex.size.range = c(.05, .05),
    vertex.nsides = 30,
    vertex.label = NA
  ) |>
  set_vertex_gpar(
    element = "core",
    fill = lapply(c("steelblue", "tomato")[x], \(i) {
      radialGradient(c("white", i), cx1=.8, cy1=.8, r1=0)
      }))
}

# Scale-free network
g_sf <- sample_pa(n, m = 1, directed = FALSE, power = 1.5)

# Erdos-renyi graph with same density as g
g_er <- igraph::sample_gnm(n = vcount(g_sf), m = ecount(g_sf))

# Ergm simulation
g_er <- igraph::set_vertex_attr(g_er, "group", value = x)
y <- intergraph::asNetwork(g_er)

g_ergm <- simulate_formula(y ~ edges + nodematch("group"), coef = c(-3, 5), constraints = ~ edges)

# Degree sequence preserving
g_degseq <- igraph::rewire(
  intergraph::asIgraph(g_ergm),
  keeping_degseq(niter = ecount(g_sf)*20)
) 

# Small world
g_sw <- igraph::sample_smallworld(
  n,
  dim = 1,
  nei = ceiling(mean(sna::degree(g_ergm, gmode="graph", cmode = "indegree")))/2,
  p   = 0.1
  )

g_sw2 <- igraph::sample_smallworld(
  n,
  dim = 1,
  nei = ceiling(mean(sna::degree(g_ergm, gmode="graph", cmode = "indegree")))/2,
  p   = 0.2
  )

np_sf <- myplot(g_sf) 

np_er <- myplot(g_er) 

np_ergm <- myplot(g_ergm) 

np_degseq <- myplot(g_degseq) 

np_sw <- myplot(g_sw)
np_sw2 <- myplot(g_sw2)

np <- gridExtra::grid.arrange(
  grobs = list(np_ergm, np_sf, np_er, np_degseq, np_sw, np_sw2),
  nrow = 2, 
  )

ggplot2::ggsave(
  filename = "graph-drawings.svg",
  plot     = np,
  width    = 10,
  height = 5, dpi = 300
  )
