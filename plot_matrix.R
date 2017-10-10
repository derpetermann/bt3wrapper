plot_matrix <- function(p_matrix, margin = 0, label = 'X', arrow_size = .6, arrow_width = 1,
                        scale_factor = 5, add = FALSE) {
  graph <- graph_from_adjacency_matrix(adjmatrix = p_matrix, weighted = TRUE, mode = "directed")
  plot.igraph(graph,ylim = c(-1,0), 
              xlab = label, asp = 1, margin = margin, add = add,
              vertex.shape = 'square', vertex.color=c('black','white'),
              edge.color = 'black', vertex.size=20,
              vertex.label = NA, layout = layout_on_grid,
              edge.curved = .4, edge.width = E(graph)$weight*scale_factor,
              edge.arrow.size = arrow_size, edge.arrow.width = arrow_width,
              edge.loop.angle = -80)}
