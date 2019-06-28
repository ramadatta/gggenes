#' A 'ggplot2' geom to draw genes as arrows
#'
#' `geom_gene_arrow()` draws genes as arrows, allowing gene maps to be drawn.
#'
#' This geom draws genes as arrows along a horizontal line representing the
#' molecule. The start and end locations of the gene are expressed with the
#' `xmin` and `xmax` aesthetics, while the molecule can be specified with the
#' `y` aesthetic. Optionally, an additional `forward` aesthetic can be used to
#' reverse the orientation of some or all genes from that implied by `xmin` and
#' `xmax`.
#'
#' Unless the plot is faceted with a free x scale, all the molecules will share
#' a common x axis. This means that if the locations are very different across
#' different molecules, the genes might appear very small and squished together
#' with a lot of unnecessary empty space. To get around this, either facet the
#' plot with `scales = "free_x"`, or normalise the gene locations if their
#' exact locations are not important.
#'
#' See `make_alignment_dummies()` for a method to align genes between molecules.
#'
#' @section Aesthetics:
#'
#' - xmin,xmax (start and end of the gene; will be used to determine gene
#' orientation)
#' - y (molecule)
#' - forward (if any value that is not TRUE, or coercible to TRUE, the gene
#' arrow will be drawn in the opposite direction to that determined by `xmin`
#' and `xmax`)
#' - alpha
#' - colour
#' - fill
#' - linetype
#' - size
#'
#' @param mapping,data,stat,position,na.rm,show.legend,inherit.aes,... As
#' standard for ggplot2.
#' @param arrowhead_width `grid::unit()` object giving the width of the
#' arrowhead.  Defaults to 4 mm. If the gene is drawn smaller than this width,
#' only the arrowhead will be drawn, compressed to the length of the gene.
#' @param arrowhead_height `grid::unit()` object giving the height of the
#' arrowhead.  Defaults to 4 mm.
#' @param height `grid::unit()` object giving the height of the body
#' of the arrow. Defaults to 3 mm.
#'
#' @examples
#'
#' ggplot2::ggplot(example_genes, ggplot2::aes(xmin = start, xmax = end,
#'                                             y = molecule, fill = gene)) +
#' geom_gene_arrow() +
#' ggplot2::facet_wrap(~ molecule, scales = "free")
#'
#' @seealso [theme_genes()], [make_alignment_dummies()], [geom_gene_label()]
#'
#' @export
geom_gene_arrow <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  arrowhead_width = grid::unit(4, "mm"),
  arrowhead_height = grid::unit(4, "mm"),
  ...
) {
  ggplot2::layer(
    geom = GeomGeneArrow, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      arrowhead_width = arrowhead_width,
      arrowhead_height = arrowhead_height,
      ...
    )
  )
}

#' GeomGeneArrow
#' @noRd
GeomGeneArrow <- ggplot2::ggproto("GeomGeneArrow", ggplot2::Geom,
  required_aes = c("xmin", "xmax", "y"),
  default_aes = ggplot2::aes(
    forward = TRUE,
    alpha = 1,
    colour = "black",
    fill = "white",
    linetype = 1,
    size = 0.3
  ),
  draw_key = function(data, params, size) {
    grid::rectGrob(
      width = grid::unit(1, "npc") - grid::unit(1, "mm"),
      height = grid::unit(1, "npc") - grid::unit(1, "mm"),
      gp = grid::gpar(
        col = data$colour,
        fill = ggplot2::alpha(data$fill, data$alpha),
        lty = data$linetype,
        lwd = data$size * ggplot2::.pt
      )
    )
  },

  setup_data = function(data, params) {

    # Reverse non-forward genes
    if ("forward" %in% names(data)) {
      data[data$forward == 1, c("xmin", "xmax")] <- 
        data[data$forward == 1, c("xmax", "xmin")]
    }

    # Set arrow body height
    data$height <- data$height %||%
      params$height %||% (resolution(data$y, FALSE) * 0.9)

    transform(data,
      ymin = y - height / 2, ymax = y + height / 2, height = NULL
    )
  },

  draw_panel = function(
    data,
    panel_scales,
    coord,
    arrowhead_width,
    arrowhead_height,
    height
  ) {

    print(data)

    data <- data.frame(
      x = as.vector(rbind(
        data$xmin,
        data$xmin,
        data$xmin + (0.8 * (data$xmax - data$xmin)),
        data$xmin + (0.8 * (data$xmax - data$xmin)),
        data$xmax,
        data$xmin + (0.8 * (data$xmax - data$xmin)),
        data$xmin + (0.8 * (data$xmax - data$xmin)),
        data$xmin
      )),
      y = as.vector(rbind(
        data$ymin,
        data$ymax,
        data$ymax,
        data$ymax + (0.2 * (data$ymax - data$ymin)),
        data$ymin + (0.5 * (data$ymax - data$ymin)),
        data$ymin - (0.2 * (data$ymax - data$ymin)),
        data$ymin,
        data$ymin
      )),
      fill = rep(data$fill, each = 8),
      # alpha = rep(data$alpha, each = 8),
      size = rep(data$size, each = 8),
      linetype = rep(data$linetype, each = 8),
      group = rep(1:(nrow(data)), each = 8),
      row.names = 1:(nrow(data) * 8)
    )

    GeomPolygon$draw_panel(data, panel_scales, coord)
  }
)

#' @importFrom grid makeContent
#' @export
makeContent.genearrowtree <- function(x) {

  data <- x$data
  print(data)

  # Prepare grob for each gene
  grobs <- lapply(1:nrow(data), function(i) {

    gene <- data[i, ]

    # Determine orientation
    orientation <- ifelse(gene$xmax > gene$xmin, 1, -1)

    # Calculate x coordinate of flange
    flangex <- (-orientation * arrowhead_width) + gene$xmax

    # Set arrow and arrowhead heights; it's convenient to divide these by two
    # for calculating y coordinates on the polygon
    arrowhead_height <- as.numeric(grid::convertHeight(x$arrowhead_height, "native")) / 2

    # Arrowhead defaults to 4 mm, unless the gene is shorter in which case the
    # gene is 100% arrowhead
    params$arrowhead_width <- as.numeric(
      grid::convertWidth(params$arrowhead_width, "native")
    )

    # Create polygon grob
    pg <- grid::polygonGrob(
      x = c(
        gene$xmin,
        gene$xmin,
        flangex,
        flangex,
        gene$xmax,
        flangex,
        flangex
      ),
      y = c(
        gene$y + arrow_body_height,
        gene$y - arrow_body_height,
        gene$y - arrow_body_height,
        gene$y - arrowhead_height,
        gene$y,
        gene$y + arrowhead_height,
        gene$y + arrow_body_height
      ),
      gp = grid::gpar(
        fill = ggplot2::alpha(gene$fill, gene$alpha),
        col = ggplot2::alpha(gene$colour, gene$alpha),
        lty = gene$linetype,
        lwd = gene$size * ggplot2::.pt
      )
    )

    # Return the polygon grob
    pg
  })

  class(grobs) <- "gList"
  grid::setChildren(x, grobs)
}
