#' @export gsoap_plot
#' @import ggplot2
#' @importFrom ggforce geom_circle geom_arc_bar
#' @importFrom ggrepel geom_text_repel


geom_radius_legend = function(radius, x, y, no=4, font.size = 10, labeller) {

  # Create sequence of n radius placeholders
  if (length(radius) > no) {
    radius = seq(min(radius), max(radius), length.out = no)
    #radius = unique(sapply(radius, round_digit))
  }

  # Check if proper labeller is provided
  have_labeller = FALSE
  if (!missing(labeller)) {
    if (!inherits(labeller, "function")) {
      stop("Radius labeller must be a function for converting radius")
    }
    have_labeller = TRUE
  }

  # Legend data frame
  dd = data.frame(r=radius,
                  start=0,
                  end=2*pi,
                  x = x,
                  y = y + radius - max(radius),
                  maxr = max(radius))

  # Label placeholders
  if(have_labeller) {
    dd$label = labeller(dd$r)
  } else {
    dd$label = as.character(dd$r)
  }

  # Create legend circles
  legend_circles = geom_arc_bar(aes_(x0=~x, y0=~y, r0=~r, r=~r, start=~start, end=~end),
                                data=dd,
                                inherit.aes=FALSE)

  # Create legend segments
  legend_segments = geom_segment(aes_(x=~x, xend=~x+maxr*1.5, y=~y+r, yend=~y+r),
                                 data=dd,
                                 inherit.aes=FALSE)

  # Create legend text/labels
  legend_text = geom_text(aes_(x=~x+maxr*1.6, y=~y+r, label=~label),
                          data=dd,
                          hjust='left',
                          size = font.size * 0.352777778,
                          inherit.aes=FALSE)

  # Return legend as list
  list(legend_circles, legend_segments, legend_text)
}

#' A function to create Gene Set Overrepresentation Analysis Plot
#'
#' Some description
#'
#' More description
#'
#' @param layout a data frame, a layout extracted from the gsoap object obtained from the \code{\link{create_gsoap_layout}}.
#' @param as.color a character or integer, indicating name or index of the column with color aesthetic.
#' @param as.alpha a character or integer, indicating name or index of the column with alpha aesthetic.
#' @param which.labels a character or integer vector, indicating gene sets (rows)
#' whose names (rownames) should be show in the plot
#' @param show.size.guide a boolean indicating whether to show size guide.
#' @param show.color.guide a boolean indicating whether to show color guide.
#' @param size.guide.loc a numeric vector of length 2 with the size guide x and y coordinates.
#' @param size.guide.no an integer indicating number of circles of the size guide.
#' @param default.color a character indicating color to be used if no color aesthetics is used.
#' @param default.alpha a value between 0 and 1 indicating alpha to be used if no alpha aesthetics is used.
#' @param range.alpha a numeric vector of length 2 with lower and upper limit of alpha range.
#' @param viridis.option a character vector indicating Viridis color scale to be used,
#' options include \emph{viridis} (default), \emph{magma}, \emph{plasma}, \emph{inferno}, \emph{cividis}.
#' @param viridis.direction a direction of the Viridis color palette.
#' @param viridis.range a numeric vector of length 2 with lower and upper limit of the Viridis color palette.
#' @param label.alpha  a value between 0 and 1 indicating alpha of the labels.
#' @param segment.alpha  a value between 0 and 1 indicating alpha of the labels repels.
#' @param repel.xynudges a numeric vector of length 2 with horizontal and vertical adjustments to nudge the position of labels.
#' @param repel.hvjust a numeric vector of length 2 with horizontal and vertical justification of labels.
#' @param repel.direction a character indicating direction in which to adjust position of labels, options are \emph{both} (default), \emph{x} and \emph{y}.
#' @param base.fontsize an integer indicating base font size.
#' @param size.guide.fontsize an integer indicating size guide font size.
#' @param label.fontsize an integer indicating labels font size.
#' @param title a text for the title.
#' @param subtitle a text for the subtitle.
#' @param xlabel a text for the x-axis label.
#' @param ylabel a text for the y-axis label.
#' @param xlimits a numeric vector of length 2 indicating lower and upper x-axis limits
#' @param ylimits a numeric vector of length 2 indicating lower and upper y-axis limits
#' @param void a boolean indicating whether theme void should be used.
#'
#' @return A gsoap plot, where individual gene sets are represented by circles,
#' whose size reflects the number of query gene members of the given set,
#' and whose mutual proximities reflect the number common gene members.
#'
#' @author Tomas Tokar <tomastokar@gmail.com>
#'
#' @examples
#' data(pxgenes)
#'
#' l = gsoap_layout(pxgenes, 'Members', 'p.value')
#' gsoap_plot(l$layout, as.color = 'Cluster', as.alpha = 'Centrality')
#'
gsoap_plot = function(layout,
                      as.color = NULL,
                      as.alpha = NULL,
                      which.labels = NULL,
                      show.size.guide = TRUE,
                      show.color.guide = TRUE,
                      show.alpha.guide = TRUE,
                      size.guide.loc = c(1., 0),
                      size.guide.no = 4,
                      default.color = 'black',
                      default.alpha = 0.5,
                      range.alpha = c(0.1, 0.9),
                      viridis.option = 'viridis',
                      viridis.direction = 1,
                      viridis.range = c(0, 1),
                      label.alpha = 0.6,
                      segment.alpha = 0.5,
                      repel.xynudges = c(0.1, 0.1),
                      repel.hvjust = c(0, 0),
                      repel.direction = 'both',
                      base.fontsize = 8,
                      size.guide.fontsize = base.fontsize - 2,
                      label.fontsize = base.fontsize - 1,
                      title = NULL,
                      subtitle = NULL,
                      xlabel = 'proj. 1',
                      ylabel = 'proj. 2',
                      xlimits = c(-0.1, 1.1),
                      ylimits = c(-0.1, 1.1),
                      void = FALSE){ # font size of repel
  # -------------
  # Check inputs
  #--------------
  if (missing(layout)){
    stop('Layout is missing')
  }
  if (!is.data.frame(layout)){
    stop('Layout is not data frame')
  }
  if (!is.character(rownames(layout))){
    stop('Input has missing or improper rownames')
  }
  if (any(is.na(layout))){
    warning("Layout contains missing values")
  }

  # ----------------------------
  # Set default color and alpha
  # ----------------------------
  # Set defaults
  update_geom_defaults(ggforce::GeomCircle,
                       list(colours = default.color,
                            fill = default.color,
                            alpha = default.alpha))
  # --------------------
  # Initialize the plot
  # --------------------
  p = ggplot(data = layout)
  # Add circles
  p = p + geom_circle(aes(x0 = x, y0 = y, r = radius))
  # Set fixed coords
  p = p + coord_fixed()
  # Set theme
  if (void){
    p = p + theme_void(base_size = base.fontsize)
  } else {
    p = p + theme_light(base_size = base.fontsize)
  }
  # Set x and y labels
  p = p + labs(title = title,
               subtitle = subtitle,
               xlab = xlabel,
               ylab = ylabel)
  # --------------
  # Resolve sizer
  # --------------
  # Set sizer
  if(show.size.guide){
    size.labeller = function(r){
      s = max(layout[,'size']) / max(layout[,'radius'])^2
      l = as.character(round(r^2 * s))
      return(l)
    }
    p = p + geom_radius_legend(layout$radius,
                               size.guide.loc[1],
                               size.guide.loc[2],
                               no = size.guide.no,
                               font.size = size.guide.fontsize,
                               labeller = size.labeller)
  }
  # --------------
  # Resolve color
  # --------------
  # Set color if color is used
  if(!is.null(as.color)){
    if (is.numeric(layout[,as.color])){
      p = p + aes_string(color = as.color, fill = as.color)
      p = p + scale_color_viridis_c(begin = viridis.range[1],
                                    end = viridis.range[2],
                                    direction = viridis.direction,
                                    option = viridis.option)
      p = p + scale_fill_viridis_c(begin = viridis.range[1],
                                   end = viridis.range[2],
                                   direction = viridis.direction,
                                   option = viridis.option)
    }
    if (is.factor(layout[,as.color]) | is.character(layout[,as.color])){
      p = p + aes_string(color = as.color, fill = as.color)
      p = p + scale_color_viridis_d(begin = viridis.range[1],
                                    end = viridis.range[2],
                                    direction = viridis.direction,
                                    option = viridis.option)
      p = p + scale_fill_viridis_d(begin = viridis.range[1],
                                   end = viridis.range[2],
                                   direction = viridis.direction,
                                   option = viridis.option)
    }
  }
  # --------------
  # Resolve alpha
  # --------------
  # Set alpha if alpha is used
  if(!is.null(as.alpha)){
    p = p + aes_string(alpha = as.alpha) # !!!!
    p = p + scale_alpha(range = range.alpha,
                        guide = guide_legend(override.aes = list(fill=default.color)))
  }
  # ---------------
  # Resolve guides
  # ---------------
  if (!show.color.guide){
    p = p + guides(color = FALSE, fill = FALSE)
  }
  if (!show.alpha.guide){
    p = p + guides(alpha = FALSE)
  }
  # ---------------
  # Add annotation
  # ---------------
  if(!is.null(which.labels)){
    layout.subset = layout[which.labels,]
    p = p + geom_text_repel(data = layout[which.labels,],
                            nudge_x = repel.xynudges[1],
                            nudge_y = repel.xynudges[2],
                            hjust = repel.hvjust[1],
                            vjust = repel.hvjust[2],
                            force = 10,
                            direction = repel.direction,
                            segment.alpha = segment.alpha,
                            alpha = label.alpha,
                            color = 'black',
                            size = label.fontsize * 0.352777778,
                            max.iter = 1e+4,
                            aes(x = x,
                                y = y,
                                label = rownames(layout[which.labels,])))
  }
  # -------------
  # Resolve rest
  # -------------
  # Return the plot
  return(p)
}
