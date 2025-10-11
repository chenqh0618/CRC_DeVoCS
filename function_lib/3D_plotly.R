plot_wireframe_sphere <- function(sphere_list = list(sphere1 = c(0, 0, 0, 1, 15),
                                                     sphere2 = c(1,2,2, 3, 15)) ,
                                  line_color = 'rgb(100, 100, 100)',
                                  line_width = 1.5,
                                  points_df = NULL,
                                  x_col = "x", y_col = "y", z_col = "z",
                                  title = NULL) {
  
  require(plotly)
  fig <- plot_ly()
  
  for (i in 1:length(sphere_list)) {
    center <- sphere_list[[i]][1:3]
    radius <- sphere_list[[i]][4]
    grid_density_deg <- sphere_list[[i]][5]
    
    if (length(center) != 3 || !is.numeric(center)) stop("center should be a vector with 3 values")
    if (!is.numeric(radius) || radius <= 0) stop("radius must > 0")
    if (!is.numeric(grid_density_deg) || grid_density_deg <= 0) stop("grid_density_deg must > 0")
    if (grid_density_deg > 180) warning("when grid_density_deg > 180，the grid will be too sparse")
    
    n_phi <- max(2, ceiling(180 / grid_density_deg) + 1)
    n_theta <- max(2, ceiling(360 / grid_density_deg) + 1)
    
    phi_vals <- seq(0, pi, length.out = n_phi)
    theta_vals <- seq(0, 2 * pi, length.out = n_theta)
    
    x_rel <- radius * outer(sin(phi_vals), cos(theta_vals))
    y_rel <- radius * outer(sin(phi_vals), sin(theta_vals))
    z_rel <- radius * outer(cos(phi_vals), rep(1, n_theta))
    
    x_coords_sphere <- center[1] + x_rel
    y_coords_sphere <- center[2] + y_rel
    z_coords_sphere <- center[3] + z_rel
    
    for (i in 1:n_phi) {
      fig <- add_trace(fig, x = x_coords_sphere[i, ], y = y_coords_sphere[i, ], z = z_coords_sphere[i, ],
                       type = 'scatter3d', mode = 'lines',
                       line = list(color = line_color, width = line_width),
                       showlegend = FALSE, hoverinfo = 'none') 
    }
    
    for (j in 1:n_theta) {
      fig <- add_trace(fig, x = x_coords_sphere[, j], y = y_coords_sphere[, j], z = z_coords_sphere[, j],
                       type = 'scatter3d', mode = 'lines',
                       line = list(color = line_color, width = line_width),
                       showlegend = FALSE, hoverinfo = 'none')
    }
    
  }
  

  for (n in 1:length(points_df)) {
    min_coords_data <- c(Inf, Inf, Inf)
    max_coords_data <- c(-Inf, -Inf, -Inf)
    has_data_points = FALSE
    
    data_points_df <- points_df[[n]]$df
    size_col <- points_df[[n]]$size_col
    color_col <- points_df[[n]]$color_col
    colorscale <- points_df[[n]]$colorscale
    point_default_color <- points_df[[n]]$point_default_color
    marker_opacity <- points_df[[n]]$marker_opacity
    marker_size_min <- points_df[[n]]$marker_size_min
    marker_size_max <- points_df[[n]]$marker_size_max
    show_colorbar <- points_df[[n]]$show_colorbar
    show_legend_points <- points_df[[n]]$show_legend_points
    point_trace_name <- points_df[[n]]$point_trace_name
    
    
    if (!is.null(data_points_df) && inherits(data_points_df, "data.frame") && nrow(data_points_df) > 0) {
      has_data_points = TRUE
      coord_cols <- c(x_col, y_col, z_col)
      if (!all(coord_cols %in% names(data_points_df))) {
        stop("the x,y,z is not found in plot_data")
      }
      x_pts <- data_points_df[[x_col]]
      y_pts <- data_points_df[[y_col]]
      z_pts <- data_points_df[[z_col]]
      
      if (!all(sapply(list(x_pts, y_pts, z_pts), is.numeric))) {
        stop("x, y, z must be numric")
      }
      
      min_coords_data <- c(min(x_pts, na.rm = TRUE), min(y_pts, na.rm = TRUE), min(z_pts, na.rm = TRUE))
      max_coords_data <- c(max(x_pts, na.rm = TRUE), max(y_pts, na.rm = TRUE), max(z_pts, na.rm = TRUE))
      
      marker_list <- list(opacity = marker_opacity)
      show_points_legend_flag = FALSE 
      
      color_values_formula = NULL
      if (!is.null(color_col)) {
        if (!color_col %in% names(data_points_df)) {
          warning("the color_col is not found in plot_data, a default color will be use")
          marker_list$color <- point_default_color
        } else {
          color_vec <- data_points_df[[color_col]]
          color_values_formula = as.formula(paste0("~`", color_col, "`"))  
          
          if (is.numeric(color_vec)) {
            marker_list$colorscale <- colorscale
            marker_list$showscale <- show_colorbar
            marker_list$colorbar <- list(title = gsub("_", " ", color_col), 
                                         len=0.6, y=0.5,
                                         titlefont=list(size=10), tickfont=list(size=8))
          } else {
            color_values_formula = as.formula(paste0("~as.factor(`", color_col, "`)"))
            if(show_legend_points) {
              show_points_legend_flag = TRUE 
            }
          }
          marker_list$color <- color_values_formula
        }
      } else {
        marker_list$color <- point_default_color
      }
      
      size_values_formula = NULL
      if (!is.null(size_col)) {
        if (!size_col %in% names(data_points_df)) {
          warning("the size_col is not found in plot_data, a default value will be use")
          marker_list$size <- (marker_size_min + marker_size_max) / 2
        } else {
          size_vec <- data_points_df[[size_col]]
          size_vec[is.na(size_vec)] = min(size_vec, na.rm = T)
          if (!is.numeric(size_vec)) {
            warning("the size_col is not numric, a deafult value will be use")
            marker_list$size <- (marker_size_min + marker_size_max) / 2
          } else {
            fact_dist = max(size_vec) - min(size_vec)
            set_dist = marker_size_max - marker_size_min
            size_vec <-  (size_vec-min(size_vec))/fact_dist*set_dist + marker_size_min
            marker_list$size <- size_vec
            marker_list$sizemin <- marker_size_min
            marker_list$sizemax <- marker_size_max
            marker_list$sizemode <- 'diameter' 
          }
        }
      } else {
        marker_list$size <- (marker_size_min + marker_size_max) / 2
      }
      
      hover_text_list <- list(
        paste("X:", round(x_pts, 2)),
        paste("Y:", round(y_pts, 2)),
        paste("Z:", round(z_pts, 2))
      )
      if (!is.null(color_col) && color_col %in% names(data_points_df)) {
        col_val <- data_points_df[[color_col]]
        hover_text_list[[length(hover_text_list) + 1]] <- paste(
          gsub("_", " ", color_col), ":",
          if(is.numeric(col_val)) round(col_val, 2) else col_val
        )
      }

      if (!is.null(size_col) && size_col %in% names(data_points_df) && is.numeric(data_points_df[[size_col]])) {
        hover_text_list[[length(hover_text_list) + 1]] <- paste(
          gsub("_", " ", size_col), ":", round(data_points_df[[size_col]], 2)
        )
      }
      hover_text_combined <- apply(do.call(cbind, hover_text_list), 1, paste, collapse = "<br>")
      
      
      fig <- add_trace(fig,
                       data = data_points_df,
                       x = as.formula(paste0("~`", x_col, "`")), 
                       y = as.formula(paste0("~`", y_col, "`")),
                       z = as.formula(paste0("~`", z_col, "`")),
                       type = 'scatter3d',
                       mode = 'markers',
                       marker = marker_list,
                       text = hover_text_combined, 
                       hoverinfo = 'text',         
                       name = point_trace_name,    
                       showlegend = show_points_legend_flag
      )
      
    } 
    }  }

  
  all_x <- c(as.vector(x_coords_sphere), if(has_data_points) x_pts else NULL)
  all_y <- c(as.vector(y_coords_sphere), if(has_data_points) y_pts else NULL)
  all_z <- c(as.vector(z_coords_sphere), if(has_data_points) z_pts else NULL)
  

  buffer_factor = 0.05
  x_range <- range(all_x, na.rm = TRUE) + c(-1, 1) * diff(range(all_x, na.rm = TRUE)) * buffer_factor
  y_range <- range(all_y, na.rm = TRUE) + c(-1, 1) * diff(range(all_y, na.rm = TRUE)) * buffer_factor
  z_range <- range(all_z, na.rm = TRUE) + c(-1, 1) * diff(range(all_z, na.rm = TRUE)) * buffer_factor
  if(diff(x_range) == 0) x_range <- x_range + c(-max(1, abs(x_range[1]*0.1)), max(1, abs(x_range[1]*0.1)))
  if(diff(y_range) == 0) y_range <- y_range + c(-max(1, abs(y_range[1]*0.1)), max(1, abs(y_range[1]*0.1)))
  if(diff(z_range) == 0) z_range <- z_range + c(-max(1, abs(z_range[1]*0.1)), max(1, abs(z_range[1]*0.1)))
  
  
  final_title <- if(has_data_points) title else gsub(" with Data Points", "", title)  
  
  fig <- fig %>% layout(
    # title = final_title,
    # scene = list(
    #   xaxis = list(title = "X", range = x_range, showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, ticks=""),
    #   yaxis = list(title = "Y", range = y_range, showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, ticks=""),
    #   zaxis = list(title = "Z", range = z_range, showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, ticks=""),
    #   aspectmode = 'data' # 保持数据纵横比
    # ),
    scene =  list(
      xaxis = list(showticklabels = T, 
                   zeroline = T, showgrid = T, 
                   showline = T, showaxis = T, 
                   showlable =T,  linewidth =4,
                   titlefont = list(size = 42), tickfont = 36),
      yaxis = list(showticklabels = T, 
                   zeroline = T, showgrid = T, 
                   showline = T, showaxis = T, 
                   showlable =T, linewidth =4,
                   titlefont = list(size = 42), tickfont = 36),
      zaxis = list(showticklabels = T, 
                   zeroline = T, showgrid = T, 
                   showline = T, showaxis = T, 
                   showlable =T, linewidth =4,
                   titlefont = list(size = 42), tickfont = 36))
  )
  
  return(fig)
}
