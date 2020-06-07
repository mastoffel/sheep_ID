# theme for most plots

theme_simple <- function(grid_lines = TRUE, axis_lines = FALSE, base_size = 12,
                        base_family='Lato', line_size = 0.2) {
        
        # base 
        th <- ggplot2::theme_minimal(base_family = base_family, 
                                     base_size = base_size)
        
        # start with blank
        th <- th + theme(panel.grid=element_blank())
        
        if (grid_lines) {
                th <- th + theme(panel.grid.minor = element_blank())  
                th <- th + theme(panel.grid.major = element_line(color = "#C0C0C0",
                                                                 size =  line_size))
                th <- th + theme(axis.ticks = element_blank())
                
        } 
        
        if (axis_lines) {
                th <- th + theme(axis.line = element_line(size =  line_size))
                th <- th + theme(axis.ticks = element_line(size = line_size))
        }
       
        th <- th + theme(axis.text.x=element_text(margin=margin(t=5)))
        th <- th + theme(axis.text.y=element_text(margin=margin(r=5)))
        
        th <- th + theme(axis.title.x = element_text(margin=margin(t=10)))
        th <- th + theme(axis.title.y = element_text(margin=margin(r=10)))
        
        # legend
        th <- th + theme(legend.background = element_blank())
        th <- th + theme(legend.key = element_blank())
        
       # th <- th + theme(axis.line = element_line(color = '#333333',  size = 0.5))
        #th <- th + theme(axis.text = element_text(family=highlight_family))
        # if (!is.na(legend.position)) th <- th + theme(legend.position = legend.position)
        
        return (th)
        
}
