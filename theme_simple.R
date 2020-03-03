theme_simple <- function(axis_lines = TRUE,   # do we want axis lines or not
                         base_size = 12,      # base size of all text
                         base_family= 'Avenir', # font family 
                         line_size = 0.2) {   # line thickness
        
        # lets start with the minimal theme because it is really simple.
        th <- ggplot2::theme_minimal(base_family = base_family, 
                                     base_size = base_size)
        
        # remove the grid lines 
        th <- th + theme(panel.grid=element_blank())
        
        # if we want axis lines
        if (axis_lines) {
                # We add axis lines and give them our preferred thickness
                th <- th + theme(axis.line = element_line(size =  line_size))
                th <- th + theme(axis.ticks = element_line(size = line_size))
        }
        
        # Let's give the axis text a bit more space
        th <- th + theme(axis.text.x=element_text(margin=margin(t=5)))
        th <- th + theme(axis.text.y=element_text(margin=margin(r=5)))
        
        # And also a bit more space for the axis titles
        th <- th + theme(axis.title.x = element_text(margin=margin(t=10)))
        th <- th + theme(axis.title.y = element_text(margin=margin(r=10)))
        
        return (th)
}