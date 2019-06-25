theme_clean <- function(grid=TRUE, legend.position=NA, base_size = 12,
                        base_family='Helvetica', highlight_family='Helvetica') {
        
        # credit to Henrik Lindberg https://github.com/halhen
        th <- ggplot2::theme_minimal(base_family = base_family, base_size = base_size)
        
        th <- th + theme(text = element_text(color='#333333'))
        th <- th + theme(legend.background = element_blank())
        th <- th + theme(legend.key = element_blank())
 
        
        th <- th + theme(panel.grid=element_blank())
        
        th <- th + theme(axis.text = element_text(family=highlight_family))
        th <- th + theme(axis.ticks = element_blank())
        
        th <- th + theme(axis.text.x=element_text(margin=margin(t=5)))
        th <- th + theme(axis.text.y=element_text(margin=margin(r=5)))
        th <- th + theme(axis.line = element_line(color = '#333333',  size = 0.5))
        if (!is.na(legend.position)) th <- th + theme(legend.position = legend.position)
        
        return (th)
        
}
