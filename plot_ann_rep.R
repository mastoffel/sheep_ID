p_marginal_effs <- ggplot(inla_preds, aes(froh, prediction)) +
        
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.05, linetype = 2, size = 0.2) +
        geom_line(aes(color = age), size = 0.8) +
        scale_color_viridis_d("Life stage (age)", labels = c("Lamb (0)", "Early life (1,2)", "Mid Life (3,4)", "Late life (5+)")) +
        scale_fill_viridis_d("Life stage (age)", labels = c("Lamb (0)", "Early life (1,2)", "Mid Life (3,4)", "Late life (5+)")) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14,
                     base_family = "Avenir Next") +
        theme(axis.line.y = element_blank(),
             # legend.position = "none",
              axis.text.x = element_text(size = 11)) +
        xlab(expression(Inbreeding~coefficient~F[ROH])) +
        ylab("Predicted\nsurvival probability\n(% per year)")
p_marginal_effs

ggsave(p_marginal_effs, file = "figs/ann_rep.jpg", width = 6.5, height = 4)
ggsave(p_marginal_effs, file = "figs/ann_rep.pdf", device = cairo_pdf, width = 6.5, height = 4)
