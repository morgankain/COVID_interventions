fig1_data$county <- factor(fig1_data$county)
fig1_data %<>% mutate(name = mapvalues(name, from = c("cases", "deaths", "Reff", "detect"),
                                       to = c("Daily cases", "Daily deaths", "Reff", "detect")))
{
 fig1.1 <- {fig1_data %>% 
   filter(county != "Los Angeles") %>%
     mutate(county = mapvalues(county,
                                   from = c("Santa Clara",
                                            "King",
                                            "Fulton",
                                            "Miami-Dade"),
                                    to  = c(
                                       "Santa Clara County, CA"
                                     , "Seattle, WA"
                                     , "Atlanta, GA"
                                     , "Miami, FL"
                                   ))) %>%
   filter(date < as.Date("2020-06-19")) %>%
   filter(date >= as.Date("2020-02-10")) %>%
  # filter(!(name %in% c("Reff", "Detect"))) %>%
   group_by(county) %>% 
   mutate(nparams = 0.5/length(unique(paramset))) %>% 
   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
              fill  = county, color = county, 
              group = interaction(county, paramset))) +
   geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
   geom_line() +
   geom_point(aes(x = date, y = data), 
              color = "black", size = 0.75) + 
   scale_y_continuous(trans = "sqrt") + 
   scale_fill_manual(guide = F, values = fig1_colors) +
   scale_color_manual(guide = F, values = fig1_colors) +
   facet_grid(name ~ county, scales = "free_y", switch = "y")  +
   ylab("") + 
   theme(
     strip.background = element_blank()
     , strip.placement = "outside"
     , strip.text.x = element_text(size = 11)
     , strip.text.y = element_text(size = 16)
     , axis.text.x = element_text(size = 12)
     , axis.title.y = element_text(margin = margin(l = 1))
     , plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")
     ) + 
     xlab("") #+
 #    labs(tag = "A.") +
  #   theme(plot.tag.position = "topleft"
  #         , plot.tag = element_text(face = "bold" 
  #                                  , margin = margin(b = -16, 
  #                                                    r = -35)))
   }
 
 fig1.2 <- {fig1_data %>% 
   arrange(county, paramset, date) %>%
   filter(county == "Los Angeles") %>%
 #  mutate(county = paste0(county, " County, ", state)) %>%
     mutate(county = mapvalues(county,
                                   from = c("Los Angeles"),
                                    to  = c("Los Angeles, CA")
                                   )) %>%
   filter(date < as.Date("2020-06-19")) %>%
   filter(date >= as.Date("2020-02-10")) %>%
   filter(!(name %in% c("Reff", "Detect"))) %>%
   group_by(county) %>% 
   mutate(nparams = 0.5/length(unique(paramset))) %>% 
   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
              fill  = county, color = county, 
              group = interaction(county,paramset))) +
   geom_ribbon(aes(alpha = I(nparams)), colour = NA) +
   geom_line() +
   geom_point(aes(x = date, y = data), 
              color = "black", size = 0.75) + 
   scale_y_continuous(trans = "sqrt", position = "right") + 
   scale_fill_manual(guide = F, values = fig1_colors[5]) +
   scale_color_manual(guide = F, values = fig1_colors[5]) +
   facet_grid(name ~ county, scales = "free_y", switch = "y")  +
   ylab("") + 
   theme(
     strip.background = element_blank()
     , strip.placement = "outside"
     , strip.text.y = element_blank()
     , strip.text.x = element_text(size = 11)
     , axis.text.x = element_text(size = 12)
     , plot.margin = unit(c(0.25,0.25,0.25,0), "cm")) +
    xlab("")
   }
 
 fig1.3 <- {fig1_data %>% 
   filter(date < as.Date("2020-06-19")) %>%
   filter(date >= as.Date("2020-02-10")) %>%
   filter(name == "Reff") %>%
     mutate(county = mapvalues(county,
                                   from = c("Santa Clara",
                                            "King",
                                            "Fulton",
                                            "Miami-Dade"),
                                    to  = c(
                                       "Santa Clara County, CA"
                                     , "Seattle, WA"
                                     , "Atlanta, GA"
                                     , "Miami, FL"
                                   ))) %>%
# ugly way to calculate 7 day smoothing 
   arrange(county, paramset, date) %>%
   group_by(county, paramset) %>%
   mutate(
     lag3 = lag(mid, 3),
     lag2 = lag(mid, 2),
     lag1 = lag(mid, 1),
     lead1 = lead(mid, 1),
     lead2 = lead(mid, 2),
     lead3 = lead(mid, 3)
   ) %>% 
   rowwise() %>% 
   mutate(mid = mean(c(lag1, lag2, lag3, mid, lead1, lead2, lead3), na.rm = T)) %>% 
   select(-starts_with("lag"), -starts_with("lead")) %>% 
   mutate(paramset = as.character(paramset)) %>%
   group_by(county, state, date) %>% 
   {rbind(., 
          dplyr::summarise(., 
                           med = mean(mid), 
                           lwr = min(mid),
                           upr = max(mid),
                           .groups = "drop") %>% 
            rename(mid = med) %>% 
            mutate(paramset = "summary",
                   name = "Detect",
                   intervention = "Reality",
                   data = NA))} %>% 
   mutate(width = ifelse(paramset == "summary", 1.5, 0.0),
          alpha = ifelse(paramset == "summary", 1, 0.0)) %>%
   ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
              fill  = county, color = county, 
              group = interaction(county, paramset))) +
   geom_ribbon(alpha = 0.25, color = NA) +
   geom_line(aes(size = I(width), alpha = I(alpha))) +
   scale_color_manual(guide = F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
   scale_fill_manual(guide  = F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
   geom_hline(yintercept    = 1, linetype = "dashed") + 
   ylab(expression("R"[e])) + 
   xlab("") +
   scale_y_continuous(position = "left") + 
 #  labs(tag = "B.") +
   ggtitle("Reproduction number") +
   theme(plot.title = element_text(face = "bold")
      , plot.margin = unit(c(0.00,0,0.00,2.01), "cm"))
    # , plot.margin = unit(c(0.00,2.6,0.00,2.01), "cm"))
     # , plot.tag.position = "topleft"
    #  , plot.tag = element_text(face = "bold"
    #    ))
 }
 
 fig1.4 <- {fig1_data %>% 
       filter(date < as.Date("2020-06-19")) %>%
       filter(date >= as.Date("2020-04-15")) %>%
       filter(name == "Reff") %>%
       mutate(county = mapvalues(county,
                                 from = c("Santa Clara",
                                          "King",
                                          "Fulton",
                                          "Miami-Dade"),
                                 to  = c(
                                    "Santa Clara County, CA"
                                    , "Seattle, WA"
                                    , "Atlanta, GA"
                                    , "Miami, FL"
                                 ))) %>%
       # ugly way to calculate 7 day smoothing 
       arrange(county, paramset, date) %>%
       group_by(county, paramset) %>%
       mutate(
          lag3 = lag(mid, 3),
          lag2 = lag(mid, 2),
          lag1 = lag(mid, 1),
          lead1 = lead(mid, 1),
          lead2 = lead(mid, 2),
          lead3 = lead(mid, 3)
       ) %>% 
       rowwise() %>% 
       mutate(mid = mean(c(lag1, lag2, lag3, mid, lead1, lead2, lead3), na.rm = T)) %>% 
       select(-starts_with("lag"), -starts_with("lead")) %>% 
       mutate(paramset = as.character(paramset)) %>%
       group_by(county, state, date) %>% 
       {rbind(., 
              dplyr::summarise(., 
                               med = mean(mid), 
                               lwr = min(mid),
                               upr = max(mid),
                               .groups = "drop") %>% 
                 rename(mid = med) %>% 
                 mutate(paramset = "summary",
                        name = "Detect",
                        intervention = "Reality",
                        data = NA))} %>% 
       mutate(width = ifelse(paramset == "summary", 1.5, 0.0),
              alpha = ifelse(paramset == "summary", 1, 0.0)) %>%
       ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
                  fill  = county, color = county, 
                  group = interaction(county, paramset))) +
       geom_ribbon(alpha = 0.25, color = NA) +
       geom_line(aes(size = I(width), alpha = I(alpha))) +
       scale_color_manual(guide = F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
       scale_fill_manual(guide  = F, values = fig1_colors[c(1, 2, 5, 3, 4)]) +
       geom_hline(yintercept    = 1, linetype = "dashed") + 
       # ylab(expression("R"[e])) + 
       ylab("") +
       xlab("") +
       scale_y_continuous(position = "right") + 
       #  labs(tag = "B.") +
       # ggtitle("Reproduction number") +
       ggtitle(" ") +
       theme(plot.title = element_text(face = "bold")
             , plot.margin = unit(c(0.00,0.75,0,0.6), "cm"))
    # , plot.tag.position = "topleft"
    #  , plot.tag = element_text(face = "bold"
    #    ))
 }
 
# fig1 <- gridExtra::arrangeGrob(fig1.1, fig1.2, fig1.3,  
#                         layout_matrix = matrix(c(1, 1, 2, 3, 3, 3), 
#                                                byrow  = T, nrow = 2),
#                         widths = c(2.6, 1.3, 1.3), heights = c(2.5, 1.5)) 
fig1 <- gridExtra::arrangeGrob(fig1.1, fig1.2, fig1.3,  fig1.4,
                               layout_matrix = matrix(c(1, 1, 2, 3, 4, 4), 
                                                      byrow  = T, nrow = 2),
                               widths = c(2.9, 1.0, 1.3), heights = c(2.5, 1.5)) 

  # gridExtra::grid.arrange(fig1)
  ggsave("figures/Manuscript2/Figure1.pdf",
         fig1,
         device = "pdf",
         width = 12,
         height = 8)
}
