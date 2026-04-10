library(shiny)
library(ggplot2)
library(prosail)

ui <- fluidPage(
  titlePanel("PROSAIL BRDF Explorer"),
  sidebarLayout(
    sidebarPanel(
      h4("Canopy parameters"),
      numericInput("N",      "N (leaf structure)",       1.5,   min = 1.0,  max = 3.5,  step = 0.1),
      numericInput("cab",    "cab (chlorophyll, ug/cm2)", 40,   min = 0,    max = 80,   step = 1),
      numericInput("car",    "car (carotenoids, ug/cm2)", 8,    min = 0,    max = 30,   step = 1),
      numericInput("cbrown", "cbrown (brown pigment)",   0.0,   min = 0,    max = 1,    step = 0.01),
      numericInput("cw",     "cw (water thickness, cm)", 0.015, min = 0,    max = 0.05, step = 0.001),
      numericInput("cm",     "cm (dry matter, g/cm2)",   0.009, min = 0,    max = 0.03, step = 0.001),
      numericInput("lai",    "lai (leaf area index)",    3.0,   min = 0.1,  max = 8.0,  step = 0.1),
      numericInput("lidfa",  "lidfa (leaf angle, deg)",  50,    min = 0,    max = 90,   step = 1),
      numericInput("hspot",  "hspot (hotspot)",          0.1,   min = 0,    max = 1,    step = 0.01),
      h4("Soil"),
      numericInput("rsoil", "rsoil (brightness)",        1.0,   min = 0.0,  max = 2.0,  step = 0.05),
      numericInput("psoil", "psoil (dry fraction 1=dry)", 1.0,  min = 0.0,  max = 1.0,  step = 0.05),
      h4("Solar illumination"),
      checkboxInput("manual_solar", "Manual override (sza/saa)", FALSE),
      conditionalPanel(
        condition = "input.manual_solar == false",
        dateInput("date_utc", "Date (UTC)", value = Sys.Date()),
        fluidRow(
          column(6, numericInput("hour_utc",   "Hour (UTC)",   12, min = 0, max = 23, step = 1)),
          column(6, numericInput("minute_utc", "Minute (UTC)",  0, min = 0, max = 59, step = 1))
        ),
        numericInput("lat", "Latitude",  0.0, min = -90,  max = 90,  step = 0.1),
        numericInput("lon", "Longitude", 0.0, min = -180, max = 180, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.manual_solar == true",
        numericInput("sza", "sza (deg)", 30, min = 0, max = 90,  step = 1),
        numericInput("saa", "saa (deg)",  0, min = 0, max = 360, step = 1)
      ),
      h4("Geometry grid"),
      numericInput("vza_max",  "vza_max (deg)", 70, min = 0, max = 85, step = 1),
      numericInput("step_deg", "step_deg",       10, min = 1, max = 20, step = 1),
      h4("Vegetation index"),
      selectInput("vi", "Index", c("NDVI", "EVI", "SAVI", "NDRE")),
      actionButton("run", "Run simulation")
    ),
    mainPanel(
      plotOutput("brdf_plot", height = "650px"),
      verbatimTextOutput("solar_text")
    )
  )
)

server <- function(input, output, session) {

  # ---------- Solar geometry (NOAA algorithm, no external package) ------------
  solar_angles <- reactive({
    if (isTRUE(input$manual_solar)) {
      list(sza = input$sza, saa = input$saa)
    } else {
      d         <- as.Date(input$date_utc)
      doy       <- as.integer(strftime(d, format = "%j"))
      local_time <- input$hour_utc + input$minute_utc / 60
      gamma <- 2 * pi / 365 * (doy - 1 + (local_time - 12) / 24)
      eqtime <- 229.18 * (0.000075 +
                            0.001868 * cos(gamma)  - 0.032077 * sin(gamma) -
                            0.014615 * cos(2*gamma) - 0.040849 * sin(2*gamma))
      decl <- 0.006918  - 0.399912 * cos(gamma)  + 0.070257 * sin(gamma) -
        0.006758 * cos(2*gamma) + 0.000907 * sin(2*gamma) -
        0.002697 * cos(3*gamma) + 0.00148  * sin(3*gamma)
      tst     <- local_time * 60 + eqtime + 4 * input$lon
      ha      <- (tst / 4) - 180
      lat_rad <- input$lat * pi / 180
      ha_rad  <- ha        * pi / 180
      cos_zen <- max(-1, min(1,
                             sin(lat_rad) * sin(decl) +
                               cos(lat_rad) * cos(decl) * cos(ha_rad)))
      zen_rad <- acos(cos_zen)
      sza     <- max(0, min(90, zen_rad * 180 / pi))
      sin_az  <- -sin(ha_rad) * cos(decl) / sin(zen_rad + 1e-12)
      cos_az  <- (sin(decl) - sin(lat_rad) * cos_zen) /
        (cos(lat_rad) * sin(zen_rad) + 1e-12)
      saa <- (atan2(sin_az, cos_az) * 180 / pi + 180) %% 360
      list(sza = sza, saa = saa)
    }
  })

  # ---------- Main simulation --------------------------------------------------
  run_sim <- eventReactive(input$run, {
    angles <- solar_angles()
    vza    <- seq(0, input$vza_max, by = input$step_deg)
    raa    <- seq(0, 360,           by = input$step_deg)
    grid   <- expand.grid(vza = vza, raa = raa)

    # PRO4SAIL output: 2101 bands, 400-2500 nm, index = wavelength - 400 + 1
    compute_vi <- function(R) {
      if (is.null(R) || length(R) < 401) return(NA_real_)
      blue    <- R[71]   # 470 nm
      red     <- R[271]  # 670 nm
      rededge <- R[306]  # 705 nm
      nir     <- R[401]  # 800 nm
      switch(input$vi,
        "NDVI" = (nir - red)  / (nir + red),
        "EVI"  = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1),
        "SAVI" = 1.5 * (nir - red) / (nir + red + 0.5),
        "NDRE" = (nir - rededge) / (nir + rededge),
        NA_real_
      )
    }

    # Soil spectrum: psoil=1 -> dry/bright, psoil=0 -> wet/dark
    rsoil_vec <- input$rsoil * (
      input$psoil        * prosail::SpecSOIL$Dry_Soil +
      (1 - input$psoil)  * prosail::SpecSOIL$Wet_Soil
    )

    # PROSPECT-PRO leaf parameter data.frame required by PRO4SAIL
    Input_PROSPECT <- data.frame(
      CHL = input$cab,    # chlorophyll a+b (ug/cm2)
      CAR = input$car,    # carotenoids     (ug/cm2)
      ANT = input$cbrown, # brown pigments  (ug/cm2)
      EWT = input$cw,     # equivalent water thickness (cm)
      LMA = input$cm,     # leaf mass per area (g/cm2)
      N   = input$N       # leaf structure parameter
    )

    values <- vapply(seq_len(nrow(grid)), function(i) {
      res <- tryCatch(
        prosail::PRO4SAIL(
          Input_PROSPECT = Input_PROSPECT,
          TypeLidf       = 2,           # ellipsoidal LIDF
          LIDFa          = input$lidfa,
          lai            = input$lai,
          q              = input$hspot,
          tts            = angles$sza,
          tto            = grid$vza[i],
          psi            = grid$raa[i],
          rsoil          = rsoil_vec
        ),
        error = function(e) NULL
      )
      compute_vi(if (!is.null(res)) as.numeric(res$rso) else NULL)
    }, numeric(1))

    grid$vi  <- values
    grid$sza <- angles$sza
    grid$saa <- angles$saa
    grid
  })

  # ---------- Polar BRDF plot --------------------------------------------------
  output$brdf_plot <- renderPlot({
    df <- run_sim()
    validate(need(nrow(df) > 0, "No simulation results yet. Press Run simulation."))
    sza_val <- unique(df$sza)[1]
    saa_val <- unique(df$saa)[1]
    ggplot(df, aes(x = raa, y = vza, fill = vi)) +
      geom_tile() +
      coord_polar(theta = "x") +
      annotate("point", x = saa_val %% 360, y = sza_val,
               colour = "orange", size = 4) +
      scale_x_continuous(
        limits = c(0, 360),
        breaks = seq(0, 315, by = 45),
        labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
      ) +
      scale_y_continuous(limits = c(0, max(df$vza, na.rm = TRUE))) +
      scale_fill_viridis_c(option = "magma", na.value = "grey80") +
      labs(
        title = paste("BRDF polar plot -", input$vi),
        x     = "Relative azimuth angle (deg)",
        y     = "View zenith angle (deg)",
        fill  = input$vi
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  # ---------- Solar angle display ----------------------------------------------
  output$solar_text <- renderText({
    angles <- solar_angles()
    sprintf("SZA: %.2f deg\nSAA: %.2f deg", angles$sza, angles$saa)
  })
}

shinyApp(ui, server)
