library(shiny)
library(ggplot2)
library(suncalc)
library(prosail)

ui <- fluidPage(
  titlePanel("PROSAIL BRDF Explorer"),
  sidebarLayout(
    sidebarPanel(
      h4("Canopy parameters"),
      numericInput("N", "N (leaf structure)", 1.5, min = 1.0, max = 3.5, step = 0.1),
      numericInput("cab", "cab (chlorophyll)", 40, min = 0, max = 80, step = 1),
      numericInput("car", "car (carotenoids)", 8, min = 0, max = 30, step = 1),
      numericInput("cbrown", "cbrown (brown pigment)", 0.0, min = 0, max = 1, step = 0.01),
      numericInput("cw", "cw (water thickness)", 0.015, min = 0, max = 0.05, step = 0.001),
      numericInput("cm", "cm (dry matter)", 0.009, min = 0, max = 0.03, step = 0.001),
      numericInput("lai", "lai (leaf area index)", 3.0, min = 0.1, max = 8.0, step = 0.1),
      numericInput("lidfa", "lidfa (leaf angle)", 50, min = 0, max = 90, step = 1),
      numericInput("hspot", "hspot (hotspot)", 0.1, min = 0, max = 1, step = 0.01),
      h4("Soil"),
      numericInput("rsoil", "rsoil (brightness)", 1.0, min = 0.0, max = 2.0, step = 0.05),
      numericInput("psoil", "psoil (moisture mix)", 0.5, min = 0.0, max = 1.0, step = 0.05),
      h4("Solar illumination"),
      checkboxInput("manual_solar", "Manual override (sza/saa)", FALSE),
      conditionalPanel(
        condition = "input.manual_solar == false",
        dateInput("date_utc", "Date (UTC)", value = Sys.Date()),
        fluidRow(
          column(6, numericInput("hour_utc", "Hour (UTC)", 12, min = 0, max = 23, step = 1)),
          column(6, numericInput("minute_utc", "Minute (UTC)", 0, min = 0, max = 59, step = 1))
        ),
        numericInput("lat", "Latitude", 0.0, min = -90, max = 90, step = 0.1),
        numericInput("lon", "Longitude", 0.0, min = -180, max = 180, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.manual_solar == true",
        numericInput("sza", "sza (deg)", 30, min = 0, max = 90, step = 1),
        numericInput("saa", "saa (deg)", 0, min = 0, max = 360, step = 1)
      ),
      h4("Geometry grid"),
      numericInput("vza_max", "vza_max (deg)", 70, min = 0, max = 85, step = 1),
      numericInput("step_deg", "step_deg", 10, min = 1, max = 20, step = 1),
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
  solar_angles <- reactive({
    if (isTRUE(input$manual_solar)) {
      list(sza = input$sza, saa = input$saa)
    } else {
      dt <- as.POSIXct(paste(input$date_utc, sprintf("%02d:%02d:00", input$hour_utc, input$minute_utc)), tz = "UTC")
      pos <- suncalc::getSunlightPosition(dt = dt, lat = input$lat, lon = input$lon)
      altitude_deg <- pos$altitude * 180 / pi
      azimuth_deg <- pos$azimuth * 180 / pi
      sza <- max(0, min(90, 90 - altitude_deg))
      saa <- (azimuth_deg + 180) %% 360
      list(sza = sza, saa = saa)
    }
  })

  run_sim <- eventReactive(input$run, {
    angles <- solar_angles()
    vza <- seq(0, input$vza_max, by = input$step_deg)
    raa <- seq(0, 360, by = input$step_deg)
    grid <- expand.grid(vza = vza, raa = raa)

    band_index <- function(wl) wl - 400 + 1
    idx_blue <- band_index(470)
    idx_red <- band_index(670)
    idx_rededge <- band_index(705)
    idx_nir <- band_index(800)

    compute_vi <- function(refl) {
      blue <- refl[idx_blue]
      red <- refl[idx_red]
      rededge <- refl[idx_rededge]
      nir <- refl[idx_nir]
      if (input$vi == "NDVI") {
        (nir - red) / (nir + red)
      } else if (input$vi == "EVI") {
        2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
      } else if (input$vi == "SAVI") {
        1.5 * (nir - red) / (nir + red + 0.5)
      } else {
        (nir - rededge) / (nir + rededge)
      }
    }

    values <- vapply(seq_len(nrow(grid)), function(i) {
      vza_i <- grid$vza[i]
      raa_i <- grid$raa[i]
      refl <- PROSAIL::run_prosail(
        N = input$N,
        Cab = input$cab,
        Car = input$car,
        Cbrown = input$cbrown,
        Cw = input$cw,
        Cm = input$cm,
        LAI = input$lai,
        lidfa = input$lidfa,
        hspot = input$hspot,
        rsoil = input$rsoil,
        psoil = input$psoil,
        tts = angles$sza,
        tto = vza_i,
        psi = raa_i
      )
      compute_vi(refl)
    }, numeric(1))

    grid$vi <- values
    grid$sza <- angles$sza
    grid$saa <- angles$saa
    grid
  })

  output$brdf_plot <- renderPlot({
    df <- run_sim()
    ggplot(df, aes(x = raa, y = vza, fill = vi)) +
      geom_tile() +
      coord_polar(theta = "x") +
      geom_point(aes(x = 0, y = unique(sza)), color = "orange", size = 3) +
      scale_x_continuous(limits = c(0, 360)) +
      scale_y_continuous(limits = c(0, max(df$vza))) +
      labs(x = "Relative azimuth angle (deg)", y = "View zenith angle (deg)", fill = input$vi) +
      theme_minimal()
  })

  output$solar_text <- renderText({
    angles <- solar_angles()
    sprintf("SZA: %.2f deg\nSAA: %.2f deg", angles$sza, angles$saa)
  })
}

shinyApp(ui, server)
