# ==============================================================================
# PROSAIL Hemispherical BRDF Analysis - Shiny App
# Correct API based on jbferet/prosail package (jbferet.gitlab.io/prosail)
#
# KEY FACTS about jbferet prosail package:
#   - Exported function: PRO4SAIL()  (NOT prosail() - that does NOT exist)
#   - Input_PROSPECT is a data.frame with cols: CHL, CAR, ANT, EWT, LMA, N
#   - rsoil must be provided (use SpecSOIL$Dry_Soil from the package)
#   - Output is a list; bi-directional reflectance (BRF) is in $rso
#   - Wavelength range: 400-2500 nm, 1 nm steps = 2101 bands
#   - Band index formula: wavelength_nm - 400 + 1  (1-based R indexing)
#     e.g. 470 nm = index 71, 670 nm = 271, 720 nm = 321, 800 nm = 401
# ==============================================================================

# --- Package installation -----------------------------------------------------
required_pkgs <- c("prosail", "suncalc", "shiny", "leaflet", "ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(shiny)
library(leaflet)
library(prosail)   # provides PRO4SAIL() and SpecSOIL
library(suncalc)
library(ggplot2)
library(dplyr)

# ==============================================================================
# Helper functions
# ==============================================================================

# ------------------------------------------------------------------------------
# run_pro4sail(): single PRO4SAIL call with correct argument structure.
# Returns the bi-directional reflectance vector ($rso), or NULL on error.
# ------------------------------------------------------------------------------
run_pro4sail <- function(cab, lai, sza, vza, raa) {
  
  # Input_PROSPECT must be a data.frame (one row = one simulation)
  Input_PROSPECT <- data.frame(
    CHL = cab,
    CAR = 8,
    ANT = 0.0,
    EWT = 0.01,
    LMA = 0.009,
    N   = 1.5
  )
  
  tryCatch({
    res <- prosail::PRO4SAIL(
      Input_PROSPECT = Input_PROSPECT,
      TypeLidf       = 2,       # ellipsoidal LIDF
      LIDFa          = 30,      # mean leaf inclination angle
      lai            = lai,
      q              = 0.1,     # hot-spot parameter
      tts            = sza,     # solar zenith angle (degrees)
      tto            = vza,     # view  zenith angle (degrees)
      psi            = raa,     # relative azimuth  (degrees)
      rsoil          = prosail::SpecSOIL$Dry_Soil  # required soil spectrum
    )
    # $rso = bi-directional reflectance factor (BRF) — what remote sensors see
    return(as.numeric(res$rso))
  }, error = function(e) {
    message("PRO4SAIL error: ", conditionMessage(e))
    return(NULL)
  })
}

# ------------------------------------------------------------------------------
# compute_vi(): derive a vegetation index from a 2101-band reflectance vector.
# Band index = wavelength_nm - 400 + 1
# ------------------------------------------------------------------------------
compute_vi <- function(R, index_type) {
  if (is.null(R) || length(R) < 401) return(NA_real_)
  
  blue <- R[71]    # 470 nm
  red  <- R[271]   # 670 nm
  redE <- R[321]   # 720 nm
  nir  <- R[401]   # 800 nm
  
  denom_ndvi <- nir + red
  denom_evi  <- nir + 6 * red - 7.5 * blue + 1
  denom_ndre <- nir + redE
  
  switch(index_type,
         "ndvi" = if (denom_ndvi != 0) (nir - red)  / denom_ndvi else NA_real_,
         "evi"  = if (denom_evi  != 0) 2.5 * (nir - red) / denom_evi  else NA_real_,
         "ndre" = if (denom_ndre != 0) (nir - redE) / denom_ndre else NA_real_,
         NA_real_
  )
}

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  titlePanel("PROSAIL Hemispherical BRDF Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Location & Time"),
      numericInput("lat",  "Latitude:",  value = -31.95, step = 0.0001),
      numericInput("lng",  "Longitude:", value = 115.86, step = 0.0001),
      dateInput("date", "Date:", value = Sys.Date()),
      fluidRow(
        column(6, numericInput("hour12", "Hour (1-12):",
                               value = 12, min = 1, max = 12)),
        column(6, br(), checkboxInput("is_pm", "PM?", value = TRUE))
      ),
      
      hr(),
      h4("2. Model Parameters"),
      selectInput("index_type", "Vegetation Index:",
                  choices = c("NDVI" = "ndvi", "EVI" = "evi", "NDRE" = "ndre")),
      sliderInput("lai_val", "Leaf Area Index (LAI):",
                  min = 0.1, max = 7, value = 3, step = 0.1),
      sliderInput("cab", "Chlorophyll a+b (Cab, µg/cm2):",
                  min = 10, max = 100, value = 40),
      
      br(),
      actionButton("run_all", "Run Simulation",
                   class = "btn-primary", width = "100%")
    ),
    
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Map Selection",
                           br(),
                           p("Click the map to update Latitude/Longitude."),
                           leafletOutput("mymap", height = "500px")),
                  tabPanel("Polar VI Analysis",
                           plotOutput("polarPlot", height = "600px")),
                  tabPanel("Full Spectrum (Nadir)",
                           plotOutput("specPlot",  height = "500px"))
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================
server <- function(input, output, session) {
  
  nadir_ref <- reactiveVal(NULL)
  
  # ---------- Base map ---------------------------------------------------------
  output$mymap <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = 115.86, lat = -31.95, zoom = 5)
  })
  
  # Update lat/lng when user clicks map
  observeEvent(input$mymap_click, {
    click <- input$mymap_click
    updateNumericInput(session, "lat", value = round(click$lat, 4))
    updateNumericInput(session, "lng", value = round(click$lng, 4))
  })
  
  # ---------- Main simulation --------------------------------------------------
  observeEvent(input$run_all, {
    
    # -- Solar geometry ---------------------------------------------------------
    # FIX: 12 %% 12 == 0 breaks noon. Use explicit mapping instead.
    hour24 <- ifelse(input$hour12 == 12, 0L, as.integer(input$hour12)) +
      ifelse(input$is_pm, 12L, 0L)
    
    time_posix <- as.POSIXct(
      paste(input$date, sprintf("%02d:00:00", hour24)),
      tz = "UTC"
    )
    
    sun <- suncalc::getSunlightPosition(
      date = time_posix,
      lat  = input$lat,
      lon  = input$lng
    )
    
    sza_deg <- 90 - (sun$altitude * 180 / pi)   # altitude in radians -> SZA
    
    if (sza_deg >= 85) {
      showNotification(
        paste0("Sun too low for simulation (SZA = ", round(sza_deg, 1),
               "degrees). Choose a different time or location."),
        type = "warning", duration = 8
      )
      return()
    }
    
    # -- BRDF angular sampling grid --------------------------------------------
    vza_vec <- seq(0, 60, by = 10)    # 7 view-zenith angles
    raa_vec <- seq(0, 345, by = 15)   # 24 relative-azimuth angles
    grid    <- expand.grid(VZA = vza_vec, RAA = raa_vec)
    n_pts   <- nrow(grid)
    
    # -- Loop over grid with progress bar --------------------------------------
    withProgress(message = "Simulating BRDF...", value = 0, {
      
      vi_vals <- numeric(n_pts)
      
      for (i in seq_len(n_pts)) {
        incProgress(1 / n_pts,
                    detail = paste0(i, " / ", n_pts, " angular directions"))
        
        R <- run_pro4sail(
          cab = input$cab,
          lai = input$lai_val,
          sza = sza_deg,
          vza = grid$VZA[i],
          raa = grid$RAA[i]
        )
        
        vi_vals[i] <- compute_vi(R, input$index_type)
      }
      
      # Nadir reflectance for the spectrum tab (tto = 0, psi = 0)
      R_nadir <- run_pro4sail(
        cab = input$cab,
        lai = input$lai_val,
        sza = sza_deg,
        vza = 0,
        raa = 0
      )
      nadir_ref(R_nadir)
    })
    
    grid$VI_Value <- vi_vals
    
    # -- Polar BRDF plot -------------------------------------------------------
    output$polarPlot <- renderPlot({
      df <- na.omit(grid)
      validate(need(nrow(df) > 0, "No valid simulation results to display."))
      
      ggplot(df, aes(x = RAA, y = VZA, fill = VI_Value)) +
        geom_tile(width = 15, height = 10) +
        coord_polar(start = -pi / 2, direction = -1) +
        scale_fill_viridis_c(
          option = "magma",
          name   = toupper(input$index_type)
        ) +
        scale_x_continuous(
          breaks = seq(0, 315, by = 45),
          labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
        ) +
        labs(
          title    = paste("Hemispherical", toupper(input$index_type),
                           "Distribution"),
          subtitle = paste0("SZA: ", round(sza_deg, 1),
                            " deg  |  LAI: ", input$lai_val,
                            "  |  Cab: ", input$cab, " ug/cm2"),
          x = "Relative Azimuth Angle",
          y = "View Zenith Angle (degrees)"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title      = element_text(face = "bold"),
          legend.position = "right"
        )
    })
    
    # -- Nadir spectrum plot ---------------------------------------------------
    output$specPlot <- renderPlot({
      R <- nadir_ref()
      validate(need(!is.null(R) && length(R) >= 401,
                    "Nadir spectrum not available."))
      
      wvl <- seq(400, 400 + length(R) - 1)
      df  <- data.frame(Wavelength = wvl, Reflectance = as.numeric(R))
      
      ggplot(df, aes(x = Wavelength, y = Reflectance)) +
        geom_line(color = "darkgreen", linewidth = 1) +  # linewidth replaces deprecated size
        coord_cartesian(xlim = c(400, 1000)) +
        labs(
          title    = "Nadir Canopy Reflectance Spectrum",
          subtitle = paste0("LAI: ", input$lai_val,
                            "  |  Cab: ", input$cab, " ug/cm2"),
          x = "Wavelength (nm)",
          y = "Reflectance"
        ) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold"))
    })
    
    # Switch to results tab automatically
    updateTabsetPanel(session, "tabs", selected = "Polar VI Analysis")
  })
}

# ==============================================================================
shinyApp(ui, server)
