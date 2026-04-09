from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import folium
import io
import json

from shiny import App, ui, render, reactive

try:
    from suncalc import get_position
except Exception:
    get_position = None

try:
    import prosail
except Exception:
    prosail = None


def compute_solar_angles(
    manual: bool,
    sza: float,
    saa: float,
    date_val,
    hour: int,
    minute: int,
    lat: float,
    lon: float,
) -> Tuple[float, float]:
    if manual or get_position is None:
        return float(sza), float(saa)

    dt = datetime(
        date_val.year,
        date_val.month,
        date_val.day,
        int(hour),
        int(minute),
        tzinfo=timezone.utc,
    )
    pos = get_position(dt, lat, lon)
    altitude_deg = np.degrees(pos["altitude"])
    azimuth_deg = np.degrees(pos["azimuth"])
    sza_out = float(np.clip(90.0 - altitude_deg, 0.0, 90.0))
    saa_out = float((azimuth_deg + 180.0) % 360.0)
    return sza_out, saa_out


def run_prosail_spectrum(params: Dict[str, float], sza: float, vza: float, raa: float) -> np.ndarray:
    if prosail is None:
        raise RuntimeError("Python package 'prosail' is not installed.")

    try:
        # Call PROSAIL with explicit keyword arguments for rsoil and psoil
        return prosail.run_prosail(
            params["N"],           # N (positional)
            params["cab"],         # cab (positional)
            params["car"],         # car (positional)
            params["cbrown"],      # cbrown (positional)
            params["cw"],          # cw (positional)
            params["cm"],          # cm (positional)
            params["lai"],         # lai (positional)
            params["lidfa"],       # lidfa (positional)
            params["hspot"],       # hspot (positional)
            sza,                   # tts (positional)
            vza,                   # tto (positional)
            raa,                   # psi (positional)
            rsoil=params["rsoil"], # rsoil (keyword)
            psoil=params["psoil"], # psoil (keyword)
        )
    except ValueError as e:
        raise RuntimeError(f"PROSAIL error: {str(e)}")


def compute_vi(refl: np.ndarray, vi: str) -> float:
    def wl_idx(wl: int) -> int:
        return wl - 400

    blue = refl[wl_idx(470)]
    red = refl[wl_idx(670)]
    rededge = refl[wl_idx(705)]
    nir = refl[wl_idx(800)]

    if vi == "NDVI":
        return (nir - red) / (nir + red)
    if vi == "EVI":
        return 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
    if vi == "SAVI":
        return 1.5 * (nir - red) / (nir + red + 0.5)
    return (nir - rededge) / (nir + rededge)


app_ui = ui.page_fluid(
    ui.h2("PROSAIL BRDF Explorer (Python)"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4("Canopy parameters"),
            ui.input_numeric("N", "N (leaf structure)", 1.5, min=1.0, max=3.5, step=0.1),
            ui.input_numeric("cab", "cab (chlorophyll)", 40, min=0, max=80, step=1),
            ui.input_numeric("car", "car (carotenoids)", 8, min=0, max=30, step=1),
            ui.input_numeric("cbrown", "cbrown (brown pigment)", 0.0, min=0, max=1, step=0.01),
            ui.input_numeric("cw", "cw (water thickness)", 0.015, min=0, max=0.05, step=0.001),
            ui.input_numeric("cm", "cm (dry matter)", 0.009, min=0, max=0.03, step=0.001),
            ui.input_numeric("lai", "lai (leaf area index)", 3.0, min=0.1, max=8.0, step=0.1),
            ui.input_numeric("lidfa", "lidfa (leaf angle)", 50, min=0, max=90, step=1),
            ui.input_numeric("hspot", "hspot (hotspot)", 0.1, min=0, max=1, step=0.01),
            ui.h4("Soil"),
            ui.input_numeric("rsoil", "rsoil (brightness)", 1.0, min=0.0, max=2.0, step=0.05),
            ui.input_numeric("psoil", "psoil (moisture mix)", 0.5, min=0.0, max=1.0, step=0.05),
            ui.input_text("rsoil0", "rsoil0 (optional)", ""),
            ui.h4("Solar illumination"),
            ui.input_checkbox("manual_solar", "Manual override (sza/saa)", False),
            ui.input_date("date_utc", "📅 Date (UTC)", value=datetime.now(timezone.utc).date()),
            ui.layout_columns(
                ui.input_numeric("hour_utc", "⏰ Hour (UTC)", 12, min=0, max=23, step=1),
                ui.input_numeric("minute_utc", "⏰ Minute (UTC)", 0, min=0, max=59, step=1),
                col_widths=(6, 6),
            ),
            ui.h4("📍 Location (click map or enter manually)"),
            ui.input_numeric("lat", "Latitude", 0.0, min=-90, max=90, step=0.01),
            ui.input_numeric("lon", "Longitude", 0.0, min=-180, max=180, step=0.01),
            ui.input_numeric("sza", "sza (deg)", 30, min=0, max=90, step=1),
            ui.input_numeric("saa", "saa (deg)", 0, min=0, max=360, step=1),
            ui.h4("Geometry grid"),
            ui.input_numeric("vza_max", "vza_max (deg)", 70, min=0, max=85, step=1),
            ui.input_numeric("step_deg", "step_deg", 10, min=1, max=20, step=1),
            ui.h4("Vegetation index"),
            ui.input_select("vi", "Index", ["NDVI", "EVI", "SAVI", "NDRE"]),
            ui.input_action_button("run", "Run simulation"),
        ),
        ui.layout_columns(
            ui.div(
                ui.h4("BRDF Simulation"),
                ui.output_plot("brdf_plot", height="650px"),
                ui.output_text_verbatim("solar_text"),
            ),
            ui.div(
                ui.h4("📍 Location Reference Map"),
                ui.output_ui("location_map"),
            ),
            col_widths=(6, 6),
        ),
    ),
)


def server(input, output, session):
    @reactive.calc
    def solar_angles() -> Tuple[float, float]:
        return compute_solar_angles(
            input.manual_solar(),
            input.sza(),
            input.saa(),
            input.date_utc(),
            input.hour_utc(),
            input.minute_utc(),
            input.lat(),
            input.lon(),
        )

    @output
    @render.ui
    def location_map():
        """Generate a folium map displayed directly."""
        try:
            lat = float(input.lat())
            lon = float(input.lon())
        except:
            lat, lon = 0.0, 0.0
        
        try:
            # Create map centered on current location
            m = folium.Map(
                location=[lat, lon],
                zoom_start=5,
                tiles="CartoDB positron",
            )
            
            # Add marker at current location
            folium.Marker(
                location=[lat, lon],
                popup=f"<b>Coordinates:</b><br>Lat: {lat:.4f}<br>Lon: {lon:.4f}",
                tooltip="Click on map to update",
                icon=folium.Icon(color="blue", icon="info-sign"),
            ).add_to(m)
            
            # Save to HTML string
            html_buf = io.BytesIO()
            m.save(html_buf, close_file=False)
            html_str = html_buf.getvalue().decode('utf-8')
            
            # Extract just the body content and embed in divwith styling
            import re
            body_match = re.search(r'<body>(.*?)</body>', html_str, re.DOTALL)
            body_content = body_match.group(1) if body_match else html_str
            
            # Extract head content for scripts and styles
            head_match = re.search(r'<head>(.*?)</head>', html_str, re.DOTALL)
            head_content = head_match.group(1) if head_match else ''
            
            # Create wrapped HTML with click handler
            wrapped_html = f'''
            <div style="width: 100%; height: 400px; position: relative; overflow: hidden; border: 1px solid #ddd; border-radius: 5px;">
            {head_content}
            <style>
                .leaflet-container {{ height: 400px !important; }}
            </style>
            {body_content}
            <script>
            (function() {{
                function waitForMap(attempt = 0) {{
                    if (attempt > 50) return;
                    var mapDiv = document.querySelector('.leaflet-container');
                    if (!mapDiv || !mapDiv._leaflet_map) {{
                        setTimeout(() => waitForMap(attempt + 1), 100);
                        return;
                    }}
                    var map = mapDiv._leaflet_map;
                    map.on('click', function(e) {{
                        var newLat = parseFloat(e.latlng.lat.toFixed(4));
                        var newLon = parseFloat(e.latlng.lng.toFixed(4));
                        if (typeof Shiny !== 'undefined' && Shiny.setInputValue) {{
                            Shiny.setInputValue('lat', newLat, {{priority: 'event'}});
                            Shiny.setInputValue('lon', newLon, {{priority: 'event'}});
                        }}
                    }});
                }}
                waitForMap();
            }})();
            </script>
            </div>
            '''
            
            return ui.HTML(wrapped_html)
        except Exception as e:
            import traceback
            error_msg = traceback.format_exc()
            return ui.HTML(f"<p style='color:red;'>Map error: {error_msg}</p>")

    @reactive.event(input.run)
    def run_sim() -> pd.DataFrame:
        sza, saa = solar_angles()
        vza = np.arange(0, input.vza_max() + 0.1, input.step_deg())
        raa = np.arange(0, 360 + 0.1, input.step_deg())
        grid = pd.MultiIndex.from_product([vza, raa], names=["vza", "raa"]).to_frame(index=False)

        rsoil0_raw = input.rsoil0().strip()
        rsoil0_val = float(rsoil0_raw) if rsoil0_raw else None

        params = {
            "N": float(input.N()),
            "cab": float(input.cab()),
            "car": float(input.car()),
            "cbrown": float(input.cbrown()),
            "cw": float(input.cw()),
            "cm": float(input.cm()),
            "lai": float(input.lai()),
            "lidfa": float(input.lidfa()),
            "hspot": float(input.hspot()),
            "rsoil": float(input.rsoil()),
            "psoil": float(input.psoil()),
            "rsoil0": rsoil0_val,
        }

        vi_values = []
        for _, row in grid.iterrows():
            spectrum = run_prosail_spectrum(params, sza, row.vza, row.raa)
            vi_values.append(compute_vi(spectrum, input.vi()))

        grid["vi"] = vi_values
        grid["sza"] = sza
        grid["saa"] = saa
        return grid

    @output
    @render.plot
    def brdf_plot():
        df = run_sim()
        if df.empty:
            return None

        fig, ax = plt.subplots(figsize=(8, 7), subplot_kw={"projection": "polar"})
        raa_rad = np.radians(df["raa"].to_numpy())
        contour = ax.tricontourf(raa_rad, df["vza"].to_numpy(), df["vi"].to_numpy(), levels=30, cmap="YlGn")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks(np.radians(np.arange(0, 360, 45)))
        ax.set_xticklabels(["0 (Sun)", "45", "90", "135", "180", "225", "270", "315"])
        ax.set_ylim(0, max(df["vza"]))

        sza, _ = solar_angles()
        ax.plot(0, sza, marker="*", markersize=14, color="gold", markeredgecolor="black")
        fig.colorbar(contour, ax=ax, orientation="vertical", pad=0.1, shrink=0.7, label=input.vi())
        ax.set_title("PROSAIL BRDF " + input.vi(), pad=20)
        return fig

    @output
    @render.text
    def solar_text():
        sza, saa = solar_angles()
        return f"SZA: {sza:.2f} deg\nSAA: {saa:.2f} deg"


app = App(app_ui, server)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8000)
