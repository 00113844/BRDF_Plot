import streamlit as st
import numpy as np
import pandas as pd
import prosail
import matplotlib.pyplot as plt
from datetime import datetime, date
import math
from streamlit_folium import st_folium
import folium

# Page configuration
st.set_page_config(layout="wide", page_title="PROSAIL BRDF Explorer")

# --- NOAA Solar Geometry Calculation ---
def calculate_solar_geometry(lat, lon, date_obj, hour, minute):
    """Manual calculation of solar zenith and azimuth based on NOAA equations."""
    doy = date_obj.timetuple().tm_yday
    time_dec = hour + minute / 60.0
    
    # Fractional year in radians
    gamma = 2 * math.pi / 365 * (doy - 1 + (time_dec - 12) / 24)
    
    # Equation of time (minutes) — Spencer (1971) corrected constant: 0.0000075 not 0.000075
    eqtime = 229.18 * (0.0000075 + 0.001868 * math.cos(gamma) - 0.032077 * math.sin(gamma) -
                      0.014615 * math.cos(2 * gamma) - 0.040849 * math.sin(2 * gamma))
    
    # Declination (radians) — Spencer (1971) full 3rd-order series
    decl = (0.006918 - 0.399912 * math.cos(gamma) + 0.070257 * math.sin(gamma)
            - 0.006758 * math.cos(2 * gamma) + 0.000907 * math.sin(2 * gamma)
            - 0.002697 * math.cos(3 * gamma) + 0.00148 * math.sin(3 * gamma))
    
    # True Solar Time
    time_offset = eqtime + 4 * lon
    tst = time_dec * 60 + time_offset
    ha = (tst / 4) - 180  # Hour angle in degrees
    
    lat_rad = math.radians(lat)
    ha_rad = math.radians(ha)
    
    # Zenith Angle
    cos_sza = (math.sin(lat_rad) * math.sin(decl) + 
               math.cos(lat_rad) * math.cos(decl) * math.cos(ha_rad))
    sza = math.degrees(math.acos(max(-1, min(1, cos_sza))))
    
    # Azimuth Angle (Clockwise from North)
    num = (math.sin(decl) - math.sin(lat_rad) * math.cos(math.radians(sza)))
    den = (math.cos(lat_rad) * math.sin(math.radians(sza)))
    
    if abs(den) < 1e-6:
        saa = 180.0
    else:
        saa = math.degrees(math.acos(max(-1, min(1, num / den))))
        if ha > 0:
            saa = 360 - saa
            
    # Return true SZA — let the caller handle nighttime / extreme angles
    return sza, saa

# --- Session State for map-driven lat/lon ---
# Apply any map click from the previous run BEFORE widgets are created
if st.session_state.get("_map_clicked"):
    st.session_state["lat_val"] = st.session_state.pop("_map_lat")
    st.session_state["lon_val"] = st.session_state.pop("_map_lon")
    st.session_state["_map_clicked"] = False

if "lat_val" not in st.session_state:
    st.session_state["lat_val"] = -31.95
if "lon_val" not in st.session_state:
    st.session_state["lon_val"] = 115.86

# --- UI Layout ---
st.title("PROSAIL BRDF & Spectrum Explorer")

with st.sidebar:
    st.header("Canopy Parameters")
    n_param = st.slider("N (leaf structure)", 1.0, 3.5, 1.5)
    cab = st.slider("cab (chlorophyll)", 0, 80, 40)
    car = st.slider("car (carotenoids)", 0, 30, 8)
    lai = st.slider("LAI (leaf area index)", 0.1, 8.0, 3.0)
    lidfa = st.slider("LIDFa (mean leaf angle °, 57=spherical)", 0.0, 90.0, 57.0)
    hspot = st.slider("hspot (hotspot)", 0.0, 1.0, 0.1)
    rsoil = st.slider("rsoil (soil brightness)", 0.0, 2.0, 1.0)
    psoil = st.slider("psoil (soil dryness 0=wet, 1=dry)", 0.0, 1.0, 1.0)
    
    st.header("Location & Time")
    input_date = st.date_input("Date (UTC)", date.today())
    c1, c2 = st.columns(2)
    hour = c1.number_input("Hour", 0, 23, 12)
    minute = c2.number_input("Min", 0, 59, 0)
    
    lat_val = st.number_input("Latitude", -90.0, 90.0, key="lat_val", step=0.01)
    lon_val = st.number_input("Longitude", -180.0, 180.0, key="lon_val", step=0.01)
    
    vi_choice = st.selectbox("Vegetation Index", ["NDVI", "NDRE", "SAVI"])
    run_btn = st.button("Run Simulation", type="primary", use_container_width=True)

# --- Logic & Simulation ---
actual_sza, saa = calculate_solar_geometry(lat_val, lon_val, input_date, hour, minute)
# Safe cap for PROSAIL — model is not validated beyond ~85°
prosail_sza = min(actual_sza, 85.0)

if actual_sza >= 90.0:
    st.sidebar.error(f"Sun below horizon (SZA = {actual_sza:.1f}°). Adjust time or location.")
elif actual_sza > 80.0:
    st.sidebar.warning(f"High solar zenith: {actual_sza:.2f}° — results may be unreliable.")
else:
    st.sidebar.success(f"Solar Zenith: {actual_sza:.2f}°")
st.sidebar.markdown(f"**Solar Azimuth:** `{saa:.2f}°` (N-relative)")

tab_map, tab_polar, tab_spec = st.tabs(["Location Map", "Polar BRDF Plot", "Full Spectrum"])

with tab_map:
    st.write("Click the map to select coordinates for solar geometry calculation.")
    m = folium.Map(location=[lat_val, lon_val], zoom_start=4)
    m.add_child(folium.LatLngPopup())
    map_data = st_folium(m, height=450, width=800)
    
    if map_data and map_data.get("last_clicked"):
        clicked_lat = map_data["last_clicked"]["lat"]
        clicked_lng = map_data["last_clicked"]["lng"]
        if (abs(clicked_lat - st.session_state["lat_val"]) > 1e-6 or
                abs(clicked_lng - st.session_state["lon_val"]) > 1e-6):
            st.session_state["_map_lat"] = clicked_lat
            st.session_state["_map_lon"] = clicked_lng
            st.session_state["_map_clicked"] = True
            st.rerun()
        st.success(f"Location: {clicked_lat:.4f}°, {clicked_lng:.4f}°")

if run_btn:
    if actual_sza >= 90.0:
        st.error(f"Cannot simulate: the sun is below the horizon (SZA = {actual_sza:.1f}°). "
                 "Select a daytime hour for the chosen location, or change the location.")
    else:
        with st.spinner("Simulating Radiative Transfer..."):
            # 1. Full Spectrum at Nadir
            # Python prosail: n, cab, car, cbrown, cw, cm, lai, lidfa, hspot, tts, tto, psi
            # Wavelengths 400-2500nm (2101 bands)
            refl_nadir = prosail.run_prosail(n_param, cab, car, 0.0, 0.015, 0.009, lai, lidfa, hspot, prosail_sza, 0, 0,
                                            rsoil=rsoil, psoil=psoil)
            wvl = np.arange(400, 2501)

            # 2. BRDF Grid
            vzas = np.arange(0, 81, 5)
            raas = np.arange(0, 361, 10)
            V, R = np.meshgrid(raas, vzas)
            # V[i,j] = raas[j] (azimuth), R[i,j] = vzas[i] (VZA)
            raa_f = V.flatten()
            vza_f = R.flatten()
            vi_vals = []

            for v, r in zip(vza_f, raa_f):
                res = prosail.run_prosail(n_param, cab, car, 0.0, 0.015, 0.009, lai, lidfa, hspot, prosail_sza, v, r,
                                         rsoil=rsoil, psoil=psoil)
                # Band indices (0-based): 670nm=270, 800nm=400, 705nm=305, 470nm=70
                red, nir, rededge, blue = res[270], res[400], res[305], res[70]

                if vi_choice == "NDVI":
                    vi_vals.append((nir - red) / (nir + red))
                elif vi_choice == "NDRE":
                    vi_vals.append((nir - rededge) / (nir + rededge))
                elif vi_choice == "SAVI":
                    vi_vals.append(1.5 * (nir - red) / (nir + red + 0.5))

            Z = np.array(vi_vals).reshape(V.shape)

            with tab_polar:
                fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
                ax.set_theta_zero_location("N")
                ax.set_theta_direction(-1)  # Clockwise

                c = ax.pcolormesh(np.radians(V), R, Z, cmap='magma', shading='auto')
                plt.colorbar(c, label=vi_choice, ax=ax, pad=0.1)

                ax.plot(np.radians(saa), actual_sza, marker='o', color='orange', markersize=12,
                        label='Sun Position', markeredgecolor='white')

                ax.set_title(f"{vi_choice} Angular Variation (SZA={prosail_sza:.1f}°)", pad=20, fontsize=15)
                ax.set_rlabel_position(135)
                if actual_sza > 80:
                    st.warning(f"Note: actual SZA={actual_sza:.1f}°, capped to {prosail_sza:.1f}° for PROSAIL.")
                st.pyplot(fig)

            with tab_spec:
                fig_s, ax_s = plt.subplots(figsize=(10, 5))
                ax_s.plot(wvl, refl_nadir, color='forestgreen', lw=2)
                ax_s.set_xlabel("Wavelength (nm)")
                ax_s.set_ylabel("Reflectance")
                ax_s.set_title("Nadir Spectral Signature")
                ax_s.grid(True, alpha=0.3)
                st.pyplot(fig_s)
else:
    with tab_polar:
        st.info("Adjust parameters in the sidebar and click **Run Simulation** to generate the BRDF polar plot.")
    with tab_spec:
        st.info("Adjust parameters in the sidebar and click **Run Simulation** to generate the spectrum plot.")