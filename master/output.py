# -*- coding: utf-8 -*-
from pcraster.framework import *
from copy import deepcopy

save_path = "tss\\output\\"

# HYDRO
def defineHydroTSS(model):
    # Rain
    model.resW_accRain_m3_tss = TimeoutputTimeseries("resW_accRain_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # GA infiltration
    model.resW_accGAinf_m3_tss = TimeoutputTimeseries("resW_accGAinf_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # GA runoff
    model.resW_accGAroff_m3_tss = TimeoutputTimeseries("resW_accGAroff_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # Runoff
    model.resW_accRunoff_m3_tss = TimeoutputTimeseries("resW_accRunoff_m3", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # z0 Percolation
    model.resW_accDPz0_m3_tss = TimeoutputTimeseries("resW_accDPz0_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # z1 Percolation
    model.resW_accDPz1_m3_tss = TimeoutputTimeseries("resW_accDPz1_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.tot_runoff_m3_tss = TimeoutputTimeseries("resW_totRunoff_m3", model, nominal("maps\\static\\outlet"),
                                                   noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # Percolation Basement
    model.resW_accPercol_Bsmt_m3_tss = TimeoutputTimeseries("resW_accPercol_Bsmt_m3", model,
                                                            nominal("maps\\static\\outlet"),
                                                            noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # Deep percolation Basement
    model.tot_perc_z3_m3_tss = TimeoutputTimeseries("resW_totPercol_z3_m3", model, nominal("maps\\static\\outlet"),
                                                    noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # ETP
    model.resW_accEtp_m3_tss = TimeoutputTimeseries("resW_accEtp_m3", model, nominal("maps\\static\\outlet"),
                                                    noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accEvap_m3_tss = TimeoutputTimeseries("resW_accEvap_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accTransp_m3_tss = TimeoutputTimeseries("resW_accTransp_m3", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.tot_etp_m3_tss = TimeoutputTimeseries("resW_totEtp_m3", model, nominal("maps\\static\\outlet"),
                                                noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # Baseflow
    model.out_baseflow_m3_tss = TimeoutputTimeseries("resW_accBaseflow_m3", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.tot_baseflow_m3_tss = TimeoutputTimeseries("resW_totBaseflow_m3", model, nominal("maps\\static\\outlet"),
    #                                                 noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # LF Drainage
    model.resW_accDrain_m3_tss = TimeoutputTimeseries("resW_accDrain_m3", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_o_cumDrain_m3_tss = TimeoutputTimeseries("resW_o_cumDrain_m3", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Cumulative ADR
    # LF options
    model.sat_accu_overflow_m3_tss = TimeoutputTimeseries("resW_of_accLatflow_m3", model,
                                                          nominal("maps\\static\\outlet"),
                                                          noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.tot_accu_of_latflow_m3_tss = TimeoutputTimeseries("resW_of_totLatflow_m3", model,
                                                            nominal("maps\\static\\outlet"),
                                                            noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # Inflow
    model.out_cell_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_cellLatflow_m3", model,
                                                           nominal("maps\\static\\outlet"),
                                                           noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.out_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_accLatflow_m3", model,
                                                           nominal("maps\\static\\outlet"),
                                                           noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.tot_accu_i_latflow_m3_tss = TimeoutputTimeseries("resW_i_totLatflow_m3", model,
                                                           nominal("maps\\static\\outlet"),
                                                           noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # Outflow
    model.resW_outLatflow_m3_tss = TimeoutputTimeseries("resW_outLatflow_m3", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Outlet LF
    model.resW_o_cumLatflow_m3_tss = TimeoutputTimeseries("resW_o_cumLatflow_m3", model,
                                                          nominal("maps\\static\\outlet"),
                                                          noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resW_accChStorage_m3_tss = TimeoutputTimeseries("resW_accChStorage_m3", model,
                                                          nominal("maps\\static\\outlet"),
                                                          noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.global_mb_water_tss = TimeoutputTimeseries("resW_global_waterMB", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.storage_m3_tss = TimeoutputTimeseries("resW_accStorage_m3", model, nominal("maps\\static\\outlet"),
                                                noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # Basement layer analysis
    # model.resW_accBAL_z3_m3_tss = TimeoutputTimeseries("resW_accBAL_z3_m3", model, nominal("maps\\static\\outlet"),
    #                                                   noHeader=False)
    # model.resW_accInfil_z3_m3_tss = TimeoutputTimeseries("resW_accInfil_z3_m3", model, nominal("maps\\static\\outlet"),
    #                                                     noHeader=False)
    # model.resW_accLF_z3_m3_tss = TimeoutputTimeseries("resW_accLF_z3_m3", model, nominal("maps\\static\\outlet"),
    #                                                  noHeader=False)
    # model.resW_accETP_z3_m3_tss = TimeoutputTimeseries("resW_accETP_z3_m3", model, nominal("maps\\static\\outlet"),
    #                                                   noHeader=False)

    model.resW_accBAL_Bsmt_m3_tss = TimeoutputTimeseries("resW_accBAL_Bsmt_m3", model,
                                                         nominal("maps\\static\\outlet"),
                                                         noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accInfil_Bsmt_m3_tss = TimeoutputTimeseries("resW_accInfil_Bsmt_m3", model,
                                                           nominal("maps\\static\\outlet"),
                                                           noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accLF_Bsmt_m3_tss = TimeoutputTimeseries("resW_accLF_Bsmt_m3", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accETP_Bsmt_m3_tss = TimeoutputTimeseries("resW_accETP_Bsmt_m3", model,
                                                         nominal("maps\\static\\outlet"),
                                                         noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # Storage
    model.resW_accVOL_z0_m3_tss = TimeoutputTimeseries("resW_accVOL_z0_m3", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accVOL_z1_m3_tss = TimeoutputTimeseries("resW_accVOL_z1_m3", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accVOL_z2_m3_tss = TimeoutputTimeseries("resW_accVOL_z2_m3", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accVOL_z3_m3_tss = TimeoutputTimeseries("resW_accVOL_z3_m3", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accVOL_Bsmt_m3_tss = TimeoutputTimeseries("resW_accVOL_Bsmt_m3", model,
                                                         nominal("maps\\static\\outlet"),
                                                         noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # Analysis
    # This is 'q' as time series.
    model.i_Q_m3_tss = TimeoutputTimeseries("resW_i_accVol_m3", model, nominal("maps\\static\\outlet"),
                                            noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.o_Q_m3_tss = TimeoutputTimeseries("resW_o_accVol_m3", model, nominal("maps\\static\\outlet"),
                                            noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_accQ_m3_tss = TimeoutputTimeseries("resW_accQ_m3", model, nominal("maps\\static\\outlet"),
                                                  noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.q_obs_cum_tss = TimeoutputTimeseries("resW_cum_q_obs_m3", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Equivalent to net_Q
    model.rain_obs_cum_tss = TimeoutputTimeseries("resW_cum_rain_obs_m3", model, nominal("maps\\static\\outlet"),
                                                  noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Equivalent to net_Q
    model.rest_obs_tss = TimeoutputTimeseries("resW_q_restit_obs_m3", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # = rain/q_obs
    model.q_sim_cum_tss = TimeoutputTimeseries("resW_cum_q_sim_m3", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Sum sim discharge (if obs available).
    model.q_sim_ave_tss = TimeoutputTimeseries("resW_q_sim_ave_m3", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # This is 'Nash_q' as time series.


def definePestTSS(model):
    # PESTI
    # Pesticide
    model.global_mb_pest_tss = TimeoutputTimeseries("resM_global_mb_pest", model, nominal("maps\\static\\outlet"),
                                                    noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accAPP_g_tss = TimeoutputTimeseries("resM_accAPP", model, nominal("maps\\static\\outlet"),
                                                   noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accVOLATz0_tss = TimeoutputTimeseries("resM_accVOLATz0", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accVOLATz0nor_tss = TimeoutputTimeseries("resM_accVOLATz0nor", model, nominal("maps\\sampling\\north_ave"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accVOLATz0val_tss = TimeoutputTimeseries("resM_accVOLATz0val", model, nominal("maps\\sampling\\valley_ave"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accVOLATz0sou_tss = TimeoutputTimeseries("resM_accVOLATz0sou", model, nominal("maps\\sampling\\south_ave"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_accROz0_tss = TimeoutputTimeseries("resM_accROz0", model, nominal("maps\\static\\outlet"),
                                                  noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accROz0nor_tss = TimeoutputTimeseries("resM_accROz0nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accROz0val_tss = TimeoutputTimeseries("resM_accROz0val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accROz0sou_tss = TimeoutputTimeseries("resM_accROz0sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_accDEGz0_tss = TimeoutputTimeseries("resM_accDEGz0", model, nominal("maps\\static\\outlet"),
                                                   noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accDEGzX_tss = TimeoutputTimeseries("resM_accDEGzX", model, nominal("maps\\static\\outlet"),
                                                   noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accDEGz0nor_tss = TimeoutputTimeseries("resM_accDEGz0nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accDEGz0val_tss = TimeoutputTimeseries("resM_accDEGz0val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accDEGz0sou_tss = TimeoutputTimeseries("resM_accDEGz0sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_accAGEDz0_tss = TimeoutputTimeseries("resM_accAGEDz0", model, nominal("maps\\static\\outlet"),
                                                    noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accAGEDzX_tss = TimeoutputTimeseries("resM_accAGEDzX", model, nominal("maps\\static\\outlet"),
                                                    noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accAGED_DEGz0_tss = TimeoutputTimeseries("resM_accAGED_DEGz0", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accAGED_DEGzX_tss = TimeoutputTimeseries("resM_accAGED_DEGzX", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_accLCHz0_tss = TimeoutputTimeseries("resM_accLCHz0", model, nominal("maps\\static\\outlet"),
                                                   noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accLCHz0nor_tss = TimeoutputTimeseries("resM_accLCHz0nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accLCHz0val_tss = TimeoutputTimeseries("resM_accLCHz0val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accLCHz0sou_tss = TimeoutputTimeseries("resM_accLCHz0sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_accLCHz1_tss = TimeoutputTimeseries("resM_accLCHz1", model, nominal("maps\\static\\outlet"),
                                                   noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_accDP_tss = TimeoutputTimeseries("resM_accDP", model, nominal("maps\\static\\outlet"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accADR_tss = TimeoutputTimeseries("resM_accADR", model, nominal("maps\\static\\outlet"),
                                                 noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accLF_tss = TimeoutputTimeseries("resM_accLF", model, nominal("maps\\static\\outlet"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accBF_tss = TimeoutputTimeseries("resM_accBF", model, nominal("maps\\static\\outlet"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accCHS_tss = TimeoutputTimeseries("resM_accCHS", model, nominal("maps\\static\\outlet"),
                                                 noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_accCHS_AGED_tss = TimeoutputTimeseries("resM_accCHS_AGED", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resM_EXP_light_g_tss = TimeoutputTimeseries("resM_EXP_light_g", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Total outlet mass (g) exports (light fraction only)
    model.resM_EXP_heavy_g_tss = TimeoutputTimeseries("resM_EXP_heavy_g", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Total outlet mass (g) exports (heavy fraction only)
    # Concentrations outlet
    model.resM_oCONC_ugL_tss = TimeoutputTimeseries("resM_oCONC_ugL", model, nominal("maps\\static\\outlet"),
                                                    noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Total outlet conc (ug/L)
    model.resM_oCONC_ROFF_ugL_tss = TimeoutputTimeseries("resM_oCONC_ROFF_ugL", model,
                                                         nominal("maps\\static\\outlet"),
                                                         noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Runoff outlet conc (ug/L)
    model.resM_oCONC_LF_ugL_tss = TimeoutputTimeseries("resM_oCONC_LF_ugL", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Latflow outlet conc (ug/L)

    model.resM_oCONC_ADR_ugL_tss = TimeoutputTimeseries("resM_oCONC_ADR_ugL", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Artificial drainage outlet conc (ug/L)
    # Isotopes outlet
    model.resM_outISO_d13C_tss = TimeoutputTimeseries("resM_outISO_d13C", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  #
    model.resM_outISO_ROFF_d13C_tss = TimeoutputTimeseries("resM_outISO_ROFF_d13C", model,
                                                           nominal("maps\\static\\outlet"),
                                                           noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Runoff outlet
    model.resM_outISO_LF_d13C_tss = TimeoutputTimeseries("resM_outISO_LF_d13C", model,
                                                         nominal("maps\\static\\outlet"),
                                                         noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Latflow outlet
    model.resM_outISO_ADR_d13C_tss = TimeoutputTimeseries("resM_outISO_ADR_d13C", model,
                                                          nominal("maps\\static\\outlet"),
                                                          noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Artificial drainage outlet

    # Cumulative Pesticide
    model.cum_degZ0_g_tss = TimeoutputTimeseries("resM_cumDEGz0", model, nominal("maps\\static\\outlet"),
                                                 noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Deg z0
    model.cum_deg_L_g_tss = TimeoutputTimeseries("resM_cumDEG_L", model, nominal("maps\\static\\outlet"),
                                                 noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Deg z0
    model.cum_aged_deg_L_g_tss = TimeoutputTimeseries("resM_cumAGE_DEG_L", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resM_cumLCHz0_L_g_tss = TimeoutputTimeseries("resM_cumLCHz0_L", model, nominal("maps\\static\\outlet"),
                                                       noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Leaching z0
    model.cum_roZ0_L_g_tss = TimeoutputTimeseries("resM_cumROz0_L", model, nominal("maps\\static\\outlet"),
                                                  noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Runoff

    model.cum_volatZ0_L_g_tss = TimeoutputTimeseries("resM_cumVOLATz0_L", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Runoff

    model.cum_adr_L_g_tss = TimeoutputTimeseries("resM_cumADR_L", model, nominal("maps\\static\\outlet"),
                                                 noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Art. drainage
    model.cum_latflux_L_g_tss = TimeoutputTimeseries("resM_cumLF_L", model, nominal("maps\\static\\outlet"),
                                                     noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Soil column, outlet cells

    model.resM_cumEXP_Smet_g_tss = TimeoutputTimeseries("resM_cumEXP_Smet_g", model, nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)  # Total cum. outlet mass (g) exports


def getLayerAnalysis(model, layer,
                     percolation, latflow_out, evap, transp,
                     root_depth, out_baseflow_m3=None):
    # Report Basement layer fluxes
    # Infil
    infil_m3 = percolation[layer - 1] * cellarea() / 1000  # m3
    infil_m3 = areatotal(infil_m3, model.is_catchment)
    # Latflow
    latflow_m3 = latflow_out[layer] * cellarea() / 1000
    latflow_m3 = areatotal(latflow_m3, model.outlet_multi)  # Only outlet cells
    # Baseflow
    # Already recorded...
    if out_baseflow_m3 is None:
        out_baseflow_m3 = deepcopy(model.zero_map)
    # Evapotransp
    evap_m3 = evap[layer] * cellarea() / 1000  # m3
    transp_m3 = transp[layer] * cellarea() / 1000  # m3
    evapotransp_m3 = evap_m3 + transp_m3
    # model.report(evap_m3, 'z' + str(layer) + 'EVA')  # Check which cells
    # model.report(transp_m3, 'z' + str(layer) + 'TRA')  # Check which cells
    # model.report(root_depth[layer], 'z' + str(layer) + 'ROOT')
    evapotransp_m3 = areatotal(evapotransp_m3, model.is_catchment)
    # Storage
    storage_m3 = model.theta[layer] * model.layer_depth[layer] * cellarea() / 1000
    storage_m3 = areatotal(storage_m3, model.is_catchment)

    # Change In storage
    ch_storage_m3 = ((model.theta[layer] * model.layer_depth[layer] * cellarea() / 1000) -
                     (model.theta_ini[layer] * model.layer_depth[layer] * cellarea() / 1000))
    ch_storage_m3 = areatotal(ch_storage_m3, model.is_catchment)
    balance_m3 = infil_m3 - latflow_m3 - out_baseflow_m3 - evapotransp_m3 - ch_storage_m3

    if layer == 0:
        model.resW_accVOL_z0_m3_tss.sample(storage_m3)
    elif layer == 1:
        model.resW_accVOL_z1_m3_tss.sample(storage_m3)
    elif layer == 2:
        model.resW_accVOL_z2_m3_tss.sample(storage_m3)
    elif layer == 3:
        # model.resW_accInfil_z3_m3_tss.sample(infil_m3)
        # model.resW_accLF_z3_m3_tss.sample(latflow_m3)
        # model.resW_accETP_z3_m3_tss.sample(evapotransp_m3)
        model.resW_accVOL_z3_m3_tss.sample(storage_m3)
        # model.resW_accBAL_z3_m3_tss.sample(balance_m3)
    elif layer == (model.num_layers - 1):
        model.resW_accInfil_Bsmt_m3_tss.sample(infil_m3)
        model.resW_accLF_Bsmt_m3_tss.sample(latflow_m3)
        model.resW_accETP_Bsmt_m3_tss.sample(evapotransp_m3)
        model.resW_accVOL_Bsmt_m3_tss.sample(storage_m3)
        model.resW_accBAL_Bsmt_m3_tss.sample(balance_m3)
        # aguila --scenarios='{2}' --timesteps=[2,300,2] z3EVA z3TRA z3ROOT


def getCatchmentStorage(model):
    cell_vol_tot_m3 = deepcopy(model.zero_map)
    for layer in range(model.num_layers):
        cell_vol_tot_m3 += model.theta[layer] * model.layer_depth[layer] * cellarea() / 1000

    vol_tot_m3 = accuflux(model.ldd_subs, cell_vol_tot_m3)
    multi_vol_tot_m3 = areatotal(vol_tot_m3, model.outlet_multi)
    model.storage_m3_tss.sample(multi_vol_tot_m3)


def defineAverageMoistTSS(model):
    # Theta average proportion to saturation
    model.resW_z0_thetaPropSat = TimeoutputTimeseries("resW_z0_thetaPropSat", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resW_z1_thetaPropSat = TimeoutputTimeseries("resW_z1_thetaPropSat", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z2_thetaPropSat = TimeoutputTimeseries("resW_z2_thetaPropSat", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z3_thetaPropSat = TimeoutputTimeseries("resW_z3_thetaPropSat", model, nominal("maps\\static\\outlet"),
                                                      noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_Bsmt_thetaPropSat = TimeoutputTimeseries("resW_Bsmt_thetaPropSat", model,
                                                        nominal("maps\\static\\outlet"),
                                                        noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)


def defineTopSoilConditions(model):
    # Catchment Theta
    model.resW_z0_theta = TimeoutputTimeseries("resW_z0_theta", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z1_theta = TimeoutputTimeseries("resW_z1_theta", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z2_theta = TimeoutputTimeseries("resW_z2_theta", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z3_theta = TimeoutputTimeseries("resW_z3_theta", model, nominal("maps\\static\\outlet"),
                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # model.resW_z0_theta_max = TimeoutputTimeseries("resW_z0_theta_max", model, nominal("maps\\static\\outlet"),
    #                                            noHeader=False)
    # model.resW_z1_theta_max = TimeoutputTimeseries("resW_z1_theta_max", model, nominal("maps\\static\\outlet"),
    #                                            noHeader=False)
    # model.resW_z2_theta_max = TimeoutputTimeseries("resW_z2_theta_max", model, nominal("maps\\static\\outlet"),
    #                                            noHeader=False)

    # Transect theta
    # model.resW_z0_theta_nor = TimeoutputTimeseries("resW_z0_theta_nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_theta_val = TimeoutputTimeseries("resW_z0_theta_val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_theta_sou = TimeoutputTimeseries("resW_z0_theta_sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    #
    # model.resW_z1_theta_nor = TimeoutputTimeseries("resW_z1_theta_nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z1_theta_val = TimeoutputTimeseries("resW_z1_theta_val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z1_theta_sou = TimeoutputTimeseries("resW_z1_theta_sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    #
    # model.resW_z2_theta_nor = TimeoutputTimeseries("resW_z2_theta_nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z2_theta_val = TimeoutputTimeseries("resW_z2_theta_val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z2_theta_sou = TimeoutputTimeseries("resW_z2_theta_sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # Temperature
    model.resW_z0_temp = TimeoutputTimeseries("resW_z0_Temp", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resW_z1_temp = TimeoutputTimeseries("resW_z1_Temp", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resW_z2_temp = TimeoutputTimeseries("resW_z2_Temp", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    model.resW_z3_temp = TimeoutputTimeseries("resW_z3_Temp", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # model.resW_z0_temp_max = TimeoutputTimeseries("resW_z0_temp_max", model, nominal("maps\\static\\outlet"),
    #                                                noHeader=False)
    # model.resW_z1_temp_max = TimeoutputTimeseries("resW_z1_temp_max", model, nominal("maps\\static\\outlet"),
    #                                                noHeader=False)
    # model.resW_z2_temp_max = TimeoutputTimeseries("resW_z2_temp_max", model, nominal("maps\\static\\outlet"),
    #                                                noHeader=False)

    # Transect Temp
    # model.resW_z0_temp_nor = TimeoutputTimeseries("resW_z0_temp_nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_temp_val = TimeoutputTimeseries("resW_z0_temp_val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_temp_sou = TimeoutputTimeseries("resW_z0_temp_sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    #
    # model.resW_z1_temp_nor = TimeoutputTimeseries("resW_z1_temp_nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z1_temp_val = TimeoutputTimeseries("resW_z1_temp_val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z1_temp_sou = TimeoutputTimeseries("resW_z1_temp_sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    #
    # model.resW_z2_temp_nor = TimeoutputTimeseries("resW_z2_temp_nor", model, nominal("maps\\sampling\\north_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z2_temp_val = TimeoutputTimeseries("resW_z2_temp_val", model, nominal("maps\\sampling\\valley_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z2_temp_sou = TimeoutputTimeseries("resW_z2_temp_sou", model, nominal("maps\\sampling\\south_ave"), noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # DT50
    model.resW_z0_DT50 = TimeoutputTimeseries("resW_z0_DT50", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z1_DT50 = TimeoutputTimeseries("resW_z1_DT50", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    model.resW_z2_DT50 = TimeoutputTimeseries("resW_z2_DT50", model, nominal("maps\\static\\outlet"),
                                              noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # Transect DT50s
    # Reported in pesti_vXX.py module
    # model.resW_z0_DT50_min = TimeoutputTimeseries("resW_z0_DT50_min", model, nominal("maps\\static\\outlet"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_DT50_max = TimeoutputTimeseries("resW_z0_DT50_max", model, nominal("maps\\static\\outlet"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)

    # model.resW_z0_DT50_nor = TimeoutputTimeseries("resW_z0_DT50_nor", model, nominal("maps\\sampling\\north_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_DT50_val = TimeoutputTimeseries("resW_z0_DT50_val", model, nominal("maps\\sampling\\valley_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z0_DT50_sou = TimeoutputTimeseries("resW_z0_DT50_sou", model, nominal("maps\\sampling\\south_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    #
    # model.resW_z1_DT50_nor = TimeoutputTimeseries("resW_z1_DT50_nor", model, nominal("maps\\sampling\\north_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z1_DT50_val = TimeoutputTimeseries("resW_z1_DT50_val", model, nominal("maps\\sampling\\valley_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z1_DT50_sou = TimeoutputTimeseries("resW_z1_DT50_sou", model, nominal("maps\\sampling\\south_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    #
    # model.resW_z2_DT50_nor = TimeoutputTimeseries("resW_z2_DT50_nor", model, nominal("maps\\sampling\\north_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z2_DT50_val = TimeoutputTimeseries("resW_z2_DT50_val", model, nominal("maps\\sampling\\valley_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)
    # model.resW_z2_DT50_sou = TimeoutputTimeseries("resW_z2_DT50_sou", model, nominal("maps\\sampling\\south_ave"),
    #                                               noHeader=False, save_path=save_path, period=model.period, sample_nr=model.sample_nr)


def getTopSoilConditions(model, layer=0):
    # Catchment
    theta_ave = areaaverage(model.theta[layer], model.is_catchment)
    theta_ave_nor = areaaverage(model.theta[layer], model.is_north)
    theta_ave_val = areaaverage(model.theta[layer], model.is_valley)
    theta_ave_sou = areaaverage(model.theta[layer], model.is_south)

    prop_temp = model.temp_fin[layer] / 27.  # Max air temp
    temp_ave = areaaverage(prop_temp, model.is_catchment)
    temp_ave_nor = areaaverage(prop_temp, model.is_north)
    temp_ave_val = areaaverage(prop_temp, model.is_valley)
    temp_ave_sou = areaaverage(prop_temp, model.is_south)

    if layer == 0:
        model.resW_z0_theta.sample(theta_ave)
        # model.resW_z0_theta_nor.sample(theta_ave_nor)
        # model.resW_z0_theta_val.sample(theta_ave_val)
        # model.resW_z0_theta_sou.sample(theta_ave_sou)

        model.resW_z0_temp.sample(temp_ave)
        # model.resW_z0_temp_nor.sample(temp_ave_nor)
        # model.resW_z0_temp_val.sample(temp_ave_val)
        # model.resW_z0_temp_sou.sample(temp_ave_sou)

    elif layer == 1:
        model.resW_z1_theta.sample(theta_ave)
        # model.resW_z1_theta_nor.sample(theta_ave_nor)
        # model.resW_z1_theta_val.sample(theta_ave_val)
        # model.resW_z1_theta_sou.sample(theta_ave_sou)

        model.resW_z1_temp.sample(temp_ave)
        # model.resW_z1_temp_nor.sample(temp_ave_nor)
        # model.resW_z1_temp_val.sample(temp_ave_val)
        # model.resW_z1_temp_sou.sample(temp_ave_sou)

    elif layer == 2:
        model.resW_z2_theta.sample(theta_ave)
        # model.resW_z2_theta_nor.sample(theta_ave_nor)
        # model.resW_z2_theta_val.sample(theta_ave_val)
        # model.resW_z2_theta_sou.sample(theta_ave_sou)

        model.resW_z2_temp.sample(temp_ave)
        # model.resW_z2_temp_nor.sample(temp_ave_nor)
        # model.resW_z2_temp_val.sample(temp_ave_val)
        # model.resW_z2_temp_sou.sample(temp_ave_sou)


def getAverageMoisture(model):
    prop_sat_layers = []
    for layer in range(model.num_layers):
        prop_sat = model.theta[layer] / model.theta_sat[layer]
        prop_sat_ave = areaaverage(prop_sat, model.is_catchment)
        prop_sat_layers.append(prop_sat_ave)

    model.resW_z0_thetaPropSat.sample(prop_sat_layers[0])
    model.resW_z1_thetaPropSat.sample(prop_sat_layers[1])
    model.resW_z2_thetaPropSat.sample(prop_sat_layers[2])
    model.resW_z3_thetaPropSat.sample(prop_sat_layers[3])
    model.resW_Bsmt_thetaPropSat.sample(prop_sat_layers[-1])


def reportCumHydro(model, q_obs, out_runoff_m3, out_drain_m3, tot_rain_m3,
                   out_etp_m3, outlet_latflow_m3, out_percol_m3=None):
    # Don't need it for GLUE
    model.rain_cum_m3 += tot_rain_m3
    model.rain_obs_cum_tss.sample(model.rain_cum_m3)

    # Mass balance components
    model.tot_runoff_m3 += ifthenelse(q_obs >= 0, out_runoff_m3, 0)
    model.tot_runoff_m3_tss.sample(model.tot_runoff_m3)

    if out_percol_m3 is not None:
        model.tot_perc_z3_m3 += ifthenelse(q_obs >= 0, out_percol_m3, 0)  # <- == 0
        model.tot_perc_z3_m3_tss.sample(model.tot_perc_z3_m3)

    model.tot_etp_m3 += ifthenelse(q_obs >= 0, out_etp_m3, 0)
    model.tot_etp_m3_tss.sample(model.tot_etp_m3)
    # model.tot_baseflow_m3 += ifthenelse(q_obs >= 0, out_baseflow_m3, 0)
    # model.tot_baseflow_m3_tss.sample(model.tot_baseflow_m3)

    # Cumulative drainage
    model.tot_drain_m3 += ifthenelse(q_obs >= 0, out_drain_m3, 0)  # o_drain_z1_m3
    model.resW_o_cumDrain_m3_tss.sample(model.tot_drain_m3)

    # model.tot_ilf_m3 += ifthenelse(q_obs >= 0, outlet_lat_inflow_m3, 0)
    # model.tot_accu_i_latflow_m3_tss.sample(model.tot_ilf_m3)
    model.cum_olf_m3 += ifthenelse(q_obs >= 0, outlet_latflow_m3, 0)
    model.resW_o_cumLatflow_m3_tss.sample(model.cum_olf_m3)
    # model.tot_nlf_m3 += ifthenelse(q_obs >= 0, n_latflow_m3, 0)
    # model.tot_accu_n_latflow_m3_tss.sample(model.tot_nlf_m3)
    # model.tot_of_m3 += ifthenelse(q_obs >= 0, of_latflow_m3, 0)
    # model.tot_accu_of_latflow_m3_tss.sample(model.tot_of_m3)


def reportGlobalWaterBalance(model, tot_rain_m3, out_runoff_m3, out_drain_m3,
                             outlet_latflow_m3,
                             out_etp_m3, accu_ch_storage_m3,
                             GA_add_m3,
                             GA_runoff_m3,
                             out_percol_m3=None,
                             out_baseflow_m3=None):
    model.resW_accRain_m3_tss.sample(tot_rain_m3)
    model.resW_accGAinf_m3_tss.sample(GA_add_m3)
    model.resW_accGAroff_m3_tss.sample(scalar(GA_runoff_m3))
    model.resW_accEtp_m3_tss.sample(out_etp_m3)
    model.resW_accRunoff_m3_tss.sample(out_runoff_m3)  # save to outlet
    model.resW_accDrain_m3_tss.sample(out_drain_m3)  # Outlet discharge - Drain
    model.resW_outLatflow_m3_tss.sample(outlet_latflow_m3)  # Outlet discharge - LF
    model.resW_accChStorage_m3_tss.sample(accu_ch_storage_m3)

    if out_percol_m3 is not None:
        model.resW_accPercol_Bsmt_m3_tss.sample(out_percol_m3)  # Basement

    if out_baseflow_m3 is not None:
        model.out_baseflow_m3_tss.sample(out_baseflow_m3)
    # model.out_accu_o_latflow_m3_tss.sample(catch_lat_outflow_m3)
    # model.resW_accChStorage_m3_tss.sample(out_ch_storage_m3)

    # GLOBAL Water
    if out_baseflow_m3 is not None:
        global_mb_water = (tot_rain_m3 + GA_add_m3 + scalar(GA_runoff_m3) -
                           out_etp_m3 - out_runoff_m3 -
                           # out_percol_m3 - # Basement percolation
                           outlet_latflow_m3 -
                           out_drain_m3 -
                           accu_ch_storage_m3 -
                           out_baseflow_m3
                           )
    else:
        global_mb_water = (tot_rain_m3 + GA_add_m3 + scalar(GA_runoff_m3) -
                           out_etp_m3 - out_runoff_m3 -
                           # out_percol_m3 -  # Basement percolation
                           outlet_latflow_m3 -
                           out_drain_m3 -
                           accu_ch_storage_m3
                           )

    model.global_mb_water_tss.sample(global_mb_water)


def reportGlobalPestBalance(model,
                            catch_app_z0,
                            z0_deg_catch,
                            z0_aged_deg_catch,
                            catch_volat,
                            catch_runoff_mass,
                            z0_catch_leach,
                            # catch_drain_light,
                            z0_latflux,
                            z0_ch_storage,
                            z0_ch_storage_aged):
    z0_mb_pest = (catch_app_z0 -
                  z0_deg_catch -
                  z0_aged_deg_catch -
                  catch_volat -
                  catch_runoff_mass -
                  z0_catch_leach -
                  # catch_drain_light -
                  z0_latflux -  #
                  z0_ch_storage -
                  z0_ch_storage_aged)

    model.global_mb_pest_tss.sample(z0_mb_pest)


def repCumOutMass(model, conc_outlet_obs, outlet_light_export,
                  catch_latflux_light, catch_drain_light, catch_runoff_light,
                  catch_volat_light, catch_deg_light):
    # Outlet-specific, Cumulative masses
    model.cum_exp_L_g += ifthenelse(conc_outlet_obs > 0, outlet_light_export, scalar(0))
    model.resM_cumEXP_Smet_g_tss.sample(model.cum_exp_L_g)

    model.cum_latflux_L_g += ifthenelse(conc_outlet_obs > 0, catch_latflux_light, scalar(0))
    model.cum_latflux_L_g_tss.sample(model.cum_latflux_L_g)

    model.cum_adr_L_g += ifthenelse(conc_outlet_obs > 0, catch_drain_light, scalar(0))
    model.cum_adr_L_g_tss.sample(model.cum_adr_L_g)

    model.cum_roZ0_L_g += ifthenelse(conc_outlet_obs > 0, catch_runoff_light, scalar(0))
    model.cum_roZ0_L_g_tss.sample(model.cum_roZ0_L_g)

    # Not in outlet, but relevant cumulative sinks
    model.cum_volatZ0_L_g += catch_volat_light
    model.cum_volatZ0_L_g_tss.sample(model.cum_volatZ0_L_g)

    model.cum_deg_L_g += catch_deg_light
    model.cum_deg_L_g_tss.sample(model.cum_deg_L_g)
