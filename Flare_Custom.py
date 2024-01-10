


time, flux, detrend, error = tier0('/data/whitsett.n/TESS_Light_Curves/All_Exoplanets/2MASS J01033563-5515561 A/tess2020238165205-s0029-0000000206502540-0193-s_lc.fits')

d, e = tier1(detrend, 3)
tier2(time, flux, error, d, e, '/data/whitsett.n/', host_name = 'Test')
tier3('/data/whitsett.n/Test', '/data/whitsett.n/Test_working_Dir', '/data/whitsett.n/Test_output', '/data/whitsett.n/allesfitter/templates/N_Flares/settings.csv', '/data/whitsett.n/allesfitter/templates/N_Flares/params.csv', host_name = 'Test', T = 3000, host_radius = 1.2, MCMC_CPUS=63)





# TESS_Folder_ID = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/')
# TOI_Catalog = pd.read_csv('/data/whitsett.n/Reference_Files/All_TOI_12_17_23.csv')
# total_flares = 0
# total_possible_flares = 0
# total_observation_time = 0
# TOI_ID_list = []
# flare_number = []
# peak_time = []
# amplitude = []
# time_scale = []
# Teff = []
# radius = []
# flare_phase = []
# TOI_period = []
# total_flare_energies = []
# total_flare_phases = []
# item_count = 0
# total_periastron_list = []
# total_periastron_epoch_list = []
# total_epoch_list = []
# list_index = 0

# tier0_tau = []
# tier1_tau = []
# tier2_tau = []
# for M_dwarves in TESS_Folder_ID:
#     ##Iteration Scheme
#     TOI_ID = M_dwarves
#     print(TOI_ID)

#     a = os.listdir('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves)
#     # os.mkdir('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/')
#     print(item_count, M_dwarves)
#     ##Relevant parameters from TOI catalog
#     period = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'pl_orbper'])[0]
#     epoch = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'pl_tranmid'])[0]
#     stellar_radius = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'st_rad'])[0]
#     T = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == float(M_dwarves), 'st_teff'])[0]
#     # periastron = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == TOI_ID, 'pl_orbper'])[0]
#     # epoch_periastron = np.array(TOI_Catalog.loc[TOI_Catalog['hostname'] == TOI_ID, 'pl_orbtper'])[0]
#     if T == '':
#         T = np.NAN
#     if stellar_radius == '':
#         stellar_radius = np.NAN
#     # if epoch_periastron == '':
#     #     epoch_periastron = np.NAN
#     # if periastron == '':
#     #     periastron = np.NAN
#     if epoch == '':
#         epoch = np.NAN
#     flare_count = 1
    
#     ##Trackable values per star
#     phase_folded_flare_list = []
#     flare_amplitude = []
#     flare_time_scale = []
#     flare_energy = []
#     accepted_flare_index = []
#     flare_phase = []
#     accepted_flare_number = []
#     observation_time = 0
#     possible_flares = 0
#     # periastron_list = []
#     # periastron_epoch_list = []
#     epoch_list = []
    
#     ##Total Trackable Lists for all Data
#     TOI_ID_list = []
#     flare_number = []
#     peak_time = []
#     amplitude = []
#     time_scale = []
#     Teff = []
#     radius = []
#     flare_phase = []
#     TOI_period = []
#     total_flare_energies = []
#     total_flare_phases = []
#     item_count = 0
#     # total_periastron_list = []
#     # total_periastron_epoch_list = []
#     total_epoch_list = []
#     for folders in a:
#         t1 = timer.time()
#         try:
#             b, pdcsap_flux, pdcsap_error = TESS_data_extract('/data/whitsett.n/TESS_Light_Curves/All_TOI/' + M_dwarves + '/' + folders, PDCSAP_ERR=True)
#         except:
#             continue
#         if folders.endswith('a_fast-lc.fits') == True:
#             observation_time += len(time)*(0.33333333)
#             total_observation_time += len(time)*(0.33333333)
#         elif folders.endswith('a_fast-lc.fits') == False:  
#             observation_time += len(time)*2
#             total_observation_time += len(time)*2
#         time, flux = delete_nans(b, pdcsap_flux)
#         detrend_flux = SMA_detrend(time, flux, 80, LS_Iterations=5)
#         t2 = timer.time()
#         tier0_tau.append(t2-t1)
#         flares, lengths = flare_ID(detrend_flux, 3)
#         index = 0
#         possible_flares += len(flares)
#         total_possible_flares += len(flares)

#         for flare_events in flares:
#             if flare_events >= 100 and len(flux) - flare_events > 100:
#                 new_time = time[flare_events-100:flare_events+100]
#                 new_data = flux[flare_events-100:flare_events+100]
#             elif flare_events < 100:
#                 new_time = time[0+flare_events:flare_events+100]
#                 new_data = flux[0+flare_events:flare_events+100]
#             elif len(flux) - flare_events < 100:
#                 new_time = time[flare_events:]
#                 new_data = flux[flare_events:]
#             new_error = pdcsap_error[flare_events-100:flare_events+100]
#             recenter = np.max(new_data[int(len(new_data)/2-10):int(len(new_data)/2+10)])
#             c, d =flare_ID(np.array(new_data), 3)
#             norm_time = time[flare_events]
#             events = np.where(new_data == recenter)[0][0]
#             if period != '':
#                 phase = phase_folder(new_time, period, epoch)
#                 flare_phase_value = phase[events]
#             elif period == '':
#                 phase = []
#                 flare_phase_value = np.NAN
#                 period = np.NAN
                
#             criteria1 = False
#             if recenter > np.mean(new_data)+3*(np.std(new_data)):
#                 criteria1 = True
#             t3 = timer.time()
#             tier1_tau.append(t3-t2)
#             if criteria1 == True and new_data[events+1] > np.mean(new_data)+2*(np.std(new_data)) and len(c) > 0:
#                 new_time = (new_time - new_time[events])*24*60
#                 if lengths[index] >= 25:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+30], alles_data[events:events+30], maxfev=5000)
#                     squares = (alles_data[events:events+30] - exp_decay(new_time[events:events+30], *popt))**2/(np.var(alles_data[events:events+30]))
#                     chi2_cutoff = 20.843
#                 elif lengths[index] >= 15 and lengths[index] < 25:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+20], alles_data[events:events+20], maxfev=5000)
#                     squares = (alles_data[events:events+20] - exp_decay(new_time[events:events+20], *popt))**2/(np.var(alles_data[events:events+20]))
#                     chi2_cutoff = 11.912
#                 elif lengths[index] > 5 and lengths[index] < 15:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+10], alles_data[events:events+10], maxfev=5000)
#                     squares = (alles_data[events:events+10] - exp_decay(new_time[events:events+10], *popt))**2/(np.var(alles_data[events:events+10]))
#                     chi2_cutoff = 3.455
#                 elif lengths[index] <= 5:
#                     alles_data = new_data/np.median(new_data)
#                     error = new_error/np.median(new_data)
#                     popt, pcov = curve_fit(exp_decay, new_time[events:events+7], alles_data[events:events+7], maxfev=5000)
#                     squares = (alles_data[events:events+7] - exp_decay(new_time[events:events+7], *popt))**2/(np.var(alles_data[events:events+7]))
#                     chi2_cutoff = 1.3
#                 chi_squared = np.sum(squares)
#                 if chi_squared < chi2_cutoff and popt[1] > 0 and popt[0] > 0:
#                     half_max = (alles_data[8:15].max()-np.median(alles_data[0:8]))/2
#                     time_scale.append(popt[1])
#                     flare_time_scale.append(popt[1])
#                     amplitude.append(popt[0])
#                     flare_amplitude.append(popt[0])
#                     peak_time.append(norm_time)
                    
#                     TOI_ID_list.append(TOI_ID)
#                     flare_number.append(flare_count)
#                     try:
#                         energy = bolo_flare_energy(popt, T, stellar_radius, pl.planck_integrator(600e-9, 1000e-9, T)/pl.planck_integrator(600e-9, 1000e-9, 9000), t_unit='days')
#                     except:
#                         energy = np.NaN
#                     flare_energy.append(energy)
#                     total_flare_energies.append(energy)
#                     Teff.append(T)
#                     radius.append(stellar_radius)
#                     # try:
#                     #     X = np.column_stack((new_time[events-len(new_time):events+len(new_time)], alles_data[events-len(new_time):events+len(new_time)], error[events-len(new_time):events+len(new_time)]))
#                     # except:
#                     #     X = np.column_stack(([0],[0],[0]))
#                     baseline = st.median(new_data)*(lengths[index])*2
#                     median = st.median(new_data)
#                     flare_phase.append(flare_phase_value)
#                     accepted_flare_index.append(flares[index])
#                     accepted_flare_number.append(flare_count)
#                     total_flare_phases.append(flare_phase_value)
#                     TOI_period.append(period)
#                     # total_periastron_list.append(periastron)
#                     # total_periastron_epoch_list.append(epoch_periastron)
#                     # periastron_list.append(periastron)
#                     # periastron_epoch_list.append(epoch_periastron)
#                     total_epoch_list.append(epoch)
#                     epoch_list.append(epoch)
#                     print(flare_count)
#                     if lengths[index] > 5:
#                         print('Flare ' + str(flare_count) + ' length: ' + str(lengths[index]))
#                     # np.savetxt('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/Flare' + str(flare_count) + '.csv', X, delimiter=',')
#                     flare_count += 1
#                     total_flares += 1
#                     t4 = timer.time()
#                     tier2_tau.append(t4-t3)
                    
                        
#             index += 1
#         ZZ = np.column_stack((np.array(np.array(tier0_tau).mean()), np.array(np.array(tier1_tau).mean()/len(flares)), np.array(np.array(tier2_tau).mean())))
#         with open("/data/whitsett.n/Tier_3/Time_Stats.csv", "a") as f:
#             np.savetxt(f, ZZ, delimiter=",", fmt='%s')
#             f.close()
#         tier0_tau = []
#         tier1_tau = []
#         tier2_tau = []
#     if observation_time == 0:
#         observation_time = -1
#     if possible_flares == 0 :
#         possible_flares = -1
    # Y = np.column_stack((flare_amplitude, flare_time_scale, flare_energy))
    # np.savetxt('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/All_Flares.csv', Y, delimiter=',')
    # f = open('/data/whitsett.n/Flare_Csvs/All_TOI/' + M_dwarves + '/Host_Statistics.txt', 'w')
    # f.write('Flares/Day: ' + str(flare_count*60*24/observation_time) + '\n' + 'Possible Flares: ' + str(possible_flares) + '\n' +
    #         'Accepted Flares: ' + str(flare_count) + '\n' + 'Accepted/Possible Flares: ' + str(flare_count/possible_flares) + '\n' +
    #         'Observation Time (min): ' + str(observation_time))
    # f.close()
    # item_count += 1
    # ZZ = np.column_stack((TOI_ID_list, flare_number, peak_time, amplitude, time_scale, total_flare_energies, Teff, radius, TOI_period))
    # with open("/data/whitsett.n/Tier_3/All_TOI/All_TOI_Flares2.csv", "a") as f:
    #     np.savetxt(f, ZZ, delimiter=",", fmt='%s')
    #     f.close()
    # print(len(TOI_ID_list), len(flare_number), len(peak_time), len(amplitude), len(time_scale), len(total_flare_energies), len(Teff), len(radius), len(TOI_period))
    # print(list_index)
# g = open('/data/whitsett.n/Tier_3/All_TOI/All_statistic_TOI.txt', 'w')
# g.write('Total Flares: ' + str(total_flares) + '\n' + 'Net Observation Time (Days): ' + str(total_observation_time/(60*24)))
# g.close()