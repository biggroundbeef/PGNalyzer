import pdb

from pgnalyzer.io.model import get_model
from pgnalyzer.mods.data import smooth
from pgnalyzer.mods.dates import update_str, update_raw
from pgnalyzer.astro.sza import get_solar_zenith_angle
from numpy import abs, arange, argmax, array, copy, correlate, datetime64, exp, timedelta64, interp, linspace, log, polyfit, polyval, min, max, ndarray, union1d, where
from tqdm import tqdm
from pandas import DataFrame
import matplotlib.pyplot as plt

# Eddy IDL -> python
def get_channel(
    wl_model_buff, 
    irr_model_buff, 
    wref_buff, 
    dw, 
    chn, 
    rad_buff, 
    a_buff, 
    b_buff, 
    cref0_buff
) -> int:
    cmax_buff = 0
    nn = 100 # Searching channel range +/- nn

    # Create search range in wavelengths
    idx_1_chn = where(abs(wl_model_buff - wref_buff) < dw)
    wl_m_buff = wl_model_buff[idx_1_chn]

    # Using min/max normalization on irradiance model
    want_model_buff = (
        (irr_model_buff[idx_1_chn] - min(irr_model_buff[idx_1_chn]))
        / (max(irr_model_buff[idx_1_chn]) - min(irr_model_buff[idx_1_chn]))
    )
    wavelength_range = arange(
        start = -nn, 
        stop = nn, 
        step = 1
    )
    # Initial fit to find channel for cref_buff
    for k_buff in wavelength_range:
        cref_buff = cref0_buff + 0.5*k_buff
        wl_buff = wref_buff + log((a_buff+b_buff*chn) / (a_buff+b_buff*cref_buff)) / b_buff

        idx_2_chn = where(abs(wl_buff - wref_buff) < dw)
        # Repeating min/max normalization for rad
        try:
            want_rad_buff = (
                (rad_buff[idx_2_chn] - min(rad_buff[idx_2_chn])) * 1
                / (max(rad_buff[idx_2_chn]) - min(rad_buff[idx_2_chn]))        
            )
            want_buff = interp(wl_m_buff, wl_buff[idx_2_chn], want_rad_buff)
            corr_buff = correlate(want_model_buff, want_buff, mode = "valid")
            
            if corr_buff > cmax_buff:
                cmax_buff = corr_buff
                cfit_buff = cref_buff
        
        except:
            print(f"ERROR: Iterating through wavelengths --> k: {k_buff}")
    try: 
        if not cfit_buff: print(cfit_buff)            
    except: 
        cfit_buff = 1000
        cmax_buff = 0.001
    return cfit_buff, cmax_buff

def fit_spectrum(
    model_path:str, 
    raw_routines: dict,
    dora_data: dict,
    measurement2obs: int,
    a0:float,
    b0:float,
    wref:float,
    cref0:float,
    rmin: float = 8e4,
    fs: int = 4, 
    icut: int = 5,
    dw: int = 10,
    chn_range:tuple = (1, 2052),
    w0:ndarray = array([
        822.7, 759.4, 722.2, 686.7, 656.3, 627.7, 589.3, 527.0, 
        518.4, 495.8, 486.1, 438.6, 430.8, 410.2, 396.8, 
        382.0, 358.1, 336.1, 
        # 302.1, 299.4
    ]),
    spec = "s1",
    showPlots:bool = True
):
    """
    This function replicates Dr. Wu's work

    Default w0:
    w0 = np.array([
    822.7       , 759.4         , 722.2         , 686.7         , 656.3     , 627.7     , 589.3     , 527.0 , 
    O2-Z        , O2-A          , H2O           , O2-B          , H-a       , O2-a      , Na        , Fe    ,    
    518.4       , 495.8         , 486.1         , 438.6         , 430.8     , 410.2     , 396.8     ,        
    Mg(b1)      , Fe            , H-b           , Fe            , Ca/Fe     , H-d       , Ca+       ,       , 
	382.0       , 358.1         , 336.1         ,                                                       
	Fe          , Fe            , TI+           ,               ,           ,           ,           ,       ,   
    302.1       , 299.4                                                                               
	Fe          , Ni            ,               ,               ,           ,           ,           ,       ,
])

    """
    chn = arange(start = chn_range[0], stop = chn_range[1]+1, step = 1)
    nch = len(w0)
    nch1 = 7
    nch2 = nch

    wavelength_model, irradiance_model = get_model(model_path)
    irradiance_model_smooth = smooth(irradiance_model, icut)

    plt.close("all")
    # Plot 1
    # if showPlots:
    #     fig, axes = plt.subplots(1, 1, figsize = (2.5 * fs, fs))
    #     axes.plot(
    #         wavelength_model,
    #         irradiance_model,
    #         color = "blue",
    #         linewidth = 0.5,
    #         label = "Model"
    #     )
    #     axes.plot(
    #         wavelength_model,
    #         irradiance_model_smooth,
    #         color = "orange",
    #         linewidth = 0.5,
    #         label = "Smoothed Model"
    #     )  
    #     axes.set_title("Model")
    #     axes.set_xlabel("Wavelength (nm)")
    #     axes.set_ylabel("Intensity")
    #     axes.legend(loc = "upper right", bbox_to_anchor = (1, 1))
    #     fig.suptitle("")        
    #     fig.tight_layout()
    #     plt.show()

    param_idx = dora_data["P_IDX"]
    sun_dates = array(list(map(lambda date: update_str(update_raw(date)), raw_routines[:, param_idx["UT date and time for"]])))
    sun_data = raw_routines[:, param_idx["Pixel 1"]:param_idx["Unc 1"]].astype(float)
    lat, long = float(dora_data["LATI"]), float(dora_data["LONG"])

    # Plot 2
    # if showPlots:
    #     test_data_1 = sun_data[10]
    #     test_data_2 = sun_data[11]
    #     test_rad = test_data_1 - test_data_2
    #     fig, axes = plt.subplots(1, 1, figsize = (2.5 * fs, fs))
        # axes.plot(
        #     arange(start = 1, stop = len(test_data_1) + 1),
        #     test_rad,
        #     color = "blue",
        #     linewidth = 0.5,
        #     label = "Spec Data"
        # )
        # axes.legend(loc = "upper left", bbox_to_anchor = (1,1))
        # axes.set_title(f"Spec Data: SZA = {get_solar_zenith_angle(sun_dates[10], lat, long):.3f}")
        # axes.set_xticks(linspace(axes.get_xlim()[0], axes.get_xlim()[1], 50, dtype = int))
        # axes.set_yticks(union1d(linspace(axes.get_ylim()[0], axes.get_ylim()[1], 5, dtype = int), [0]))
        # xticks, yticks = axes.get_xticks(), axes.get_yticks()
        # # axes.xaxis.set_ticklabels(axes.get_xticklabels(), rotation = 90)
        # for lab in axes.get_xticklabels():
        #     lab.set_rotation(90)
        # axes.set_xlabel("Channel")
        # axes.set_ylabel("Rad")
        # fig.suptitle("")
        # fig.tight_layout()
        # plt.show()

    # if spec == "s1":
        # Following variables are for the first guess for s1
        # a0      = 6.3
        # b0      = 0.0013
        # a0 = 5.67
        # b0 = 0.00117
        # corr = [2.8358317e+19]
        # a0 = 5.103
        # b0 = 0.0010530000000000001
        # corr = [3.12896945e+19]
        # a0 = 5.05197
        # b0 = 0.0010951200000000002
        # corr = [3.18212983e+19]
        # wref    = 430.8     # nm
        # cref0   = 1150      # first guess of channel number of this line
        # rmin    = 8e4

    for i, timestamp in enumerate(sun_dates):
        if i != len(sun_dates) - 1:
            # if i in range(len(routes2use["dates"])):
            if i == measurement2obs:
                sza = get_solar_zenith_angle(datetime64(timestamp), lat, long)

                rad = sun_data[i] - sun_data[i + 1]  
                if max(rad) > rmin:
                    cfit, _ = get_channel(
                        wl_model_buff    = wavelength_model, 
                        irr_model_buff   = irradiance_model_smooth, 
                        wref_buff        = wref, 
                        dw               = dw, 
                        chn              = chn, 
                        rad_buff         = rad, 
                        a_buff           = a0, 
                        b_buff           = b0, 
                        cref0_buff       = cref0
                    )
                    cref = cfit

                    idx = where(abs(wavelength_model - wref < dw))
                    want_model = (
                        (irradiance_model[idx] - min(irradiance_model[idx]))
                        / (max(irradiance_model[idx]) - min(irradiance_model[idx]))
                    ) 
                    smooth_model = (
                        (irradiance_model_smooth[idx] - min(irradiance_model_smooth[idx]))
                        / (max(irradiance_model_smooth[idx]) - min(irradiance_model_smooth[idx]))
                    ) 
                    # Plot 3
                    # if showPlots:
                    #     fig, axes = plt.subplots(1, 1, figsize = (2.5 * fs, fs))
                    #     axes.plot(
                    #         wavelength_model[idx], 
                    #         want_model, 
                    #         linewidth = 0.5, 
                    #         color = "blue",
                    #         label = "Want Model"
                    #     )
                    #     axes.plot(
                    #         wavelength_model[idx],
                    #         smooth_model,
                    #         linewidth = 0.5,
                    #         color = "orange",
                    #         label = "Smooth Model"
                    #     )
                
                    afit, bfit = a0, b0
                    wl = wref + log((afit + bfit * chn) / (afit + bfit * cref)) / bfit
                    idx = where(abs(wl - wref < dw))
                    want_rad = (
                        (rad[idx] - min(rad[idx])) * 1
                        / (max(rad[idx]) - min(rad[idx]))
                    )
                    # if showPlots:
                    #     axes.plot(
                    #         wl[idx],
                    #         want_rad,
                    #         color = "red",
                    #         linewidth = 0.5,
                    #         label = "Want Rad"
                    #     )
                    #     axes.legend(loc = "upper left", bbox_to_anchor = (0, -0.1))
                    #     axes.set_xlabel("Wavelength")
                    #     axes.set_ylabel("Rad")
                    #     fig.suptitle("")   
                    #     fig.tight_layout()
                    #     plt.show()

                    c0 = copy(w0) * 0
                    a, b = a0, b0
                    for k in arange(nch1, nch2 - 1):
                        cref0 = (
                            ((a+b*cref)  * exp((w0[k] - wref) * b) - a)
                            / b
                        )
                        cfit, _ = get_channel(wavelength_model, irradiance_model_smooth, w0[k], dw, chn, rad, a, b, cref0)
                        c0[k] = cfit
                        idx = where(abs(wavelength_model - w0[k]) < dw)
                        want_model = (
                            (irradiance_model[idx] - min(irradiance_model[idx]))
                            / (max(irradiance_model[idx]) - min(irradiance_model[idx]))
                        )
                        smooth_model = (
                            (irradiance_model_smooth[idx] - min(irradiance_model_smooth[idx]) )
                            / (max(irradiance_model_smooth[idx]) - min(irradiance_model_smooth[idx])) 
                        )
                        # Plot 4
                        # if showPlots:
                        #     fig, axes = plt.subplots(1, 1, figsize = (2.5 * fs, fs))
                        #     axes.plot(
                        #         wavelength_model[idx],
                        #         want_model,
                        #         color = "blue",
                        #         linewidth = 0.5,
                        #         label = "Model"
                        #     )
                        #     axes.plot(
                        #         wavelength_model[idx],
                        #         smooth_model,
                        #         color = "orange",
                        #         linewidth = 0.5,
                        #         label = "Smooth Model"
                        #     )
                        wl = w0[k] + log((a + b*chn) / (a + b*cfit)) / b
                        idx = where(abs(wl - w0[k]) < dw)
                        want_rad = (
                            (rad[idx] - min(rad[idx]))
                            / (max(rad[idx]) - min(rad[idx]))
                        )
                        # if showPlots:
                        #     axes.plot(
                        #         wl[idx],
                        #         want_rad,
                        #         color = "red",
                        #         linewidth = 0.5,
                        #         label = "Want Rad"
                        #     )
                        #     axes.set_xlabel("Wavelength (nm)")
                        #     axes.set_ylabel("Intensity")
                        #     fig.suptitle("")   
                        #     fig.tight_layout()
                        #     plt.show()

                    mm = 10
                    tmax = 0
                    all_res = []
                    for iter_i in tqdm(range(-mm, mm), desc = "Iterating through log-linear model to compute a and b"):
                        for iter_j in range(-mm, mm):
                            a = a0 *(1 + 0.01*iter_i)
                            b = b0 *(1 + 0.01*iter_j)
                            amax = 1
                            for iter_k in range(nch1, nch2 - 1):
                                cfit, cmax = get_channel(wavelength_model, irradiance_model_smooth, w0[iter_k], dw, chn, rad, a, b, c0[iter_k])
                                amax = amax * cmax
                            if amax > tmax:
                                # print(w0[iter_k], cfit, tmax, " -> ", amax)
                                tmax = amax
                                afit = a    
                                bfit = b
                            
                                # print(iter_i, iter_j, afit, bfit, tmax)
                                all_res.append({"i" : iter_i, "j" : iter_j, "a" : afit, "b" : bfit, "corr" : tmax})
                    all_res = DataFrame(all_res)
                    best_idx = argmax(all_res["corr"])
                    print(f"New best guess:\na0 = {all_res["a"][best_idx]}\nb0 = {all_res["b"][best_idx]}\ncorr = {all_res["corr"][best_idx]}")
                    # a=afit & b=bfit
                    # wl=wref+alog((a+b*chn)/(a+b*cref))/b	; log-linear model

                    for nfit in range(1, 2+1):
                        # Plot 5
                        nonzero_idx = where(c0 != 0)
                        num_nonzero_idx = len(nonzero_idx)
                        if showPlots:
                            fig, axes = plt.subplots(1, 1, figsize = (fs, fs))
                            axes.plot(
                                c0[nonzero_idx],
                                w0[nonzero_idx],
                                label = "Poly Fit Model",
                                color = "blue",
                                ls = "",
                                marker = "D",
                                markersize = 1
                            )
                            axes.set_title("Polyfit Model")
                            axes.legend(loc = "upper left", bbox_to_anchor = (0, -.15))
                            axes.set_xlabel("Channel")
                            axes.set_ylabel("Wavelength")
                            axes.set_box_aspect(1)

                        x = c0[nonzero_idx]
                        # x = np.zeros((nfit, num_nonzero_idx))
                        # for it in range(0, nfit-1):
                            # x[it, :] = c0[nonzero_idx]**(it+1)

                        cc = polyfit(x, w0[nonzero_idx], nfit)
                        cc0 = cc[-1]
                        yfit = polyval(cc, x)
                        cc = cc[:-1]

                        # if showPlots:
                        #     axes.plot(
                        #         x,
                        #         yfit,
                        #         color = "red",
                        #         linewidth = 0.5,
                        #         alpha = 0.5
                        #     )

                        wl = cc0
                        for iter_i in range(0, nfit):
                            wl = wl + cc[iter_i] * chn**(iter_i + 1)
                        wl_poly = wl 
                        if showPlots:
                            axes.plot(
                                chn,
                                wl,
                                color = "red",
                                alpha = 0.5,
                                ls = '--'
                            )

                            fig.tight_layout()
                            plt.show()

                        want_model = irradiance_model / max(irradiance_model)
                        fig, axes = plt.subplots(1, 1, figsize = (2.5 * fs, fs))
                        smooth_model = irradiance_model_smooth / max(irradiance_model_smooth)
                        want_rad = rad * 1 / max(rad)
                        if showPlots:
                            axes.plot(
                                wavelength_model,
                                want_model,
                                color = "blue",
                                label = "Model",
                                linewidth = 0.5
                            )
                            axes.plot(
                                wavelength_model,
                                smooth_model,
                                color = "orange",
                                label = "Smooth Model",
                                linewidth = 0.5
                            )
                            axes.plot(
                                wl_poly,
                                want_rad,
                                color = "red",
                                linewidth = 0.5,
                                label = "Spec Cal Data"
                            )
                            ylim = axes.get_ylim()
                            axes.vlines(
                                wref,
                                -1e6, 
                                1e10,
                                color = "black",
                                label = f"Wavelength Guess - {wref}nm",
                                ls = "--",
                                alpha = 0.5
                            )
                            axes.legend(loc = "upper left", bbox_to_anchor = (0, -0.1))
                            axes.set_ylim(ylim)
                            axes.set_xlabel("Wavelength (nm)")
                            axes.set_ylabel(r"Intensity")
                            axes.set_title("Model VS Spectral Calibrated Data")
                            fig.tight_layout()
                            plt.show()

                        return wl_poly, afit, bfit, tmax
                    