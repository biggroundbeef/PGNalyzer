import pdb

from tqdm import tqdm
from numpy import arange, argmin, array, datetime64, timedelta64, ndarray
from pandas import DataFrame
from pgnalyzer.mods.dates import update_str, update_raw
from pgnalyzer.mods.raw_text import read_PGN, apply_dark_correction
from pgnalyzer.spec_cal.solve import fit_spectrum
import matplotlib.pyplot as plt

class Reader():

    def __init__(self):
        print("Reader Initialized")

    def read_PGN_txt(self, input_path) -> dict:
        return read_PGN(input_path)

    def get_dark_corrected_counts(self, input_routines, input_dark_routines, input_param_index) -> dict:
        return apply_dark_correction(input_routines, input_dark_routines, input_param_index)

    def get_spectral_calibrated_channels(
        self,
        in_model_path:str, 
        in_raw_routines: dict,
        in_dora_data: dict,
        in_measurement2obs: int,
        in_a0:float,
        in_b0:float,
        wave_guess: float, 
        chn_guess: int,
        in_fs: int = 4, 
        in_icut: int = 5,
        in_dw: int = 10,
        in_chn_range:tuple = (1, 2052),
        in_spec = "s1",
        in_showPlots:bool = True,
        in_enhance:bool = False
    ):
        a, b = in_a0, in_b0
        if in_enhance:
            for i in tqdm(range(10), desc = "Repeating log-linear model to find best a and b: "): 
                _, a, b, corr = fit_spectrum(
                    model_path = in_model_path, 
                    raw_routines = in_raw_routines, 
                    dora_data = in_dora_data, 
                    measurement2obs = in_measurement2obs, 
                    fs = in_fs, 
                    icut = in_icut, 
                    dw = in_dw,
                    chn_range = in_chn_range,
                    spec = in_spec,
                    showPlots = False,
                    a0 = a,
                    b0 = b,
                    wref = wave_guess,
                    cref0 = chn_guess, 
                )

        return fit_spectrum(
            model_path = in_model_path, 
            raw_routines = in_raw_routines, 
            dora_data = in_dora_data, 
            measurement2obs = in_measurement2obs, 
            fs = in_fs, 
            icut = in_icut, 
            dw = in_dw,
            chn_range = in_chn_range,
            spec = in_spec,
            showPlots = in_showPlots,
            a0 = a,
            b0 = b,
            wref = wave_guess,
            cref0 = chn_guess, 
        )

    def plot_corrected_data(self, in_chn, in_cc, raw_data = None, param_idx = None, fs = 5):
        for i, counts in enumerate(in_cc):
            fig, axes = plt.subplots(1, 1, figsize = (2 * fs, fs))
            axes2 = axes.twinx()
            if raw_data is not None and param_idx is not None:
                raw_plt_data = raw_data[i][param_idx["Pixel 1"]:param_idx["Unc 1"]].astype(float)
                axes.plot(
                    arange(start = 1, stop = len(raw_plt_data)+1, step = 1),
                    raw_plt_data,
                    color = "#75F94D",
                    label = "Raw Data",
                    # alpha = 0.4
                )
                axes.plot(
                    in_chn,
                    raw_plt_data,
                    color = "#377D22",
                    label = "Raw Data - Spec Cal Adjusted",
                    # alpha = 0.4
                )
            axes2.plot(
                in_chn,
                counts,
                color = "red",
                label = "CC - Spec Cal Adjusted"
            )
            axes.set_title("Raw BC vs Calibrated CC")
            axes.set_xlabel("Channel // Wavelength (nm)")
            axes.set_ylabel("Raw Counts")
            axes2.set_ylabel(r"$\frac{W}{\mu m^2\quadm^2}$")
            axes.set_box_aspect(0.5)
            axes.legend(loc = "upper left", bbox_to_anchor = (0, -0.2))
            axes2.legend(loc = "upper right", bbox_to_anchor = (1, -0.2))
            fig.tight_layout()
            plt.show()

def quick_test():
    """
    Default first guesses:
    a0 = 6.3
    b0 = 0.0013
    """
    ron = Reader()
    
    in_PGN_data_test = r"C:\GSFC\2025\05 May\Pandora2s1_GreenbeltMD_20250505_L0.txt"
    in_model_path = r"C:\GSFC\model_01292025_Reptranfiner0.1.txt"
    routines, dora_data = ron.read_PGN_txt(in_PGN_data_test)
    
    corrected_counts = ron.get_dark_corrected_counts(routines["S"]["SQ"], routines["DarkS"]["SQ"], dora_data["P_IDX"])
    
    spectral_calibrated_channels, a, b, corr = ron.get_spectral_calibrated_channels(
        in_model_path = in_model_path, 
        in_raw_routines = routines["S"]["SQ"],
        in_dora_data = dora_data,
        in_measurement2obs = 3,
        in_showPlots = False,
        in_a0 = 4.888887602820383,
        in_b0 = 0.0012425618959219836,
        in_enhance = True, 
        wave_guess = 486,
        chn_guess = 1565
    )
    # print(f"The final fitting results are:\na0 = {a}\nb0 = {b}\nCorr = {corr}")

    ron.plot_corrected_data(spectral_calibrated_channels, corrected_counts["CC"], routines["S"]["SQ"], dora_data["P_IDX"])
    
    print("Success!")
    pdb.set_trace()

if __name__ == "__main__":
    quick_test()