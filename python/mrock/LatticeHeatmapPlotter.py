import numpy as np
import pandas as pd
import string
from scipy.signal import find_peaks

from . import path_appender as __ap
__ap.append()
import continued_fraction as cf
import spectral_peak_analyzer as spa
from legend import *
from matplotlib import ticker

from .create_figure import create_normal_figure, create_large_figure
from make_panels_touch import make_panels_touch
from label_axes import label_axes

#MAX_EXP = 2.1
#CUT_OFF_EXP = -1.5

# Settings the importer
G_MAX_LOAD = 3

__BEGIN_OFFSET__ = 2e-3
__RANGE__ = 2e-3
__SECOND_BEGIN__ = 1e-9
__SECOND_RANGE__ = 1e-8
__PEAK_PROMINCE__ = 1.

__INITIAL_IMAG__ = 1e-5j
__FIT_COMPLEX_SHIFT__ = 1e-8j
__CONTINUUM_CUT_SHIFT__ = 1e-5

__sigma__ = 0.0005

def gaussian_bell(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))

def derivative_gaussian_bell(x, mu, sigma):
    return -((x - mu) / sigma**2) * gaussian_bell(x, mu, sigma)

def is_phase_peak(peak):
    return abs(peak) < 1.5e-3

def extract_model_settings(ds):
    return f"dos={ds['dos_name']}  omega_D={ds['omega_D']}  g={ds['g']}  U={ds['U']}  E_F={ds['E_F']}"

class BaseCFIgnore:
    def __init__(self, first=70, last=400):
        self.first = int(first)
        self.last  = int(last)
    def get_first(self, g):
        return self.first
    def get_last(self, g):
        return self.last

class GradCFIgnore(BaseCFIgnore):
    def __init__(self, first, last, first_grad, last_grad):
        super().__init__(first, last)
        self.first_grad = first_grad
        self.last_grad  = last_grad
    
    def get_first(self, g):
        return int(self.first + g * self.first_grad)
    def get_last(self, g):
        return int(self.last + g * self.last_grad)

class BracketCFIgnore(BaseCFIgnore):
    def __init__(self, firsts, lasts, upper_bounds):
        super().__init__(firsts[0], lasts[-1])
        self.firsts   = np.asarray(firsts)
        self.lasts    = np.asarray(lasts)
        self.upper_bounds = np.asarray(upper_bounds)
    
    def get_bracket_index(self, g):
        mask = self.upper_bounds > g
        bracket_index = np.argmax(mask)
        if not mask[bracket_index]:
            bracket_index -= 1
        return bracket_index
    
    def get_first(self, g):
        return int(self.firsts[self.get_bracket_index(g)])
    def get_last(self, g):
        return int(self.lasts[self.get_bracket_index(g)])

def merge_data_dicts(dicts, distinguishers, same_omega=True, same_g=False, has_higgs=True, has_phase=True):
    """Merges the data_dicts created by instances of Heatmaplotter.
    dicts: the data_dicts
    distinguishers: string appended to distinguish the data, e.g., g -> g_bcc."""
    
    merged = {}
    if same_omega:
        merged["omega"] = dicts[0]["omega"]
    if same_g:
        merged["g"] = dicts[0]["g"]
    for data, distinguish in zip(dicts, distinguishers):
        merged[f"continuum_edge_{distinguish}"] = data["continuum_edge"]
        if has_higgs:
            merged[f"spectral_functions_higgs_{distinguish}"] = data["spectral_functions_higgs"]
        if has_phase:
            merged[f"spectral_functions_phase_{distinguish}"] = data["spectral_functions_phase"]
        
        if not same_omega:
            merged[f"omega_{distinguish}"] = data["omega"]
        if not same_g:
            merged[f"g_{distinguish}"] = data["g"]
    
    return merged

class HeatmapPlotter:
    def __init__(self, data_frame_param, parameter_name, xlabel, zlabel=r'$\mathcal{A}(\omega) / W^{-1}$', xscale="linear", yscale="linear",
                 energy_range=(1e-10, 2.5), scale_energy_by_gaps=False, cf_ignore=BaseCFIgnore()):
        self.data_frame = data_frame_param.sort_values(parameter_name).reset_index(drop=True)
        
        self.y = np.linspace(energy_range[0], energy_range[1], 1000) # meV
        self.x = (self.data_frame[parameter_name]).to_numpy()
        self.scale_energy_by_gaps = scale_energy_by_gaps
        self.resolvents = [cf.ContinuedFraction(pd_row, 
                                                messages=False, 
                                                ignore_first=cf_ignore.get_first(pd_row['g']), 
                                                ignore_last=cf_ignore.get_last(pd_row['g'])) 
                           for index, pd_row in self.data_frame.iterrows()]
        self.max_gaps   = np.array([2 * gap for gap in self.data_frame["Delta_max"]]) # meV
        self.true_gaps  = np.array([float(t_gap[0]) for t_gap in self.data_frame["continuum_boundaries"]]) # meV
        self.N_data = len(self.max_gaps)

        self.xlabel = xlabel
        self.zlabel = zlabel
        self.xscale = xscale
        self.yscale = yscale
        
        self.data_dict = {
            "g": self.x,
            "omega": self.y
        }

    def __scale_if__(self, x, i):
        if self.scale_energy_by_gaps:
            return self.max_gaps[i] * x
        return x

    def identify_modes(self, spectral, pos):
        if self.max_gaps[pos] == 0:
            return np.array([])
        indizes = find_peaks(spectral, prominence=8 * self.true_gaps[pos] * __PEAK_PROMINCE__)[0]
        positions = np.array([self.__scale_if__(self.y[i], pos) for i in indizes])
        positions = positions[positions < self.true_gaps[pos] - __CONTINUUM_CUT_SHIFT__]
        positions[is_phase_peak(positions)] = 0
        return positions

    def __get_error__(self, key, i):
        __std_g__ = 0.05
        base = self.data_frame["g"].iloc[i]
        
        filter_lower = self.data_frame["g"][(self.data_frame["g"] < base - __std_g__)]
        filter_upper = self.data_frame["g"][(self.data_frame["g"] > base + __std_g__)]
        
        lower_err_idx = filter_lower.idxmax() if len(filter_lower) > 0 else -1
        upper_err_idx = filter_upper.idxmin() if len(filter_upper) > 0 else 2**31 - 1
        
        if key == "true_gap":
            return [ np.abs(self.true_gaps[i] - self.true_gaps[lower_err_idx]) if lower_err_idx >= 0                  else None,
                     np.abs(self.true_gaps[i] - self.true_gaps[upper_err_idx]) if upper_err_idx < len(self.true_gaps) else None ]
            
        lower_err = np.abs(self.data_frame[key].iloc[i] - self.data_frame[key].iloc[lower_err_idx]) if lower_err_idx >= 0                  else None
        upper_err = np.abs(self.data_frame[key].iloc[i] - self.data_frame[key].iloc[upper_err_idx]) if upper_err_idx < len(self.true_gaps) else None
        return [lower_err, upper_err]
        
    def __decide_if_to_reverse__(self, peak_idx, peaks, gap):
        if len(peaks) > peak_idx and peaks[peak_idx] - __BEGIN_OFFSET__ - __RANGE__ <= 0:
            return False
        if len(peaks) > peak_idx + 1:
            return abs(peaks[peak_idx + 1] - peaks[peak_idx]) < np.sqrt(__RANGE__)
        else:
            if len(peaks) > peak_idx:
                return abs(peaks[peak_idx] - gap) < __RANGE__
            else:
                return False
    
    def fit_goldstone_peak(self, _real, _imag, i):
        __result = spa.analyze_peak(_real, _imag, 
                                    peak_position         = 0, 
                                    range                 = __RANGE__,
                                    begin_offset          = __BEGIN_OFFSET__,
                                    reversed              = False,
                                    lower_continuum_edge  = self.true_gaps[i],
                                    peak_pos_range        = self.y[20] - self.y[0],
                                    improve_peak_position = False)

        current_range = __RANGE__
        current_offset = __BEGIN_OFFSET__
        
        __best_fit = __result

        def deviation(_current):
            return abs(_current.slope.nominal_value + 2)
        def break_condition():
            return abs(__result.slope.nominal_value + 2) > 0.2
        
        while break_condition() and (current_range >= __SECOND_RANGE__):
            current_offset = __BEGIN_OFFSET__
            while break_condition() and (current_offset >= __SECOND_BEGIN__):
                # Retry the fit with changed parameters
                __result = spa.analyze_peak(_real, _imag, 
                                            peak_position         = __result.position, 
                                            range                 = 1e2 * current_range,
                                            begin_offset          = __result.position + 1e2 * current_offset,
                                            reversed              = False,
                                            lower_continuum_edge  = self.true_gaps[i],
                                            improve_peak_position = False)
                if deviation(__result) < deviation(__best_fit):
                    __best_fit = __result
                current_offset *= 0.5
            current_range *= 0.5
            
        __result = __best_fit
        
        if abs(__result.slope.nominal_value + 2) > 0.33:
            print("WARNING in Phase! Expected slope of '-2' does not match fitted slope!")
            print(__result)
            print(extract_model_settings(self.data_frame.iloc[i]), "\n")
        
        return __result
        
    # Legacy: using the Kramers-Kronig relations to obtain the weights
    # New version programmed in cpp - works well for delta-peaks
    # KK version is needed for Goldstone boson (delta derivative peaks)
    def try_to_fit(self, _real, _imag, _intial_positions, i, higgs):
        __result = [ spa.analyze_peak(  _real, _imag, 
                                        peak_position         = peak_position, 
                                        range                 = __RANGE__,
                                        begin_offset          = __BEGIN_OFFSET__ if not is_phase_peak(peak_position) else peak_position + __BEGIN_OFFSET__,
                                        reversed              = self.__decide_if_to_reverse__(__i, _intial_positions, self.true_gaps[i]),
                                        lower_continuum_edge  = self.true_gaps[i],
                                        peak_pos_range        = self.y[20] - self.y[0],
                                        improve_peak_position = not is_phase_peak(peak_position))
                                for __i, peak_position in enumerate(_intial_positions) ]
        
        for j in range(len(__result)):
            current_range = __RANGE__
            current_offset = __BEGIN_OFFSET__
            
            __best_fit = __result[j]
            
            def deviation(_current):
                if higgs or not is_phase_peak(_current.position):
                    return abs(_current.slope.nominal_value + 1)
                else:
                    return abs(_current.slope.nominal_value + 2)
            
            def break_condition():
                if higgs or not is_phase_peak(__result[j].position):
                    return abs(__result[j].slope.nominal_value + 1) > 0.005
                else:
                    return abs(__result[j].slope.nominal_value + 2) > 0.2
            
            while break_condition() and (current_range >= __SECOND_RANGE__):
                current_offset = __BEGIN_OFFSET__
                while break_condition() and (current_offset >= __SECOND_BEGIN__):
                    # Retry the fit with changed parameters
                    __result[j] = spa.analyze_peak(_real, _imag, 
                                                peak_position         = __result[j].position, 
                                                range                 = current_range if not is_phase_peak(__result[j].position) else 1e2 * current_range,
                                                begin_offset          = current_offset if not is_phase_peak(__result[j].position) else __result[j].position + 1e2 * current_offset,
                                                reversed              = self.__decide_if_to_reverse__(j, _intial_positions, self.true_gaps[i]),
                                                lower_continuum_edge  = self.true_gaps[i],
                                                improve_peak_position = False)
                    if deviation(__result[j]) < deviation(__best_fit):
                        __best_fit = __result[j]
                    current_offset *= 0.5
                current_range *= 0.5
                
            __result[j] = __best_fit
            
            if not higgs and is_phase_peak(__result[j].position):
                if abs(__result[j].slope.nominal_value + 2) > 0.33:
                    print("WARNING in Phase! Expected slope of '-2' does not match fitted slope!")
                    print(__result[j])
                    print(extract_model_settings(self.data_frame.iloc[i]), "\nReversed=", self.__decide_if_to_reverse__(j, _intial_positions, self.true_gaps[i]), "\n")
            elif abs(__result[j].slope.nominal_value + 1) > 0.01:
                print(f"WARNING in {'Higgs' if higgs else 'Phase'}! Expected slope of '-1' does not match fitted slope!")
                print(__result[j])
                print(extract_model_settings(self.data_frame.iloc[i]), "\nReversed=", self.__decide_if_to_reverse__(j, _intial_positions, self.true_gaps[i]), "\n")
        
        return __result

    def compute_phase_peaks(self):
        phase_cpp_results = [ resolvent.classify_bound_states("phase_SC", 
                                                              weight_domega=1e-8, 
                                                              is_phase_peak=is_phase_peak) 
                                for resolvent in self.resolvents ]
        
        __phase_peak_positions = [ [data[0] for data in cpp_result] for cpp_result in phase_cpp_results ]
        __phase_peak_weights   = [ [data[1] for data in cpp_result] for cpp_result in phase_cpp_results ]

        for i, res in enumerate(self.resolvents):
            if self.max_gaps[i] < __RANGE__ + __BEGIN_OFFSET__:
                continue
            __phase_peak_positions[i].insert(0, 0.0)
            __phase_peak_weights[i].insert(0, 0.0)
            # Real part should not have an imaginary shift -> yields better fits
            # because the analytical form is then 1/x rather than x/(x^2 + delta^2)
            __phase_real = lambda x: res.continued_fraction(x, "phase_SC").real
            # Imaginary part needs an imaginary shift to resolve the delta peaks
            __phase_imag = lambda x: res.continued_fraction(x + __FIT_COMPLEX_SHIFT__, "phase_SC").imag
            
            __phase_result = self.fit_goldstone_peak(__phase_real, __phase_imag, i)
             
            __phase_peak_positions[i][0] = __phase_result.position
            __phase_peak_weights[i][0]   = __phase_result.weight
        
        return (__phase_peak_positions, __phase_peak_weights)

    def compute_higgs_peaks(self):
        higgs_cpp_results = [ resolvent.classify_bound_states("amplitude_SC", weight_domega=1e-8) 
                                for resolvent in self.resolvents ]

        __higgs_peak_positions = [ [data[0] for data in cpp_result] for cpp_result in higgs_cpp_results ]
        __higgs_peak_weights   = [ [data[1] for data in cpp_result] for cpp_result in higgs_cpp_results ]
        
        return (__higgs_peak_positions, __higgs_peak_weights)

    def __remove_data_below_continuum__(self, spectral_functions):
        if not self.scale_energy_by_gaps:
            for i in range(self.N_data):
                spectral_functions[:, i][self.y < self.true_gaps[i] - __CONTINUUM_CUT_SHIFT__] = 0
        else:
            for i in range(self.N_data):
                spectral_functions[:, i][self.y * self.max_gaps[i] < self.true_gaps[i] - __CONTINUUM_CUT_SHIFT__] = 0

    def plot_one(self, ax, cmap, which="amplitude_SC", continuum_color="cyan", min_exp=-1.5, max_exp=2.1):
        spectral_functions = np.array([res.spectral_density(self.__scale_if__(self.y, __i) + __INITIAL_IMAG__, which) for __i, res in enumerate(self.resolvents)]).transpose()

        if which == "amplitude_SC":
            _peak_positions, _peak_weights = self.compute_higgs_peaks()
            self.HiggsModes = pd.DataFrame([ {
                "resolvent_type": "Higgs",
                "energies": _peak_positions[i],
                "weights": _peak_weights[i],
                "Delta_max": self.data_frame["Delta_max"].iloc[i],
                "true_gap": self.true_gaps[i],
                "g": self.data_frame["g"].iloc[i],
                "error_g": self.__get_error__("g", i),
                "error_Delta_max": self.__get_error__("Delta_max", i),
                "error_true_gap": self.__get_error__("true_gap", i),
                "omega_D": self.data_frame["omega_D"].iloc[i],
                "E_F": self.data_frame["E_F"].iloc[i],
                "U": self.data_frame["U"].iloc[i],
                "dos_name" : self.data_frame["dos_name"]
            } for i in range(self.N_data) ])
        else:
            _peak_positions, _peak_weights = self.compute_phase_peaks()
            self.PhaseModes = pd.DataFrame([ {
                "resolvent_type": "Phase",
                "energies": _peak_positions[i],
                "weights": _peak_weights[i],
                "Delta_max": self.data_frame["Delta_max"].iloc[i],
                "true_gap": self.true_gaps[i],
                "g": self.data_frame["g"].iloc[i],
                "error_g": self.__get_error__("g", i),
                "error_Delta_max": self.__get_error__("Delta_max", i),
                "error_true_gap": self.__get_error__("true_gap", i),
                "omega_D": self.data_frame["omega_D"].iloc[i],
                "E_F": self.data_frame["E_F"].iloc[i],
                "U": self.data_frame["U"].iloc[i],
                "dos_name" : self.data_frame["dos_name"]
            } for i in range(self.N_data)])

        self.__remove_data_below_continuum__(spectral_functions)
        ## Note, that the phase peak at omega=0 is the derivative of a delta peak
        ## while the other peaks below the continuum are proper delta peaks
        for i in range(self.N_data):
            for peak_position, weight in zip(_peak_positions[i], _peak_weights[i]):
                if is_phase_peak(peak_position):
                    summand = -weight * derivative_gaussian_bell(self.__scale_if__(self.y, i), 0, __sigma__)
                else:
                    summand =  weight * gaussian_bell(self.__scale_if__(self.y, i), peak_position, __sigma__)
                spectral_functions[:, i] += summand
             
        levels = np.power(10, np.linspace(min_exp, max_exp, 255, endpoint=True))
        spectral_functions = np.where(spectral_functions <= 1e-10, 1e-10, spectral_functions)

        contour = ax.contourf(self.x, self.y, spectral_functions, cmap=cmap, 
                              locator=ticker.LogLocator(), 
                              levels=levels, extend='both', zorder=-20)
        contour.set_edgecolor('face')
        
        if not self.scale_energy_by_gaps and continuum_color is not None:
            ax.plot(self.x, self.true_gaps, color=continuum_color, ls=":")
        ax.set_rasterization_zorder(-10)
        ax.set_xscale(self.xscale)
        ax.set_yscale(self.yscale)
        
        self.data_dict[f"spectral_functions_{'higgs' if which == 'amplitude_SC' else 'phase'}"] = spectral_functions
        if "continuum_edge" not in self.data_dict:
            self.data_dict["continuum_edge"] = self.true_gaps
        
        return contour

    def plot(self, axes, cmap, labels=True, min_exp=-1.5, max_exp=2.1):
        contour_higgs = self.plot_one(axes[0], cmap, "amplitude_SC", min_exp=min_exp, max_exp=max_exp) 
        self.plot_one(axes[1], cmap, "phase_SC", min_exp=min_exp, max_exp=max_exp)
        
        if labels:
            if self.scale_energy_by_gaps:
                axes[0].set_ylabel("Higgs\n " + legend(r"\omega / (2 \Delta_\mathrm{max})"))
                axes[1].set_ylabel("Phase\n " + legend(r"\omega / (2 \Delta_\mathrm{max})"))
            else:
                axes[0].set_ylabel("Higgs\n " + legend(r"\omega / W"))
                axes[1].set_ylabel("Phase\n " + legend(r"\omega / W"))
        axes[1].set_xlabel(self.xlabel)

        return contour_higgs

def __get_cf_ignore__(cf_ignore, i):
    if isinstance(cf_ignore, list):
        return __get_cf_ignore__(cf_ignore[i], 0)
    if isinstance(cf_ignore, tuple):
        return BaseCFIgnore(*cf_ignore)
    return cf_ignore

def create_plot(tasks, xscale="linear", scale_energy_by_gaps=False, 
                cmap='inferno', 
                energy_range=None, 
                fig=None, axes=None, 
                cf_ignore=BaseCFIgnore(),
                figure_generator=None,
                height_to_width_ratio=0.55,
                special_labels=[],
                min_exp=-1.5,
                max_exp=2.1):
    if energy_range is None:
        energy_range = (1e-10, 0.29) if not scale_energy_by_gaps else (1e-10, 1.95)
    if fig is None:
        assert(axes is None)
        __figkwargs = {"nrows": 2, "ncols": len(tasks), "sharex": True, "sharey": True, "height_to_width_ratio": height_to_width_ratio}
        if figure_generator is None:
            figure_generator = create_large_figure if len(tasks) > 2 else create_normal_figure
        fig, axes = figure_generator(**__figkwargs)
        fig.subplots_adjust(wspace=0.05, hspace=0.05)
    
    plotters = []
    if len(tasks) > 1:
        label_axes(axes, color="white", weight="bold", special_labels=special_labels)

        for i, (data_query, x_column, xlabel) in enumerate(tasks):
            plotters.append(HeatmapPlotter(data_query, x_column, xlabel=xlabel, xscale=xscale, 
                                           energy_range=energy_range, scale_energy_by_gaps=scale_energy_by_gaps, 
                                           cf_ignore=__get_cf_ignore__(cf_ignore, i)))
            contour_for_colorbar = plotters[-1].plot(axes[:, i], labels=not bool(i), cmap=cmap, min_exp=min_exp, max_exp=max_exp)
    else:
        for i, ax in enumerate(axes):
            ax.annotate(
                f"({string.ascii_lowercase[i]})",
                xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize', 
                verticalalignment='top', color="white", weight="bold")
        
        for i, (data_query, x_column, xlabel) in enumerate(tasks):
            plotters.append(HeatmapPlotter(data_query, x_column, xlabel=xlabel, xscale=xscale, 
                                           energy_range=energy_range, scale_energy_by_gaps=scale_energy_by_gaps, 
                                           cf_ignore=__get_cf_ignore__(cf_ignore, i)))
            contour_for_colorbar = plotters[-1].plot(axes[:], labels=not bool(i), cmap=cmap, min_exp=min_exp, max_exp=max_exp)

    cbar = fig.colorbar(contour_for_colorbar, ax=axes.ravel().tolist(), 
                        orientation='vertical', fraction=0.1, pad=0.025, extend='max')
    
    cbar.locator = ticker.LogLocator(10)
    cbar.set_label(legend(r'\mathcal{A}(\omega) / W^{-1}'))
    
    if hasattr(axes[0], "__len__"):
        make_panels_touch(fig, axes)
    else:
        make_panels_touch(fig, axes, touch_x=True, touch_y=False)
    
    return fig, axes, plotters, cbar


def create_reduced_plot(tasks, xscale="linear", 
                        which='Higgs',
                        scale_energy_by_gaps=False, 
                        cmap='inferno', 
                        energy_range=None, 
                        fig=None, axes=None, 
                        cf_ignore=BaseCFIgnore(),
                        figure_generator=None,
                        special_labels=[],
                        min_exp=-1.5,
                        max_exp=2.1):
    if energy_range is None:
        energy_range = (1e-10, 0.29) if not scale_energy_by_gaps else (1e-10, 1.95)
    if fig is None:
        assert(axes is None)
        __figkwargs = {"ncols": len(tasks), "sharex": True, "sharey": True, "height_to_width_ratio": 0.3}
        if figure_generator is None:
            figure_generator = create_large_figure# if len(tasks) > 2 else create_normal_figure
        fig, axes_raw = figure_generator(**__figkwargs)
        axes = np.asarray(axes_raw)
        fig.subplots_adjust(wspace=0.05, hspace=0.05)
    
    if which == "Higgs":
        which_internal = "amplitude_SC"
    elif which == "Phase":
        which_internal = "phase_SC"
    else:
        raise ValueError(f"which = {which} --- not recognized. Choose 'Higgs' or 'Phase'.")
    
    plotters = []

    if len(axes) > 1:
        label_axes(axes, color="white", weight="bold", special_labels=special_labels)
        
    for i, (data_query, x_column, xlabel) in enumerate(tasks):
        plotters.append(HeatmapPlotter(data_query, x_column, xlabel=xlabel, xscale=xscale, 
                                       energy_range=energy_range, scale_energy_by_gaps=scale_energy_by_gaps, 
                                       cf_ignore=__get_cf_ignore__(cf_ignore, i)))
        contour_for_colorbar = plotters[-1].plot_one(axes[i], cmap=cmap, which=which_internal, min_exp=min_exp, max_exp=max_exp)
        axes[i].set_xlabel(xlabel)
        
    if scale_energy_by_gaps:
        axes[0].set_ylabel(legend(r"\omega / (2 \Delta_\mathrm{max})"))
    else:
        axes[0].set_ylabel(legend(r"\omega / W"))
    
    cbar = fig.colorbar(contour_for_colorbar, ax=axes.ravel(), 
                        orientation='vertical', aspect=15, pad=0.025, extend='max')
    
    cbar.locator = ticker.LogLocator(10)
    cbar.set_label(legend(r'\mathcal{A}(\omega) / W^{-1}'))
    
    make_panels_touch(fig, axes)
    
    return fig, axes, plotters, cbar