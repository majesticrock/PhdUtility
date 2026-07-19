# A set of functions to unify legends in plots
# The idea is, to always use one of the functions provided here instead of manually setting legends
def legend(quantity, unit = None, unit_exp = None):
    if unit is not None:
        __exp = f"^{{{unit_exp}}}" if unit_exp is not None else ""
        return f"${quantity}$ $(\\mathrm{{{unit}}}{__exp})$"
    return f"${quantity}$"

def log_legend(quantity, unit = None, unit_exp = None):
    if unit is not None:
        __exp = f"^{{{unit_exp}}}" if unit_exp is not None else ""
        return f"$\\ln [ {quantity} (\\mathrm{{{unit}}}{__exp}) ]$"
    return f"$\\ln [ {quantity} ]$"