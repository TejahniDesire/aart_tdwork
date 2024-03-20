import subprocess

import kgeo
from matplotlib import ticker


import EZPaths
import os

import image_tools
from aart_func import *
from image_tools import curve_params
from params import *
import importlib
import params
import astroModels
import fileloading
from movieMakerV2 import movieMakerIntensity
from astropy import units as u


def radiiThickThin(ax, ax1, xaxis, mean_radii_Thin, mean_radii_Thick,
                   poi, conv_1_style, r_outer_style,flux_peak_style, action):
    '''
    poi = {
        r_outer,
        flux_peak_thick,
        flux_peak_thin,
        conv_1,
    }
    conv_1_style = {
        "color": 'k',
        "linestyle": "--",
        "linewidth": 3
    }
    r_outer_style = {
        "color": 'dimgrey',
        "linestyle": "-",
        "linewidth": 5
    }
    '''
    # ax.axhline(r_inner, color='k', linewidth=3, linestyle=":")  # , label='Blackhole Inner Shadow'

    ax.axvline(230, color='k', linestyle=":")
    ax.axvline(poi["flux_peak_thin"], color=flux_peak_style["color"],
               linestyle=flux_peak_style["linestyle"], linewidth=flux_peak_style["linewidth"])
    ax.axhline(poi["r_outer"],color=r_outer_style["color"],
               linestyle=r_outer_style["linestyle"], linewidth=r_outer_style["linewidth"])  # , label='Blackhole Outer Shadow'

    ax.plot(xaxis, mean_radii_Thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
    ax.plot(xaxis, mean_radii_Thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
    ax.plot(xaxis, mean_radii_Thin[:, 2], '-.', label='n=2', color='tab:blue', linewidth=3)


    # Labels

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_ylabel("Ring Radii ({})".format(R'$R_g$'))

    ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax.xaxis.get_minorticklabels()) if i % n != 0]

    ax.legend(frameon=False)
    ax.set_xlim(xaxis[0], xaxis[xaxis.size - 1])

    new_ticks = [xaxis[0], 230, xaxis[xaxis.size - 1]]
    ax.set_xticks(new_ticks)
    ax.tick_params('x', length=20, width=1, which='major', labelrotation=90)
    ax.title.set_text('Optically Thin Assumption')
    # Markers

    # Optically Thick

    # ax1.axhline(r_inner, color='k', linewidth=2, linestyle=":")  # , label='Blackhole Inner Shadow'

    # ax1.scatter(xaxis[0],r_outer, marker="o", linewidth=10)
    # ax1.scatter(xaxis[0],r_inner, marker="o", linewidth=10)
    ax1.axvline(230, color='k', linestyle=":")
    ax1.axvline(poi["conv_1"],
                color=conv_1_style["color"], linestyle=conv_1_style["linestyle"], linewidth=conv_1_style["linewidth"])
    ax1.axvline(poi["flux_peak_thick"], color=flux_peak_style["color"],
                linestyle=flux_peak_style["linestyle"], linewidth=flux_peak_style["linewidth"])
    ax1.axhline(poi["r_outer"],color=r_outer_style["color"],
                linestyle=r_outer_style["linestyle"], linewidth=r_outer_style["linewidth"])  # , label='Blackhole Outer Shadow'

    ax1.plot(xaxis, mean_radii_Thick[:, 0], '-', label=R'n=0', color='tab:red', linewidth=3)
    ax1.plot(xaxis, mean_radii_Thick[:, 1], ':', label=R'n=1', color='tab:orange', linewidth=3)
    ax1.plot(xaxis, mean_radii_Thick[:, 2], '--', label=R'n=2', color='tab:blue', linewidth=3)
    ax1.plot(xaxis, mean_radii_Thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

    # Labels
    ax1.set_xlabel(astroModels.var_label[action["var"]].replace('=', '')
                   + ' (' + astroModels.units_label[action["var"]] + ')')
    ax1.set_ylabel("Ring Radii ({})".format(R'$R_g$'))
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax1.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))

    new_ticks = [xaxis[0], 230, poi["conv_1"], poi["flux_peak_thick"], xaxis[xaxis.size - 1]]
    ax1.set_xticks(new_ticks)

    # new_ticks = np.append(ax1.get_yticks(), r_outer)
    # new_ticks = np.append(new_ticks, r_inner)
    # print(new_ticks)
    # ax1.set_yticks(new_ticks)
    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax1.xaxis.get_minorticklabels()) if i % n != 0]
    ax1.legend(frameon=False)
    ax1.set_xlim(xaxis[0], xaxis[xaxis.size - 1])

    ax1.tick_params('x', length=20, width=1, which='major', labelrotation=80)
    ax1.title.set_text('Full Solution')


def fluxThickThin(ax, ax1, xaxis, janksys_thin, janksys_thick,
                  poi, conv_1_style, r_outer_style,flux_peak_style, action):
    ax.plot(xaxis, janksys_thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
    ax.plot(xaxis, janksys_thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
    ax.plot(xaxis, janksys_thin[:, 2], '--', label='n=2', color='tab:blue', linewidth=3)
    ax.plot(xaxis, janksys_thin[:, 3], '-.', label='Total', color='tab:purple', linewidth=3)

    # TODO SHOULD I MARK THE PEAK?

    ax.axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
    ax.axvline(230, color='k', linestyle=":")
    ax.axvline(poi["flux_peak_thin"], color=flux_peak_style["color"],
               linestyle=flux_peak_style["linestyle"], linewidth=flux_peak_style["linewidth"])

    # Labels
    ax.set_ylabel("Total Flux ({})".format(R'$J_y$'))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    # ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    # ax1.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.4f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))

    ax.tick_params('x', which="both", labelbottom=False)
    ax.title.set_text('Optically Thin Assumption')

    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax.xaxis.get_minorticklabels()) if i % n != 0]
    ax.tick_params('both', length=10, width=1, which='major')
    ax.set_xlim(xaxis[0], xaxis[xaxis.size - 1])
    # ax.legend(loc='lower left')

    # Optically Thick

    ax1.axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
    ax1.axvline(230, color='k', linestyle=":")
    ax1.axvline(poi["conv_1"],
                color=conv_1_style["color"], linestyle=conv_1_style["linestyle"], linewidth=conv_1_style["linewidth"])
    ax1.axvline(poi["flux_peak_thick"], color=flux_peak_style["color"],
                linestyle=flux_peak_style["linestyle"], linewidth=flux_peak_style["linewidth"])

    ax1.plot(xaxis, janksys_thick[:, 0], '-', label=R'$n=0$', color='tab:red', linewidth=3)
    ax1.plot(xaxis, janksys_thick[:, 1], ':', label=R'$n=1$', color='tab:orange', linewidth=3)
    ax1.plot(xaxis, janksys_thick[:, 2], '--', label=R'$n=2$', color='tab:blue', linewidth=3)
    ax1.plot(xaxis, janksys_thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

    # Labels
    ax1.set_ylabel("Total Flux ({})".format(R'$J_y$'))
    ax1.set_xlabel(
        astroModels.var_label[action["var"]].replace('=', '') + ' (' + astroModels.units_label[action["var"]] + ')')
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    # ax1.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.4f'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))
    ax1.title.set_text('Full Solution')

    new_ticks = [xaxis[0], 230, poi["conv_1"], poi["flux_peak_thick"], xaxis[xaxis.size - 1]]
    ax1.set_xticks(new_ticks)

    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax1.xaxis.get_minorticklabels()) if i % n != 0]
    ax1.tick_params('both', length=10, width=1, which='major')
    ax1.set_xlim(xaxis[0], xaxis[xaxis.size - 1])

    ax1.legend(loc='lower left')


def opticalDepth(ax,xaxis,mean_optical_depth,
                 poi, conv_1_style,flux_peak_style, action):

    ax.axvline(230, color='k', linestyle=":")
    ax.axvline(230, color='k', linestyle=":")
    ax.axvline(poi["conv_1"],
               color=conv_1_style["color"], linestyle=conv_1_style["linestyle"], linewidth=conv_1_style["linewidth"])
    ax.axvline(poi["flux_peak_thick"], color=flux_peak_style["color"],
               linestyle=flux_peak_style["linestyle"], linewidth=flux_peak_style["linewidth"])

    ax.plot(xaxis, mean_optical_depth[0], '-', label='n=0', color='tab:red', linewidth=3)
    ax.plot(xaxis, mean_optical_depth[1], ':', label='n=1', color='tab:orange', linewidth=3)
    ax.plot(xaxis, mean_optical_depth[2], '-.', label='n=2', color='tab:blue', linewidth=3)


    ax.set_xscale('log')
    ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))

    new_ticks = [xaxis[0], 230, poi["conv_1"], poi["flux_peak_thick"], xaxis[xaxis.size - 1]]
    ax.set_xticks(new_ticks)

    ax.set_xlabel(astroModels.var_label[action["var"]].replace('=', '')
                  + ' (' + astroModels.units_label[action["var"]] + ')')
    n = 2
    [l.set_visible(False) for (i, l) in enumerate(ax.xaxis.get_minorticklabels()) if i % n != 0]
    ax.set_ylabel("Optical Depth")
    ax.set_xlim(xaxis[0], xaxis[xaxis.size - 1])

    ax.legend()


def histogram(ax,data,xlabel,ylabel):
    ax.hist(data)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def bar(ax,xdata,ydata,xlabel,ylabel,xbarlabel):
    ax.bar(xdata,ydata)
    plt.xticks(xdata, labels=xbarlabel)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def surfacePlot(X,Y,Z,ax,xlabel,ylabel,father_model,father_value):
    ax.plot_surface(X, Y, Z)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.view_init(20, 80)
    ax.title.set_text("All Models of a = .3, " + father_model + "=" + str(father_value))
    # ax.set_xlim([np.min(X) * 1 / 10, np.max(X) * 10])
