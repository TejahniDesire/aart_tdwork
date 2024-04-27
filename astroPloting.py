import subprocess

import kgeo
import matplotlib.colors
from matplotlib import ticker
from matplotlib.lines import Line2D

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
                   poi, conv_1_style, r_outer_style,flux_peak_style, action,blurr_policy=False):
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

    if blurr_policy:
        ax.plot(xaxis, mean_radii_Thin[:, 0], '-', label='Total Radii', color='tab:red', linewidth=3)
    else:
        ax.plot(xaxis, mean_radii_Thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
        ax.plot(xaxis, mean_radii_Thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
        ax.plot(xaxis, mean_radii_Thin[:, 2], '--', label='n=2', color='tab:blue', linewidth=3)
        ax.plot(xaxis, mean_radii_Thin[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)


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
    ax.set_ylim(3,8)

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
    if not blurr_policy:
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
    ax1.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

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
    if not blurr_policy:
        ax1.set_ylim(0, 6)
    else:
        ax1.set_ylim(3, 8)


    ax1.tick_params('x', length=20, width=1, which='major', labelrotation=80)
    ax1.title.set_text('Full Solution')


def fluxThickThin(ax, ax1, xaxis, janksys_thin, janksys_thick,
                  poi, conv_1_style, r_outer_style,flux_peak_style,action,blurr_policy=False):
    ax.plot(xaxis, janksys_thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
    janksky_line_style ={
        "linestyle": ['-',':','--','-.'],
        "label": ['n=0','n=1','n=2','total'],
        "color": ['tab:red','tab:orange','tab:blue','tab:purple'],
        "linewidth": [3,3,3,3]
    }
    if blurr_policy:
        amount_to_plot = 1
    else:
        amount_to_plot = 4

    for i in range(amount_to_plot):
        ax.plot(xaxis, janksys_thin[:, i], janksky_line_style["linestyle"][i],
                label=janksky_line_style["label"][i], color=janksky_line_style["color"][i],
                linewidth=janksky_line_style["linewidth"][i])

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
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))

    ax.tick_params('x', which="both", labelbottom=False)
    ax.title.set_text('Optically Thin Assumption')

    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax.xaxis.get_minorticklabels()) if i % n != 0]
    ax.tick_params('both', length=10, width=1, which='major')
    ax.set_xlim(xaxis[0], xaxis[xaxis.size - 1])
    ax.set_ylim(10e-6, 10e2)
    # ax.legend(loc='lower left')

    # Optically Thick

    ax1.axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
    ax1.axvline(230, color='k', linestyle=":")
    ax1.axvline(poi["conv_1"],
                color=conv_1_style["color"], linestyle=conv_1_style["linestyle"], linewidth=conv_1_style["linewidth"])
    ax1.axvline(poi["flux_peak_thick"], color=flux_peak_style["color"],
                linestyle=flux_peak_style["linestyle"], linewidth=flux_peak_style["linewidth"])

    for i in range(amount_to_plot):
        ax1.plot(xaxis, janksys_thick[:, i], janksky_line_style["linestyle"][i],
                label=janksky_line_style["label"][i], color=janksky_line_style["color"][i],
                linewidth=janksky_line_style["linewidth"][i])

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
    ax1.set_ylim(10e-6, 10e2)
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



def IntensityVSRadiiType1(fig,ax0,ax1,limit,thin_intensity, thick_intensity,rmax,blurr_policy=False):
    """

    Args:
        thin_intensity: [I0,I1,I2,I0 + I1 + I2]
        thick_intensity:
        rmax:

    Returns:

    """
    rsize = image_tools.num_of_radial_points
    # rmax = I0.shape[0] * .4

    axes_0 = [ax0]
    axes_1 = [ax1]
    if not blurr_policy:
        images = [thin_intensity[3], thick_intensity[3]]
    else:
        images = [thin_intensity[0], thick_intensity[0]]

    if not blurr_policy:
        peak012, interp012 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak0, interp0 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak1, interp1 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak2, interp2 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        peakAbsorb, interpAbsorb = image_tools.radii_of_thetaV2_data(thick_intensity[3])
        vmax = [np.nanmax(thin_intensity[3])*1.2,np.nanmax(thick_intensity[3])*1.2]
    else:
        peak012, interp012 = image_tools.radii_of_thetaV2_data(thin_intensity[0])
        # peak0, interp0 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak1, interp1 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak2, interp2 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        peakAbsorb, interpAbsorb = image_tools.radii_of_thetaV2_data(thick_intensity[0])
        vmax = [np.nanmax(thin_intensity[3])*1.2,np.nanmax(thick_intensity[0])*1.2]
    peaks = [peak012, peakAbsorb]
    interps = [interp012, interpAbsorb]

    model = ["for Thin Assumption", "for Full Solution"]
    for J in range(1):
        if J == 0:
            axes_0[J].get_xaxis().set_ticks([])
            axes_1[J].get_xaxis().set_ticks([])
        x = np.linspace(0, rmax - 1, rsize) * params.dx0
        ptheta = [0, np.pi/2, np.pi]  # for frame 33
        colors = ['tab:blue', 'tab:green', 'tab:red']
        parg = []
        for L in range(len(ptheta)):
            parg += [image_tools.rad_to_arg(ptheta[L])]
            axes_0[J].plot(x, interps[J][parg[L]], linewidth=2, color=colors[L],
                           label=R"$\theta= $" + f"{ptheta[L]:.2f}")
            axes_0[J].axvline(peaks[J][parg[L]], color=colors[L])

        axes_0[J].set_xlim([2, 6])
        axes_0[J].legend()
        axes_0[J].set_xlabel(R"$R_g$")
        axes_0[J].set_ylabel(R"Flux Value " + model[J])

        im1 = axes_1[J].imshow(images[J], origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])

        axes_1[J].set_xlim(-10, 10)  # units of M
        axes_1[J].set_ylim(-10, 10)

        axes_1[J].set_xlabel(r"$\alpha$" + " " + r"($M$)")
        axes_1[J].set_ylabel(r"$\beta$" + " " + r"($M$)")

        # Plot lines
        rline1 = np.array([0, 10])
        theta1 = np.array([0, 0])

        alpha1 = rline1 * np.cos(theta1)
        beta1 = rline1 * np.sin(theta1)

        for L in range(len(ptheta)):
            rline = np.array([0, 10])
            theta = np.array([ptheta[L], ptheta[L]])

            alpha = rline * np.cos(theta)
            beta = rline * np.sin(theta)

            axes_1[J].plot(alpha, beta, color=colors[L], linestyle='--')

        colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
            vmax[J] * .8,
            vmax[J] * .6,
            vmax[J] * .4,
            vmax[J] * .2,
            vmax[J] * .05
        ],
                                 label="Brightnes Temperature (K)",
                                 ax=axes_1[J]
                                 )


def IntensityVSRadiiType2(fig,ax0,ax1,limit,thin_intensity,rmax,blurr_policy=False):
    """

    Args:
        thin_intensity: [I0,I1,I2,I0 + I1 + I2]
        thick_intensity:
        rmax:

    Returns:

    """
    rsize = image_tools.num_of_radial_points

    peak012, interp012 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
    peak0, interp0 = image_tools.radii_of_thetaV2_data(thin_intensity[0])
    peak1, interp1 = image_tools.radii_of_thetaV2_data(thin_intensity[1])
    peak2, interp2 = image_tools.radii_of_thetaV2_data(thin_intensity[2])
    # peakAbsorb, interpAbsorb = image_tools.radii_of_thetaV2_data(thick_intensity[3])

    peaks = [peak0,peak1,peak2,peak012]
    interps = [interp0,interp1,interp2,interp012]



    # ptheta = [0, np.pi / 2, np.pi]
    ptheta = 4.29403619  # for frame 33
    colors = ['tab:blue', 'tab:green', 'tab:red','tab:purple']
    ring_colors = ['tab:blue', 'tab:green', 'tab:red','tab:purple']

    model = ["for Thin Assumption", "for Full Solution"]
    nrings = ['0','1','2','Cumulative']
    for J in range(len(peaks)):
        # if J == 0:
        #     axes_0[J].get_xaxis().set_ticks([])
        #     axes_1[J].get_xaxis().set_ticks([])
        x = np.linspace(0, rmax - 1, rsize) * params.dx0
        parg = image_tools.rad_to_arg(ptheta)
        ax0.plot(x, interps[J][parg], linewidth=2, color=ring_colors[J],
                 label=R"$n= $" + nrings[J] + R", $\varphi = $" + f"{ptheta:.2f}")
        ax0.axvline(peaks[J][parg], color=ring_colors[J])

    ax0.set_xlim([2, 6])
    ax0.legend()
    ax0.set_xlabel(R"$R_g$")
    ax0.set_ylabel(R"Flux Value")
    if not blurr_policy:
        vmax = np.nanmax(thin_intensity[3]) * 1.2
        im1 = ax1.imshow(thin_intensity[3], origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])
    else:
        vmax = np.nanmax(thin_intensity[0]) * 1.2
        im1 = ax1.imshow(thin_intensity[0], origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])

    ax1.set_xlim(-10, 10)  # units of M
    ax1.set_ylim(-10, 10)

    ax1.set_xlabel(r"$\alpha$" + " " + r"($M$)")
    ax1.set_ylabel(r"$\beta$" + " " + r"($M$)")

    # Plot lines
    rline1 = np.array([0, 10])
    theta1 = np.array([0, 0])

    alpha1 = rline1 * np.cos(theta1)
    beta1 = rline1 * np.sin(theta1)

    rline = np.array([0, 10])
    theta = np.array([ptheta, ptheta])

    alpha = rline * np.cos(theta)
    beta = rline * np.sin(theta)

    ax1.plot(alpha, beta, color='red', linestyle='--')

    colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax * .8,
        vmax * .6,
        vmax * .4,
        vmax * .2,
        vmax * .05
    ],
                             label="Brightnes Temperature (K)",
                             ax=ax1
                             )


def radiiVSVarphi(fig,ax0,ax1,limit,thin_intensity,blurr_policy=False,plot_intensites=False):
    if plot_intensites:
        peak0, theta0, intent_at_peaks0 = image_tools.radii_of_thetaV2(thin_intensity[0],give_intensities=True)
        peak1, theta1, intent_at_peaks1 = image_tools.radii_of_thetaV2(thin_intensity[1],give_intensities=True)
        peak2, theta2, intent_at_peaks2 = image_tools.radii_of_thetaV2(thin_intensity[2],give_intensities=True)
        peak3, theta3, intent_at_peaks3 = image_tools.radii_of_thetaV2(thin_intensity[3],give_intensities=True)
        intent_at_peaks = [intent_at_peaks0, intent_at_peaks1, intent_at_peaks2, intent_at_peaks3]
    else:
        peak0, theta0 = image_tools.radii_of_thetaV2(thin_intensity[0])
        peak1, theta1 = image_tools.radii_of_thetaV2(thin_intensity[1])
        peak2, theta2 = image_tools.radii_of_thetaV2(thin_intensity[2])
        peak3, theta3 = image_tools.radii_of_thetaV2(thin_intensity[3])

    peaks = [peak0,peak1,peak2,peak3]

    thetas = [theta0,theta1,theta2,theta3]
    colors = ['tab:red','tab:orange','tab:blue','tab:purple']
    labels = [R"$n= 0$",R"$n= 1$",R"$n= 2$",R"$Cumulative$"]
    linewidths = [3,2,1,4]
    linestyles = ['-',':','--','-.']
    alphas = []
    betas = []
    for i in range(len(peaks)):
        # if J == 0:
        #     axes_0[J].get_xaxis().set_ticks([])
        #     axes_1[J].get_xaxis().set_ticks([])
        if not plot_intensites:
            ax0.plot(thetas[i], peaks[i], linewidth=2, color=colors[i],label=labels[i],linestyle=linestyles[i])
            alphas += [peaks[i] * np.cos(thetas[i])]
            betas += [peaks[i] * np.sin(thetas[i])]
            ax0.set_ylabel(R"$R_g$")
            ax1.set_ylim(-10, 10)
        else:
            ax0.plot(thetas[i], intent_at_peaks[i], linewidth=2, color=colors[i],label=labels[i],linestyle=linestyles[i])
            alphas += [intent_at_peaks[i] * np.cos(thetas[i])]
            betas += [intent_at_peaks[i] * np.sin(thetas[i])]
            ax0.set_ylabel(R"Brightness $(K)$")

    ax0.set_xlabel(R"$\varphi$")

    ax0.legend()

    if not blurr_policy:
        vmax = np.nanmax(thin_intensity[3]) * 1.2
        im1 = ax1.imshow(thin_intensity[3], origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])
    else:
        vmax = np.nanmax(thin_intensity[0]) * 1.2
        im1 = ax1.imshow(thin_intensity[0], origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])

    ax1.set_xlim(-10, 10)  # units of M

    ax1.set_xlabel(r"$\alpha$" + " " + r"($M$)")
    ax1.set_ylabel(r"$\beta$" + " " + r"($M$)")

    # Radii Calc
    for i in range(len(peaks)):
        ax1.plot(alphas[i], betas[i],linewidth=linewidths[i], color=colors[i], linestyle=linestyles[i])

    colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax * .8,
        vmax * .6,
        vmax * .4,
        vmax * .2,
        vmax * .05
    ],
                             label="Brightnes Temperature (K)",
                             ax=ax1
                             )



def fullImage(fig,ax0,ax1,limit,thin_intensity,thick_intensity,thin_radii,thick_radii,theta,blurr_policy=False):
    thin_alpha0 = thin_radii[0] * np.cos(theta)
    thin_beta0 = thin_radii[0] * np.sin(theta)

    if not blurr_policy:
        thin_alpha1 = thin_radii[1] * np.cos(theta)
        thin_beta1 = thin_radii[1] * np.sin(theta)
        thin_alpha2 = thin_radii[2] * np.cos(theta)
        thin_beta2 = thin_radii[2] * np.sin(theta)
        thin_alpha_full = thin_radii[3] * np.cos(theta)
        thin_beta_full = thin_radii[3] * np.sin(theta)

    # full solution radii
    thick_alpha0 = thick_radii[0] * np.cos(theta)
    thick_beta0 = thick_radii[0] * np.sin(theta)
    if not blurr_policy:
        thick_alpha1 = thick_radii[1] * np.cos(theta)
        thick_beta1 = thick_radii[1] * np.sin(theta)
        thick_alpha2 = thick_radii[2] * np.cos(theta)
        thick_beta2 = thick_radii[2] * np.sin(theta)
        thick_alpha_full = thick_radii[3] * np.cos(theta)
        thick_beta_full = thick_radii[3] * np.sin(theta)



    # vmax = 10e11
    # vmin = 10e8
    # Optically Thin
    if not blurr_policy:
        # im0im = np.log10(thin_intensity[3])
        # im0im[im0im == -np.Infinity] = 0
        vmax0 = np.nanmax(thin_intensity[3])
        im0 = ax0.imshow(thin_intensity[3], vmax=vmax0, origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])
    else:
        # im0im = np.log10(thin_intensity[0])
        # im0im[im0im == -np.Infinity] = 0
        vmax0 = np.nanmax(thin_intensity[0])
        im0 = ax0.imshow(thin_intensity[0], vmax=vmax0, origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])


    # im0 = ax0.imshow(im0im, origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit],
    #                  norm=matplotlib.colors.LogNorm(vmin,vmax))



    ax0.set_xlim(-10, 10)  # units of M
    ax0.set_ylim(-10, 10)

    ax0.set_xlabel(r"$\alpha$" + " " + r"($M$)")
    ax0.set_ylabel(r"$\beta$" + " " + r"($M$)")
    ax0.title.set_text('Optically Thin Assumption')

    # Optically thick
    if not blurr_policy:
        # im0im = np.log10(thick_intensity[3])
        # im0im[im0im == -np.Infinity] = 0
        vmax1 = np.nanmax(thick_intensity[3])
        im1 = ax1.imshow(thick_intensity[3], vmax=vmax1, origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])
    else:
        # im0im = np.log10(thick_intensity[0])
        # im0im[im0im == -np.Infinity] = 0
        vmax1 = np.nanmax(thick_intensity[0])
        im1 = ax1.imshow(thick_intensity[0], vmax=vmax1, origin="lower", cmap="afmhot", extent=[-limit, limit, -limit, limit])

    #
    ax1.set_xlim(-10, 10)  # units of M
    ax1.set_ylim(-10, 10)

    ax1.set_xlabel(r"$\alpha$" + " " + r"($M$)")

    ax1.title.set_text('Full Solution')
    vmax = vmax0
    colorbar0 = fig.colorbar(im0, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt), ticks=[
        vmax * .8,
        vmax * .6,
        vmax * .4,
        vmax * .2,
        vmax * .05
    ],
                             ax=ax0
                             )
    vmax = vmax1
    colorbar1 = fig.colorbar(im1, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt), ticks=[
        vmax * .8,
        vmax * .6,
        vmax * .4,
        vmax * .2,
        vmax * .05
    ],
                             label=R"$(Brightness Temperature) (1e9 K)$",
                             ax=ax1
                             )
    # colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt),ax=ax0)
    # colorbar1 = fig.colorbar(im1, fraction=0.046, pad=0.04, format=ticker.FuncFormatter(fmt),
    #                          label=R"$Log_{10}(Brightness Temperature) (1e9 K)$",ax=ax1)

    '''Radii Calc______________________'''
    # Thin
    lineCum_thickness = 4
    line0_thickness = 3
    line1_thickness = 2
    line2_thickness = 1

    if not blurr_policy:
        ax0.plot(thin_alpha_full, thin_beta_full, color='tab:purple', linestyle='--', linewidth=lineCum_thickness)
        ax0.plot(thin_alpha1, thin_beta1, color='tab:orange', linestyle=':', linewidth=line1_thickness)
        ax0.plot(thin_alpha2, thin_beta2, color='tab:blue', linestyle='--', linewidth=line2_thickness)

    ax0.plot(thin_alpha0, thin_beta0, color='tab:red', linestyle='-', linewidth=line0_thickness)

    # Thick
    if not blurr_policy:
        ax1.plot(thick_alpha_full, thick_beta_full, color='tab:purple', linestyle='--', linewidth=lineCum_thickness,
                 label='Cumulative')
        ax1.plot(thick_alpha1, thick_beta1, color='tab:orange', linestyle=':', linewidth=line1_thickness, label=R'n=1')
        ax1.plot(thick_alpha2, thick_beta2, color='tab:blue', linestyle='--', linewidth=line2_thickness, label=R'n=2')

    ax1.plot(thick_alpha0, thick_beta0, color='tab:red', linestyle='-', linewidth=line0_thickness, label=R'n=0')
    ax1.legend()

#
# def parameterLaws():
#     fig = plt.subplots(4, 1, sharex='col', figsize=(5, 15), height_ratios=[1, 1, 1, 2])
#     # [r, theta_e.value, n.value, b_field.value, b_nu_fluid.value,
#     #                                   acoeff_I_fluid.value, tau_curve.value, specific_intensity_thin_packed.value,
#     #                                  specific_intensity_thick_packed.value], axis=0)
#     ax = [None, None, None, None]
#
#     # Feducial ORange
#     # 	1.23e4, # n_th0
#     # 	8.1e9, # t_e0
#     # 	-.7, # p_dens
#     # 	-.84 # p_temp
#
#     # Steeper Blue
#     # 	2.9726e+05, # n_th0
#     # 	-.7, # p_dens
#     # 	-1.6 # p_temp
#
#     # Shallower RED
#     # 	2.1526e+04, # n_th0
#     # 	-.7, # p_dens
#     # 	-.3 # p_temp
#
#     i = 0
#     n = 3000
#
#     # Subplot 1-----------------------------------------
#     ax[0] = plt.subplot(4, 1, 1)
#     ax[0].plot(bp0s["bp_shallowT230"][0, :], bp0s["bp_shallowT230"][1, :], 'tab:red')
#     ax[0].plot(bp0s["bp_fiducial230"][0, :], bp0s["bp_fiducial230"][1, :], 'tab:orange')
#     ax[0].plot(bp0s["bp_steeperT230"][0, :], bp0s["bp_steeperT230"][1, :], 'tab:blue')
#
#     ax[0].set_ylabel(R'$\theta_e$', fontsize=18)
#     ax[0].set_yscale('log')
#
#     # Subplot 2-----------------------------------------
#     ax[1] = plt.subplot(4, 1, 2)
#     ax[1].plot(bp0s["bp_shallowT230"][0, :], bp0s["bp_shallowT230"][2, :], 'tab:red')
#     ax[1].plot(bp0s["bp_fiducial230"][0, :], bp0s["bp_fiducial230"][2, :], 'tab:orange')
#     ax[1].plot(bp0s["bp_steeperT230"][0, :], bp0s["bp_steeperT230"][2, :], 'tab:blue')
#
#     ax[1].set_ylabel('Density ({})'.format(R'$cm^{-3}$'), fontsize=18)
#     ax[1].set_yscale('log')
#
#     # Subplot 3-----------------------------------------
#     ax[2] = plt.subplot(4, 1, 3)
#
#     ax[2].plot(bp0s["bp_shallowT230"][0, :], bp0s["bp_shallowT230"][3, :], 'tab:red', linewidth=7)
#     ax[2].plot(bp0s["bp_fiducial230"][0, :], bp0s["bp_fiducial230"][3, :], 'tab:orange', linewidth=5)
#     ax[2].plot(bp0s["bp_steeperT230"][0, :], bp0s["bp_steeperT230"][3, :], 'tab:blue', linewidth=2)
#
#     ax[2].set_yscale('log')
#     ax[2].set_ylabel('B (Gauss)', fontsize=18)
#
#     # Subplot 4-----------------------------------------
#     ax[3] = plt.subplot(4, 1, 4)
#
#     ax[3].plot(bp0s["bp_shallowT230"][0, i::n], bp0s["bp_shallowT230"][7, i::n], 'tab:red', label="Model 3")
#     ax[3].plot(bp0s["bp_fiducial230"][0, i::n], bp0s["bp_fiducial230"][7, i::n], 'tab:orange', label="Model 1")
#     ax[3].plot(bp0s["bp_steeperT230"][0, i::n], bp0s["bp_steeperT230"][7, i::n], 'tab:blue', label="Model 2")
#
#     ax[3].plot(bp0s["bp_shallowT86"][0, i::n], bp0s["bp_shallowT86"][7, i::n], 'tab:red', linewidth=3, label="Model 3",
#                linestyle=(0, (1, 5)))
#     ax[3].plot(bp0s["bp_fiducial86"][0, i::n], bp0s["bp_fiducial86"][7, i::n], 'tab:orange', linewidth=3,
#                label="Model 1", linestyle=(0, (1, 5)))
#     ax[3].plot(bp0s["bp_steeperT86"][0, i::n], bp0s["bp_steeperT86"][7, i::n], 'tab:blue', linewidth=3, label="Model 2",
#                linestyle=(0, (1, 5)))
#
#     ax[3].plot(bp0s["bp_shallowT345"][0, i::n], bp0s["bp_shallowT345"][7, i::n], 'tab:red', label="Model 3",
#                linestyle='-.')
#     ax[3].plot(bp0s["bp_fiducial345"][0, i::n], bp0s["bp_fiducial345"][7, i::n], 'tab:orange', label="Model 1",
#                linestyle='-.')
#     ax[3].plot(bp0s["bp_steeperT345"][0, i::n], bp0s["bp_steeperT345"][7, i::n], 'tab:blue', label="Model 2",
#                linestyle='-.')
#
#     ax[3].set_ylabel(R'$J_\nu$ ($erg \cdot cm^{-3} \cdot Hz^{-1} \cdot s^{-1}$)', fontsize=18)
#     ax[3].set_yscale('log')
#     ax[3].set_xscale('log')
#     ax[3].minorticks_on()
#     ax[3].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
#     ax[3].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
#     ax[3].set_xlabel('Radial Distance ({})'.format(R'$R_g$'), fontsize=18)
#     ax[3].set_xlim([2, 35])
#
#     lines = [
#         Line2D([0], [0], marker='o', markerfacecolor='tab:orange', color='w', markersize=12),
#         Line2D([0], [0], marker='o', markerfacecolor='tab:blue', color='w', markersize=12),
#         Line2D([0], [0], marker='o', markerfacecolor='tab:red', color='w', markersize=12),
#         Line2D([0, 1], [0, 1], linestyle=(0, (1, 5)), color='k'),
#         Line2D([0, 1], [0, 1], linestyle='-', color='k'),
#         Line2D([0, 1], [0, 1], linestyle='-.', color='k')
#     ]
#     labels = [
#         "Model 1",
#         "Model 2",
#         "Model 3",
#         R'$\nu = 86$GHz',
#         R'$\nu = 230$GHz',
#         R'$\nu = 345$GHz'
#     ]
#
#     ax[3].legend(lines, labels)
#     ax[3].set_ylim([1e-31, 1e-18])
#     ax[3].tick_params('x', length=10, width=1, which='major', labelrotation=90)
#
#     plt.savefig(fig_path + "emission_profiles.png", bbox_inches='tight')


def fmt(x, pos):
    x = x / 1e9
    return '{:.2f}'.format(x)


def histogram(ax,data,xlabel,ylabel):
    ax.hist(data)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def bar(ax,xdata,ydata,xlabel,ylabel,xbarlabel):
    ax.bar(xdata,ydata)
    plt.xticks(xdata, labels=xbarlabel, rotation=90)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def surfacePlot(X,Y,Z,ax,xlabel,ylabel,zlabel,father_model,father_value):
    ax.plot_surface(X, Y, Z)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.view_init(20, 80)
    ax.title.set_text("All Models of a = .3, " + father_model + "=" + str(father_value))
    # ax.set_xlim([np.min(X) * 1 / 10, np.max(X) * 10])
