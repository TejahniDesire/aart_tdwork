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
from graphValues import *


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
    # ax.set_ylim(3,8)

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
    # ax1.set_ylim(3,8)
    # ax1.set_ylim(3, 8)
    # if not blurr_policy:
    #     ax1.set_ylim(0, 8)
    # else:
    #     ax1.set_ylim(3, 8)


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
    # ax.set_ylim(10e-6, 10e2)
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
    # ax1.set_ylim(10e-6, 10e2)
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


def IntensityVSRadii(fig,ax0,ax1,limit,thin_intensity, thick_intensity,rmax,blurr_policy=False):
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
        images = [thick_intensity[3],thin_intensity[3]]
    else:
        images = [thick_intensity[0],thin_intensity[0]]

    if not blurr_policy:
        peak012, interp012 = image_tools.radii_of_thetaV2_data(thin_intensity[3],params.dx0,rmax=rmax)
        # peak0, interp0 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak1, interp1 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak2, interp2 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        peakAbsorb, interpAbsorb = image_tools.radii_of_thetaV2_data(thick_intensity[3],params.dx0,rmax=rmax)
        vmax = [np.nanmax(thin_intensity[3])*1.2,np.nanmax(thick_intensity[3])*1.2]
    else:
        peak012, interp012 = image_tools.radii_of_thetaV2_data(thin_intensity[0],params.dx0,rmax=rmax)
        # peak0, interp0 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak1, interp1 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        # peak2, interp2 = image_tools.radii_of_thetaV2_data(thin_intensity[3])
        peakAbsorb, interpAbsorb = image_tools.radii_of_thetaV2_data(thick_intensity[0],params.dx0,rmax=rmax)
        vmax = [np.nanmax(thick_intensity[0])*1.2,np.nanmax(thin_intensity[0])*1.2]
    peaks = [peakAbsorb,peak012]
    interps = [interpAbsorb,interp012]

    model = ["for Full Solution","for Thin Assumption"]
    for J in range(1):
        # if J == 0:
        #     axes_0[J].get_xaxis().set_ticks([])
        #     axes_1[J].get_xaxis().set_ticks([])
        x = np.linspace(0, rmax - 1, rsize) * params.dx0
        ptheta = [0, np.pi/2, np.pi]
        colors = [RingStyle["color"][0],RingStyle["color"][1],RingStyle["color"][2],RingStyle["color"][3]]
        parg = []
        for L in range(len(ptheta)):
            parg += [image_tools.rad_to_arg(ptheta[L])]
            axes_0[J].plot(x, interps[J][parg[L]], linewidth=2, color=colors[L],
                           label=R"$\theta= $" + f"{ptheta[L]:.2f}")
            axes_0[J].axvline(peaks[J][parg[L]], color=colors[L])

        # axes_0[J].set_xlim([2, 6])
        axes_0[J].legend()
        axes_0[J].set_xlabel(R"$R_g$")
        axes_0[J].set_ylabel(R"Flux Value " + model[J])
        axes_0[J].set_xlim(x[0],x[len(x)-1])
        axes_0[J].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
        axes_0[J].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

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


def radiiVSVarphi(fig,ax0,ax1,limit,thin_intensity,blurr_policy=False,plot_intensites=False,average=True):
    if plot_intensites:
        peak0, theta0, intent_at_peaks0 = image_tools.radii_of_thetaV2(thin_intensity[0],
                                                                       give_intensities=True,average=average)
        peak1, theta1, intent_at_peaks1 = image_tools.radii_of_thetaV2(thin_intensity[1],
                                                                       give_intensities=True,average=average)
        peak2, theta2, intent_at_peaks2 = image_tools.radii_of_thetaV2(thin_intensity[2],
                                                                       give_intensities=True,average=average)
        peak3, theta3, intent_at_peaks3 = image_tools.radii_of_thetaV2(thin_intensity[3],
                                                                       give_intensities=True,average=average)
        intent_at_peaks = [intent_at_peaks0, intent_at_peaks1, intent_at_peaks2, intent_at_peaks3]
    else:
        if not blurr_policy:
            peak0, theta0 = image_tools.radii_of_thetaV2(thin_intensity[0],average=average)
            peak1, theta1 = image_tools.radii_of_thetaV2(thin_intensity[1],average=average)
            peak2, theta2 = image_tools.radii_of_thetaV2(thin_intensity[2],average=average)
            peak3, theta3 = image_tools.radii_of_thetaV2(thin_intensity[3],average=average)
        else:
            peak0, theta0 = image_tools.radii_of_thetaV2(thin_intensity[0],average=average)


    if blurr_policy:
        peaks = [peak0]

        thetas = [theta0]
    else:
        peaks = [peak0,peak1,peak2,peak3]

        thetas = [theta0,theta1,theta2,theta3]

    colors = [RingStyle["color"][0],RingStyle["color"][1],RingStyle["color"][2],RingStyle["color"][3]]
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
    ax0.set_ylim(2.5,np.max(peaks[i])*1.1)
    ax0.set_xlim(0,2*np.pi)

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



def DONOTUSE_IntensityVSRadiiType2(fig,ax0,ax1,limit,thin_intensity,rmax,blurr_policy=False):
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

    ring_colors = [RingStyle["color"][0], RingStyle["color"][1], RingStyle["color"][2], RingStyle["color"][3]]

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
