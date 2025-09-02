using Korg
using Makie
using GLMakie
using FITSIO
using Dierckx
using DataFrames
using CSV

function find_near_nodes(x, nodes)
    x_2 = findfirst(x -> x ≥ 0, nodes .- x)
    return x_2 - 1, x_2
end

function linear_interpolation(x, x_1, x_2, y_1, y_2)
    return y_1 + (x - x_1)*(y_2 - y_1)/(x_2 - x_1)
end

function interpolate_kurucz(T_eff, logg)
    logg_nodes = [0.0:0.5:5.0;]
    T_eff_nodes = [3500:250.0:9750; 10000:5e2:12500; 13000:1e3:35000]

    i_T_eff_1, i_T_eff_2 = find_near_nodes(T_eff, T_eff_nodes)
    i_logg_1, i_logg_2 = find_near_nodes(logg, logg_nodes)

    kurucz_fits_1 = FITS("kurucz-p00/kp00_$(round(Int, T_eff_nodes[i_T_eff_1])).fits")
    kurucz_fits_2 = FITS("kurucz-p00/kp00_$(round(Int, T_eff_nodes[i_T_eff_2])).fits")

    λs_1 = read(kurucz_fits_1[2], "WAVELENGTH") 
    λs_2 = read(kurucz_fits_2[2], "WAVELENGTH")


    kurucz_spec_11 = read(kurucz_fits_1[2], "g$(round(Int, logg_nodes[i_logg_1]*10))")
    kurucz_spec_12 = read(kurucz_fits_1[2], "g$(round(Int, logg_nodes[i_logg_2]*10))")

    kurucz_spec_21 = read(kurucz_fits_2[2], "g$(round(Int, logg_nodes[i_logg_1]*10))")
    kurucz_spec_22 = read(kurucz_fits_2[2], "g$(round(Int, logg_nodes[i_logg_2]*10))")

    kurucz_spec_1 = linear_interpolation(logg, logg_nodes[[i_logg_1, i_logg_2]]..., kurucz_spec_11, kurucz_spec_12)
    kurucz_spec_2 = linear_interpolation(logg, logg_nodes[[i_logg_1, i_logg_2]]..., kurucz_spec_21, kurucz_spec_22)

    if !isnothing(findfirst(x -> abs(x) > 1e-3, λs_1 .- λs_2))
        kurucz_spec_2_spl = Spline1D(λs_2, kurucz_spec_2, k = 1)
        kurucz_spec_2 .= kurucz_spec_2_spl.(λs_1)
    end
        
    kurucz_spec = linear_interpolation(T_eff, T_eff_nodes[[i_T_eff_1, i_T_eff_2]]..., kurucz_spec_1, kurucz_spec_2)
    return λs_1, kurucz_spec
end

function plot_kurucz_korg(T_eff, logg)
    λs, kurucz_spec = interpolate_kurucz(T_eff, logg)
    wls, flux, cntm = Korg.synth(Teff = T_eff,
                             logg = logg, 
                             m_H = 0.0, 
                             wavelengths = Korg.Wavelengths(5000, 12000, 1), 
                             linelist = Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac"), 
                             R = 100)

    fig = Figure()
    ax = Axis(fig[1,1])
    xlims!(ax, (5000, 12000)) 
    colormap = Makie.wong_colors()
    lines!(ax, λs, kurucz_spec, color = colormap[1])
    lines!(ax, wls, flux .* cntm * 1e-8, color = colormap[1], linestyle = :dash)
    fig
end

function calc_tess_kurucz_flux_no_extinction(T_eff, logg)
    λs, kurucz_spec = interpolate_kurucz(T_eff, logg)
# wls, flux, cntm = Korg.synth(Teff = T_eff,
#                              logg = logg, 
#                              m_H = 0.0, 
#                              wavelengths = Korg.Wavelengths(5000, 12000, 1), 
#                              linelist = Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac"), 
#                              R = 100)

    tess_response_df = CSV.read("tess-response-function-v2.0.csv", DataFrame, comment = "#")
    tess_response_spl = Spline1D(tess_response_df[:,1] * 10, tess_response_df[:,2], k = 1)

    Δλs_kurucz = λs[2:end] - λs[1:end-1]
    λs_mid_kurucz = (λs[2:end] + λs[1:end-1])/2
    kurucz_spec_mid = (kurucz_spec[2:end] + kurucz_spec[1:end-1])/2
    return sum(@. kurucz_spec_mid * tess_response_spl(λs_mid_kurucz) * Δλs_kurucz)
end

extinct_law(λ, Av) = Av/(2.5*log10(ℯ))/(λ/5e3)

function calc_tess_kurucz_flux_with_extinction(T_eff, logg, Av)
    λs, kurucz_spec = interpolate_kurucz(T_eff, logg)
# wls, flux, cntm = Korg.synth(Teff = T_eff,
#                              logg = logg, 
#                              m_H = 0.0, 
#                              wavelengths = Korg.Wavelengths(5000, 12000, 1), 
#                              linelist = Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac"), 
#                              R = 100)

    tess_response_df = CSV.read("tess-response-function-v2.0.csv", DataFrame, comment = "#")
    tess_response_spl = Spline1D(tess_response_df[:,1] * 10, tess_response_df[:,2], k = 1)

    Δλs_kurucz = λs[2:end] - λs[1:end-1]
    λs_mid_kurucz = (λs[2:end] + λs[1:end-1])/2
    kurucz_spec_mid = (kurucz_spec[2:end] + kurucz_spec[1:end-1])/2
    return sum(@. kurucz_spec_mid * tess_response_spl(λs_mid_kurucz) * extinct_law(λs_mid_kurucz, Av) * Δλs_kurucz)
end

tess_rel_no_Av = calc_tess_kurucz_flux_no_extinction(7800, 4.2)/calc_tess_kurucz_flux_no_extinction(4700, 2.2)*1.5^2/26.6^2
tess_rel_Av = calc_tess_kurucz_flux_with_extinction(7800, 4.2, 1.6)/calc_tess_kurucz_flux_with_extinction(4700, 2.2, 1.6)*1.5^2/26.6^2