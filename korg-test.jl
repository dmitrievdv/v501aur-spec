using Korg
using Makie
using GLMakie

sp_class_temps = Dict(("A7" => 7920),
                ("F0" => 7350),
                ("F5" => 6700),
                ("G0" => 6050),
                ("G5" => 5660),
                ("K0" => 5240),
                ("K5" => 4400),
                ("M0" => 3750),
                ("M5" => 3200))

sp_classes = keys(sp_class_temps)
wls = Korg.Wavelengths(2000,11000,1)
n_wls = length(wls)
fluxes = Dict(zip(sp_classes, [zeros(n_wls) for sp_class in sp_classes]))
continuums = Dict(zip(sp_classes, [zeros(n_wls) for sp_class in sp_classes]))

for sp_class in sp_classes
    println(sp_class)
    A_X = format_A_X()
    atm = interpolate_marcs(sp_class_temps[sp_class], 4.0, A_X; clamp_abundances = true)
    synth_result = synthesize(atm, Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac"), A_X, wls)
    fluxes[sp_class] = Korg.apply_LSF(synth_result.flux, wls, 100)
    continuums[sp_class] = synth_result.cntm
end

fig = Figure()
ax = Axis(fig[1,1], yscale = log10)
# lines!(ax, wls_kurucz, flux_kurucz .* continuum_kurucz, label = "kurucz")
for sp_class in sp_classes
    lines!(ax, wls.all_wls, fluxes[sp_class] .* continuums[sp_class], label = sp_class)
end
axislegend(ax)
fig